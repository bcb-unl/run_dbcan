from __future__ import annotations

import logging
from pathlib import Path
from abc import ABC
import csv
import psutil
import pyhmmer
import time
from typing import List

from dbcan.configs.pyhmmer_config import (
    PyHMMERConfig,
    PyHMMERDBCANConfig,
    DBCANSUBConfig,
    PyHMMERSTPConfig,
    PyHMMERTFConfig,
    PyHMMERPfamConfig
)
from dbcan.process.process_utils import process_results
from dbcan.process.process_dbcan_sub import DBCANSUBProcessor
from dbcan.utils.memory_monitor import MemoryMonitor, get_memory_monitor
import dbcan.constants.pyhmmer_search_constants as P

logger = logging.getLogger(__name__)


class PyHMMERProcessor(ABC):
    """Base PyHMMER processor: config is the single source of truth."""

    # Subclasses must set these class attributes
    # HMM_FILE: str = ""
    # OUTPUT_FILE: str = ""
    # EVALUE_ATTR: str = ""          # name of e-value attribute in config
    # COVERAGE_ATTR: str = ""        # name of coverage attribute in config
    # USE_NULL_INPUT: bool = False   # for Pfam (optional alternate input)

    def __init__(self, config: PyHMMERConfig):
        self.config = config
        self._validate_basic()

    # -------- Properties --------
    @property
    def hmm_file(self) -> Path:
        return Path(self.config.db_dir) / self.config.hmm_file

    @property
    def input_faa(self) -> Path:
        return Path(self.config.output_dir) / self.config.input_faa

    @property
    def output_file(self) -> Path:
        return Path(self.config.output_dir) / self.config.output_file

    @property
    def e_value_threshold(self) -> float:
        return float(self.config.evalue_threshold)

    @property
    def coverage_threshold(self) -> float:
        return float(self.config.coverage_threshold)

    @property
    def hmmer_cpu(self) -> int:
        return int(self.config.threads)

    # -------- Validation --------
    def _validate_basic(self):
        self.output_file.parent.mkdir(parents=True, exist_ok=True)
        # Existence checks deferred to run() for flexibility

    # -------- Batch writing helper --------
    def _flush_hit_buffer(self, hit_buffer: List, writer: csv.writer) -> None:
        """Flush hit buffer to CSV file.
        
        Args:
            hit_buffer: List of hit rows to write
            writer: CSV writer for results
        """
        if hit_buffer:
            writer.writerows(hit_buffer)
            hit_buffer.clear()

    # -------- Core search --------
    def _process_sequence_block(self, seq_block, hmms: List, cpus: int, hit_buffer: List, 
                                hmm_file_stem: str, buffer_flush_threshold: int = 5000) -> int:
        """Process a sequence block (streaming mode, no file I/O).
        
        Args:
            seq_block: DigitalSequenceBlock from read_block()
            hmms: List of HMM objects (pre-loaded from HMMFile)
            cpus: Number of CPUs to use
            hit_buffer: List to collect hits for batch writing
            hmm_file_stem: Pre-computed HMM file stem (constant)
            buffer_flush_threshold: Number of hits to collect before flushing buffer
        
        Returns:
            int: Number of hits found
        """
        hit_count = 0
        
        try:
            for hits in pyhmmer.hmmsearch(
                hmms,
                seq_block,
                cpus=cpus,
                domE=self.e_value_threshold
            ):
                for hit in hits:
                    for domain in hit.domains.included:
                        aln = domain.alignment
                        coverage = (aln.hmm_to - aln.hmm_from + 1) / aln.hmm_length
                        hmm_name = aln.hmm_name.decode('utf-8')
                        if P.GT2_PREFIX in hmm_name:
                            hmm_name = P.GT2_FAMILY_NAME
                        i_evalue = domain.i_evalue
                        if i_evalue < self.e_value_threshold and coverage > self.coverage_threshold:
                            hit_buffer.append([
                                hmm_name,
                                aln.hmm_length,
                                aln.target_name.decode('utf-8'),
                                aln.target_length,
                                i_evalue,
                                aln.hmm_from,
                                aln.hmm_to,
                                aln.target_from,
                                aln.target_to,
                                coverage,
                                hmm_file_stem
                            ])
                            hit_count += 1
        except Exception as e:
            logger.error(f"Error processing sequence block: {e}")
            raise
        
        return hit_count

    def _iter_hmmsearch_hits(self, hmms_or_hmmfile, targets, cpus: int):
        """Small wrapper so we have exactly one place to tune hmmsearch kwargs."""
        return pyhmmer.hmmsearch(
            hmms_or_hmmfile,
            targets,
            cpus=cpus,
            domE=self.e_value_threshold,
        )
    
    def _calculate_batch_size(self, input_faa: Path, memory_monitor: MemoryMonitor, retry_count: int = 0) -> int:
        """Calculate optimal batch size based on available memory.
        
        Args:
            input_faa: Input FASTA file path
            memory_monitor: Memory monitor instance
            retry_count: Number of retries (used to reduce batch size on retry)
        """
        # Use configured batch size if provided (but reduce on retry)
        if self.config.batch_size is not None and self.config.batch_size > 0:
            batch_size = self.config.batch_size
            if retry_count > 0:
                # Reduce batch size by 50% on each retry
                batch_size = max(100, int(batch_size * (0.5 ** retry_count)))
                logger.info(f"Reduced batch size for retry {retry_count}: {batch_size}")
            else:
                logger.info(f"Using configured batch size: {batch_size}")
            return batch_size
        
        # Estimate batch size based on available memory
        file_size_mb = input_faa.stat().st_size / (1024 * 1024)
        
        # Estimate average sequence size (rough estimate: file_size / estimated_sequence_count)
        # For proteins, average length ~300-400 aa, so roughly 0.01-0.02 MB per sequence
        # We'll use a conservative estimate
        estimated_seq_count = max(1, int(file_size_mb / 0.01))  # Rough estimate
        avg_seq_size_mb = file_size_mb / estimated_seq_count if estimated_seq_count > 0 else 0.01
        
        # Reduce safety factor on retry
        safety_factor = self.config.memory_safety_factor * (0.7 ** retry_count)
        batch_size = memory_monitor.estimate_batch_size(
            avg_seq_size_mb,
            safety_factor=safety_factor
        )
        
        if retry_count > 0:
            # Further reduce batch size on retry
            batch_size = max(100, int(batch_size * (0.5 ** retry_count)))
        
        logger.info(
            f"Auto-calculated batch size: {batch_size} sequences "
            f"(file_size: {file_size_mb:.1f}MB, "
            f"estimated_seqs: {estimated_seq_count}, "
            f"avg_seq_size: {avg_seq_size_mb:.3f}MB, "
            f"retry: {retry_count})"
        )
        
        return batch_size
    
    def hmmsearch(self):
        # Validate files before search
        if not self.hmm_file.exists():
            raise FileNotFoundError(f"HMM file not found: {self.hmm_file}")
        if not self.input_faa.exists():
            raise FileNotFoundError(f"Input protein file not found: {self.input_faa}")

        # Start timing
        start_time = time.time()
        
        cpus = max(1, min(self.hmmer_cpu, psutil.cpu_count() or 1))
        raw_hits_file = self.output_file.with_suffix(self.output_file.suffix + ".raw.tsv")
        
        # Initialize memory monitor
        memory_monitor = get_memory_monitor(
            max_memory_usage=getattr(self.config, 'max_memory_usage', 0.8)
        )
        
        # Start monitoring
        if getattr(self.config, 'enable_memory_monitoring', True):
            memory_monitor.start_monitoring()
            memory_monitor.log_memory_status("Before HMM search")
        
        # Statistics tracking
        stats = {
            'total_sequences': 0,
            'total_batches': 0,
            'total_hits': 0,
            'batch_size_history': [],
            'retry_count': 0,
            'memory_warnings': 0
        }
        
        logger.info(
            f"Running HMM search: hmm={self.hmm_file.name} input={self.input_faa.name} "
            f"out={self.output_file.name} evalue={self.e_value_threshold} "
            f"cov={self.coverage_threshold} cpus={cpus}"
        )

        # Memory / size heuristics
        input_size_bytes = self.input_faa.stat().st_size
        input_size_mb = input_size_bytes / (1024 * 1024)
        available_mb = memory_monitor.get_available_memory_mb()
        available_bytes = int(available_mb * 1024 * 1024)

        hmm_file_size_bytes = self.hmm_file.stat().st_size
        hmm_file_size_mb = hmm_file_size_bytes / (1024 * 1024)

        # Large mode switch:
        # - explicit: --large / large_mode=True
        # - implicit: input exceeds threshold
        large_mode = bool(getattr(self.config, "large_mode", False))
        large_threshold_mb = int(getattr(self.config, "large_input_threshold_mb", 5000) or 5000)
        if input_size_mb > large_threshold_mb:
            large_mode = True

        # Old pyhmmer recommendation: prefetch targets if input is small relative to available RAM.
        # Keep it for speed, but disable in large_mode for safety.
        preload_targets = (not large_mode) and (input_size_bytes < available_bytes * 0.10)

        # Optional: preload all HMMs into memory *only if* we have enough headroom.
        # HMM objects have significant overhead; estimate ~5x the file size (conservative).
        # Real-world observation: 2.35GB HMM file can use 10-15GB+ in memory.
        estimated_hmm_mem_bytes = int(hmm_file_size_bytes * 5.0)
        preload_hmms = (estimated_hmm_mem_bytes < available_bytes * 0.70)
        # In large_mode, be very conservative or disable preloading entirely for safety:
        if large_mode:
            # For large inputs, disable HMM preloading if HMM file is >1GB to avoid OOM risk
            # Even with abundant RAM, large HMM files can cause memory spikes during loading
            if hmm_file_size_mb > 1000:  # 1GB threshold
                preload_hmms = False
                logger.info(
                    f"Large mode + large HMM file ({hmm_file_size_mb:.1f}MB): "
                    f"disabling HMM preload for safety (will use streaming mode)"
                )
            else:
                # Only preload if we have massive headroom (10% threshold) to avoid OOM
                preload_hmms = (estimated_hmm_mem_bytes < available_bytes * 0.10)

        # Pre-compute constants to avoid repeated calculations
        hmm_file_stem = self.hmm_file.stem

        # Batch write buffer configuration (rows, not bytes)
        buffer_flush_threshold = int(getattr(self.config, 'csv_buffer_size', 5000) or 5000)
        
        max_retries = int(getattr(self.config, "max_retries", 3) or 3)
        last_oom: MemoryError | None = None

        for attempt in range(max_retries + 1):
            stats["retry_count"] = attempt
            if attempt > 0:
                # On retry: force the safest mode.
                effective_large_mode = True
                effective_preload_targets = False
                effective_preload_hmms = False
                # Reset counters because we will overwrite the raw hits file on retry
                stats["total_sequences"] = 0
                stats["total_batches"] = 0
                stats["total_hits"] = 0
                stats["memory_warnings"] = 0
            else:
                effective_large_mode = large_mode
                effective_preload_targets = preload_targets
                effective_preload_hmms = preload_hmms

            try:
                # Use larger buffer for file I/O (8MB buffer)
                with raw_hits_file.open("w", newline="", buffering=8 * 1024 * 1024) as raw_handle:
                    writer = csv.writer(raw_handle, delimiter="\t")
                    hit_buffer = []

                    alphabet = pyhmmer.easel.Alphabet.amino()

                    # Decide whether to force batching:
                    # - always in large_mode (safety)
                    # - or if user explicitly configured batch_size
                    force_batching = effective_large_mode or (getattr(self.config, "batch_size", None) is not None)

                    logger.info(
                        f"Input size: {input_size_mb:.1f}MB, HMM size: {hmm_file_size_mb:.1f}MB, "
                        f"Available: {available_mb:.1f}MB, large_mode={effective_large_mode}, "
                        f"preload_targets={effective_preload_targets}, preload_hmms={effective_preload_hmms}, "
                        f"force_batching={force_batching}, csv_buffer_size={buffer_flush_threshold}, "
                        f"retry={attempt}/{max_retries}"
                    )

                    # Decide HMMs strategy (may be used across multiple batches)
                    hmms = None  # Default to streaming
                    if effective_preload_hmms:
                        logger.info("Pre-loading all HMMs into memory for speed...")
                        logger.info(f"Estimated HMM memory: {estimated_hmm_mem_bytes / (1024*1024):.1f}MB")
                        hmms_list = []
                        hmm_count = 0
                        preload_success = True
                        try:
                            with pyhmmer.plan7.HMMFile(str(self.hmm_file)) as hmm_file_handle:
                                for hmm in hmm_file_handle:
                                    hmms_list.append(hmm)
                                    hmm_count += 1
                                    # Monitor memory every 1000 HMMs in large_mode
                                    if effective_large_mode and hmm_count % 1000 == 0:
                                        if getattr(self.config, 'enable_memory_monitoring', True):
                                            current_available = memory_monitor.get_available_memory_mb()
                                            if current_available < available_mb * 0.20:  # Less than 20% remaining
                                                logger.warning(
                                                    f"Memory dropping during HMM load ({current_available:.1f}MB remaining). "
                                                    f"Stopping preload at {hmm_count} HMMs to avoid OOM."
                                                )
                                                preload_success = False
                                                break
                            if preload_success:
                                hmms = hmms_list
                                logger.info(f"Loaded {len(hmms)} HMMs into memory.")
                            else:
                                logger.warning(
                                    f"HMM preload stopped early at {hmm_count} HMMs due to memory pressure. "
                                    f"Falling back to streaming HMM mode."
                                )
                                hmms = None
                        except MemoryError:
                            logger.error(
                                f"OOM during HMM preload after {len(hmms_list)} HMMs. "
                                f"Falling back to streaming HMM mode."
                            )
                            hmms = None

                    # Targets / batching loop
                    with pyhmmer.easel.SequenceFile(str(self.input_faa), digital=True, alphabet=alphabet) as seqs:
                        # If we explicitly preload targets, keep old fast path (single call)
                        if effective_preload_targets and not force_batching:
                            targets = seqs.read_block()
                            if hmms is not None:
                                for hits in self._iter_hmmsearch_hits(hmms, targets, cpus=cpus):
                                    for hit in hits:
                                        for domain in hit.domains.included:
                                            aln = domain.alignment
                                            coverage = (aln.hmm_to - aln.hmm_from + 1) / aln.hmm_length
                                            hmm_name = aln.hmm_name.decode("utf-8")
                                            if P.GT2_PREFIX in hmm_name:
                                                hmm_name = P.GT2_FAMILY_NAME
                                            i_evalue = domain.i_evalue
                                            if i_evalue < self.e_value_threshold and coverage > self.coverage_threshold:
                                                hit_buffer.append([
                                                    hmm_name,
                                                    aln.hmm_length,
                                                    aln.target_name.decode("utf-8"),
                                                    aln.target_length,
                                                    i_evalue,
                                                    aln.hmm_from,
                                                    aln.hmm_to,
                                                    aln.target_from,
                                                    aln.target_to,
                                                    coverage,
                                                    hmm_file_stem,
                                                ])
                                                stats["total_hits"] += 1
                                                if len(hit_buffer) >= buffer_flush_threshold:
                                                    self._flush_hit_buffer(hit_buffer, writer)
                                self._flush_hit_buffer(hit_buffer, writer)
                            else:
                                with pyhmmer.plan7.HMMFile(str(self.hmm_file)) as hmm_file_handle:
                                    for hits in self._iter_hmmsearch_hits(hmm_file_handle, targets, cpus=cpus):
                                        for hit in hits:
                                            for domain in hit.domains.included:
                                                aln = domain.alignment
                                                coverage = (aln.hmm_to - aln.hmm_from + 1) / aln.hmm_length
                                                hmm_name = aln.hmm_name.decode("utf-8")
                                                if P.GT2_PREFIX in hmm_name:
                                                    hmm_name = P.GT2_FAMILY_NAME
                                                i_evalue = domain.i_evalue
                                                if i_evalue < self.e_value_threshold and coverage > self.coverage_threshold:
                                                    hit_buffer.append([
                                                        hmm_name,
                                                        aln.hmm_length,
                                                        aln.target_name.decode("utf-8"),
                                                        aln.target_length,
                                                        i_evalue,
                                                        aln.hmm_from,
                                                        aln.hmm_to,
                                                        aln.target_from,
                                                        aln.target_to,
                                                        coverage,
                                                        hmm_file_stem,
                                                    ])
                                                    stats["total_hits"] += 1
                                                    if len(hit_buffer) >= buffer_flush_threshold:
                                                        self._flush_hit_buffer(hit_buffer, writer)
                                self._flush_hit_buffer(hit_buffer, writer)
                        else:
                            # True batching: read blocks of N sequences and run hmmsearch per block.
                            batch_size = self._calculate_batch_size(self.input_faa, memory_monitor, retry_count=stats["retry_count"])
                            stats["batch_size_history"].append(batch_size)
                            while True:
                                # Proactive: if memory is already above threshold, reduce batch size before reading more.
                                if getattr(self.config, 'enable_memory_monitoring', True):
                                    if not memory_monitor.is_memory_safe() and batch_size > 100:
                                        new_bs = max(100, batch_size // 2)
                                        if new_bs != batch_size:
                                            logger.warning(
                                                f"Memory pressure detected before reading next batch; "
                                                f"reducing batch_size {batch_size} -> {new_bs}"
                                            )
                                            batch_size = new_bs
                                            stats["batch_size_history"].append(batch_size)
                                            stats["memory_warnings"] += 1

                                seq_block = seqs.read_block(sequences=batch_size)
                                if not seq_block:
                                    break

                                stats["total_batches"] += 1
                                stats["total_sequences"] += len(seq_block)

                                if getattr(self.config, 'enable_memory_monitoring', True):
                                    if not memory_monitor.check_and_warn(f"batch {stats['total_batches']}"):
                                        stats["memory_warnings"] += 1
                                    # Periodic checkpoints for visibility/debugging (avoid per-batch overhead)
                                    if stats["total_batches"] % 50 == 0:
                                        memory_monitor.record_checkpoint(f"batch {stats['total_batches']}")

                                if hmms is not None:
                                    stats["total_hits"] += self._process_sequence_block(
                                        seq_block,
                                        hmms,
                                        cpus=cpus,
                                        hit_buffer=hit_buffer,
                                        hmm_file_stem=hmm_file_stem,
                                        buffer_flush_threshold=buffer_flush_threshold,
                                    )
                                    if len(hit_buffer) >= buffer_flush_threshold:
                                        self._flush_hit_buffer(hit_buffer, writer)
                                else:
                                    # Streaming HMMFile cannot be reused across multiple hmmsearch calls safely,
                                    # so reopen per batch.
                                    with pyhmmer.plan7.HMMFile(str(self.hmm_file)) as hmm_file_handle:
                                        stats["total_hits"] += self._process_sequence_block(
                                            seq_block,
                                            hmm_file_handle,
                                            cpus=cpus,
                                            hit_buffer=hit_buffer,
                                            hmm_file_stem=hmm_file_stem,
                                            buffer_flush_threshold=buffer_flush_threshold,
                                        )
                                        if len(hit_buffer) >= buffer_flush_threshold:
                                            self._flush_hit_buffer(hit_buffer, writer)

                            self._flush_hit_buffer(hit_buffer, writer)

                # success (no OOM)
                last_oom = None
                break
            except MemoryError as e:
                last_oom = e
                if attempt < max_retries:
                    logger.warning(
                        f"OOM during HMM search (retry {attempt+1}/{max_retries}). "
                        f"Will retry with safer settings and smaller batch size."
                    )
                    continue
                error_msg = (
                    f"HMM search failed due to out of memory (OOM) after {max_retries} retries. "
                    f"Try reducing batch_size or increasing available memory. "
                    f"File size: {input_size_mb:.1f}MB, Available: {available_mb:.1f}MB."
                )
                logger.error(error_msg)
                raise MemoryError(error_msg) from e
            except Exception as e:
                error_msg = f"HMM search failed for {self.hmm_file}: {e}"
                if "memory" in str(e).lower() or "oom" in str(e).lower():
                    error_msg += (
                        f" This may be a memory issue. "
                        f"File size: {input_size_mb:.1f}MB, Available: {available_mb:.1f}MB. "
                        f"Try reducing batch_size or increasing available memory."
                    )
                logger.error(error_msg)
                raise

        # Calculate elapsed time
        elapsed_time = time.time() - start_time
        
        # Generate and log comprehensive report
        if getattr(self.config, 'enable_memory_monitoring', True):
            memory_monitor.record_checkpoint("After HMM search")
            memory_monitor.log_report("HMM search")
        
        # Log performance statistics
        logger.info(
            f"{self.hmm_file.name} search completed. "
            f"Hits: {stats['total_hits']}, "
            f"Sequences: {stats['total_sequences']}, "
            f"Batches: {stats['total_batches']}, "
            f"Time: {elapsed_time:.2f}s, "
            f"Retries: {stats['retry_count']}, "
            f"Memory warnings: {stats['memory_warnings']}"
        )
        
        if stats['batch_size_history']:
            logger.debug(
                f"Batch size history: {stats['batch_size_history']} "
                f"(initial: {stats['batch_size_history'][0]}, "
                f"final: {stats['batch_size_history'][-1]})"
            )
        
        process_results(None, str(self.output_file), temp_hits_file=raw_hits_file)

    # -------- Orchestration --------
    def run(self):
        self.hmmsearch()


class PyHMMERDBCANProcessor(PyHMMERProcessor):

    def __init__(self, config: PyHMMERDBCANConfig):
        super().__init__(config)


class PyHMMERDBCANSUBProcessor(PyHMMERProcessor):
    def __init__(self, config: DBCANSUBConfig):
        super().__init__(config)

    @property
    def mapping_file(self) -> Path:
        return Path(self.config.db_dir) / P.SUBSTRATE_MAPPING_FILE

    def run(self):
        super().run()
        # Post-processing specific to dbCAN-sub
        sub_proc = DBCANSUBProcessor(self.config)
        sub_proc.process_dbcan_sub()


class PyHMMERTFProcessor(PyHMMERProcessor):
    def __init__(self, config: PyHMMERTFConfig):
        super().__init__(config)

    def run(self):
        if self.config.fungi:
            super().run()
        else:
            logger.info("TFProcessor: fungi=False, skipping TF HMM run.")

class PyHMMERSTPProcessor(PyHMMERProcessor):
    def __init__(self, config: PyHMMERSTPConfig):
        super().__init__(config)

class PyHMMERPfamProcessor(PyHMMERProcessor):
    def __init__(self, config: PyHMMERPfamConfig):
        super().__init__(config)

    @property
    def input_faa(self) -> Path:
        fname = P.NULL_PROTEIN_FILE if self.config.null_from_gff else P.INPUT_PROTEIN_FILE
        return Path(self.config.output_dir) / fname


