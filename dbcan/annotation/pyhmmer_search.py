from __future__ import annotations

import logging
from pathlib import Path
from abc import ABC
import csv
import psutil
import pyhmmer
import json
import time
from tqdm import tqdm
from typing import Dict, List

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

    # -------- Core search --------
    def _process_sequence_block(self, seq_block, hmm_file_handle, cpus: int, writer: csv.writer) -> int:
        """Process a sequence block (streaming mode, no file I/O).
        
        Args:
            seq_block: DigitalSequenceBlock from read_block()
            hmm_file_handle: HMM file handle
            cpus: Number of CPUs to use
            writer: CSV writer for results
        
        Returns:
            int: Number of hits found
        """
        hit_count = 0
        
        try:
            for hits in pyhmmer.hmmsearch(
                hmm_file_handle,
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
                            writer.writerow([
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
                                self.hmm_file.stem
                            ])
                            hit_count += 1
        except Exception as e:
            logger.error(f"Error processing sequence block: {e}")
            raise
        
        return hit_count
    
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

        # Check if we need batch processing
        input_size_mb = self.input_faa.stat().st_size / (1024 * 1024)
        available_mb = memory_monitor.get_available_memory_mb()
        
        # Use batch processing if file is large or memory is limited
        # Threshold: if file size > 5% of available memory, use batching
        use_batching = (input_size_mb > available_mb * 0.05) or (input_size_mb > 100)  # 100MB threshold
        
        try:
            with raw_hits_file.open("w", newline="") as raw_handle:
                writer = csv.writer(raw_handle, delimiter='\t')
                
                if use_batching:
                    # Streaming batch processing mode (no temporary files)
                    logger.info("Using streaming batch processing mode to avoid OOM")
                    max_retries = 3
                    retry_count = 0
                    hit_count = 0
                    alphabet = pyhmmer.easel.Alphabet.amino()
                    
                    while retry_count <= max_retries:
                        try:
                            batch_size = self._calculate_batch_size(self.input_faa, memory_monitor, retry_count)
                            stats['batch_size_history'].append(batch_size)
                            
                            with pyhmmer.plan7.HMMFile(str(self.hmm_file)) as hmm_file_handle:
                                with pyhmmer.easel.SequenceFile(str(self.input_faa), digital=True, alphabet=alphabet) as seqs:
                                    total_hits = 0
                                    batch_num = 0
                                    
                                    # Estimate total batches for progress bar
                                    total_seqs_estimate = max(1, int(input_size_mb / 0.01))  # Rough estimate
                                    estimated_batches = (total_seqs_estimate + batch_size - 1) // batch_size
                                    
                                    logger.info(f"Processing sequences in streaming mode (batch_size={batch_size})")
                                    
                                    with tqdm(total=estimated_batches, desc="Processing sequences", unit="batch") as pbar:
                                        while True:
                                            try:
                                                # Read a block of sequences (streaming, no file I/O)
                                                seq_block = seqs.read_block(sequences=batch_size)
                                                
                                                # Check if block is empty (end of file)
                                                if seq_block is None or len(seq_block) == 0:
                                                    break
                                                
                                                batch_num += 1
                                                stats['total_batches'] = batch_num
                                                stats['total_sequences'] += len(seq_block)
                                                
                                                # Check memory before processing
                                                if getattr(self.config, 'enable_memory_monitoring', True):
                                                    memory_monitor.record_checkpoint(f"batch {batch_num}")
                                                    if not memory_monitor.check_and_warn(f"batch {batch_num}"):
                                                        stats['memory_warnings'] += 1
                                                        logger.warning(
                                                            f"Memory usage high during batch {batch_num}. "
                                                            f"Consider reducing batch_size if OOM occurs."
                                                        )
                                                
                                                # Process the sequence block directly (no file I/O)
                                                batch_hits = self._process_sequence_block(
                                                    seq_block, hmm_file_handle, cpus, writer
                                                )
                                                total_hits += batch_hits
                                                
                                                pbar.update(1)
                                                pbar.set_postfix({"hits": total_hits, "batch": batch_num})
                                                
                                                # Clear the block from memory immediately
                                                del seq_block
                                                import gc
                                                gc.collect()
                                                
                                            except MemoryError as e:
                                                logger.warning(f"Memory error in batch {batch_num}, reducing batch size and retrying")
                                                # Reduce batch size and retry
                                                batch_size = max(100, int(batch_size * 0.5))
                                                stats['batch_size_history'].append(batch_size)
                                                stats['retry_count'] += 1
                                                logger.info(f"Reduced batch size to {batch_size}, retrying...")
                                                retry_count += 1
                                                if retry_count > max_retries:
                                                    raise
                                                # Note: We can't rewind the sequence file, so we continue with reduced size
                                                continue
                                            except StopIteration:
                                                # End of file
                                                break
                                            except Exception as e:
                                                logger.error(f"Error processing batch {batch_num}: {e}")
                                                import traceback
                                                logger.error(f"Traceback: {traceback.format_exc()}")
                                                raise
                                    
                                    hit_count = total_hits
                                    stats['total_hits'] = total_hits
                                    break  # Success, exit retry loop
                        
                        except MemoryError as e:
                            retry_count += 1
                            if retry_count <= max_retries:
                                error_msg = (
                                    f"Out of memory (OOM) during streaming processing (attempt {retry_count}/{max_retries}). "
                                    f"Retrying with reduced batch size. "
                                    f"File size: {input_size_mb:.1f}MB, Available: {available_mb:.1f}MB"
                                )
                                logger.warning(error_msg)
                                # Force garbage collection
                                import gc
                                gc.collect()
                            else:
                                error_msg = (
                                    f"HMM search failed due to out of memory (OOM) after {max_retries} retries. "
                                    f"Try reducing batch_size manually or increasing available memory. "
                                    f"File size: {input_size_mb:.1f}MB, Available: {available_mb:.1f}MB. "
                                    f"Suggested batch_size: {max(100, int(batch_size * 0.25))}"
                                )
                                logger.error(error_msg)
                                raise MemoryError(error_msg) from e
                else:
                    # Single-pass mode for small files
                    logger.info("Using single-pass mode (file is small enough)")
                    alphabet = pyhmmer.easel.Alphabet.amino()
                    hit_count = 0
                    
                    with pyhmmer.plan7.HMMFile(str(self.hmm_file)) as hmm_file_handle:
                        with pyhmmer.easel.SequenceFile(str(self.input_faa), digital=True, alphabet=alphabet) as seqs:
                            for hits in pyhmmer.hmmsearch(
                                hmm_file_handle,
                                seqs,
                                cpus=cpus,
                                domE=self.e_value_threshold
                            ):
                                for hit in hits:
                                    stats['total_sequences'] += 1
                                    for domain in hit.domains.included:
                                        aln = domain.alignment
                                        coverage = (aln.hmm_to - aln.hmm_from + 1) / aln.hmm_length
                                        hmm_name = aln.hmm_name.decode('utf-8')
                                        if P.GT2_PREFIX in hmm_name:
                                            hmm_name = P.GT2_FAMILY_NAME
                                        i_evalue = domain.i_evalue
                                        if i_evalue < self.e_value_threshold and coverage > self.coverage_threshold:
                                            writer.writerow([
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
                                                self.hmm_file.stem
                                            ])
                                            hit_count += 1
                    stats['total_hits'] = hit_count
        except MemoryError as e:
            # This should only be reached if retries are exhausted or in single-pass mode
            error_msg = (
                f"HMM search failed due to out of memory (OOM). "
                f"Try reducing batch_size or increasing available memory. "
                f"File size: {input_size_mb:.1f}MB, Available: {available_mb:.1f}MB. "
                f"Consider using batch_size parameter."
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


