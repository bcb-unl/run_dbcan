from __future__ import annotations

import logging
from pathlib import Path
from abc import ABC
import psutil
import pyhmmer

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
    def hmmsearch(self):
        # Validate files before search
        if not self.hmm_file.exists():
            raise FileNotFoundError(f"HMM file not found: {self.hmm_file}")
        if not self.input_faa.exists():
            raise FileNotFoundError(f"Input protein file not found: {self.input_faa}")

        available_memory = psutil.virtual_memory().available
        target_size = self.input_faa.stat().st_size
        results = []
        cpus = max(1, min(self.hmmer_cpu, psutil.cpu_count() or 1))

        logger.info(
            f"Running HMM search: hmm={self.hmm_file.name} input={self.input_faa.name} "
            f"out={self.output_file.name} evalue={self.e_value_threshold} "
            f"cov={self.coverage_threshold} cpus={cpus}"
        )

        try:
            with pyhmmer.plan7.HMMFile(str(self.hmm_file)) as hmm_file_handle:
                with pyhmmer.easel.SequenceFile(str(self.input_faa), digital=True) as seqs:
                    targets = seqs.read_block() if target_size < available_memory * 0.1 else seqs
                    for hits in pyhmmer.hmmsearch(
                        hmm_file_handle,
                        targets,
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
                                    results.append([
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
        except Exception as e:
            logger.error(f"HMM search failed for {self.hmm_file}: {e}")
            raise

        logger.info(f"{self.hmm_file.name} search completed. Hits: {len(results)}")
        process_results(results, str(self.output_file))

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


