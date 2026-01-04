from dataclasses import dataclass
from typing import Optional
import psutil
import dbcan.constants.pyhmmer_search_constants as P
from dbcan.configs.base_config import BaseConfig

@dataclass
class PyHMMERConfig(BaseConfig):
    db_dir: str
    output_dir: str
    threads: int = psutil.cpu_count()
    # unified common fields
    input_faa: str = P.INPUT_PROTEIN_FILE
    hmm_file: str = None
    output_file: str = None
    evalue_threshold: Optional[float] = None
    coverage_threshold: Optional[float] = None
    # Memory management options
    batch_size: Optional[int] = None  # Number of sequences per batch (auto-calculated if None)
    max_memory_usage: float = 0.8  # Maximum memory usage ratio (0.0-1.0)
    enable_memory_monitoring: bool = True  # Enable memory monitoring
    memory_safety_factor: float = 0.5  # Safety factor for batch size calculation (0.0-1.0)

# dbCAN
@dataclass
class PyHMMERDBCANConfig(PyHMMERConfig):
    hmm_file: str = P.DBCAN_HMM_FILE
    output_file: str = P.DBCAN_HMM_RESULT_FILE
    evalue_threshold: float = P.DBCAN_EVALUE_DEFAULT
    coverage_threshold: float = P.DBCAN_COVERAGE_DEFAULT

# dbCAN-sub
@dataclass
class DBCANSUBConfig(PyHMMERConfig):
    hmm_file: str = P.DBCAN_SUB_HMM_FILE
    output_file: str = P.DBCAN_SUB_HMM_RAW_FILE
    output_sub_file: str = P.DBCAN_SUB_HMM_RESULT_FILE
    evalue_threshold: float = P.DBCAN_SUB_EVALUE_DEFAULT
    coverage_threshold: float = P.DBCAN_SUB_COVERAGE_DEFAULT
    mapping_file: Optional[str] = None

# TF
@dataclass
class PyHMMERTFConfig(PyHMMERConfig):
    hmm_file: str = P.TF_HMM_FILE
    output_file: str = P.TF_HMM_RESULT_FILE
    evalue_threshold: float = P.TF_EVALUE_DEFAULT
    coverage_threshold: float = P.TF_COVERAGE_DEFAULT
    fungi: bool = False

# STP
@dataclass
class PyHMMERSTPConfig(PyHMMERConfig):
    hmm_file: str = P.STP_HMM_FILE
    output_file: str = P.STP_HMM_RESULT_FILE
    evalue_threshold: float = P.STP_EVALUE_DEFAULT
    coverage_threshold: float = P.STP_COVERAGE_DEFAULT

# Pfam
@dataclass
class PyHMMERPfamConfig(PyHMMERConfig):
    hmm_file: str = P.PFAM_HMM_FILE
    output_file: str = P.PFAM_HMM_RESULT_FILE
    evalue_threshold: float = P.PFAM_EVALUE_DEFAULT
    coverage_threshold: float = P.PFAM_COVERAGE_DEFAULT
    null_from_gff: bool = False

