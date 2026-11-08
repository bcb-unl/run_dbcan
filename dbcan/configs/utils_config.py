from dataclasses import dataclass
from typing import Optional

import dbcan.constants.utils_constants as U


@dataclass
class AbundanceConfig:
    """dbcan_utils fam_abund / fam_substrate_abund / CGC_abund / CGC_substrate_abund."""
    input_dir: str                      # run_dbcan output directory
    bedtools_depth: str                 # gene read-count table written by cal_coverage (id, length, count)
    output_dir: str = "."               # every *.out file is written here
    abundance: str = U.DEFAULT_ABUNDANCE_METHOD   # FPKM/RPKM/RPM/TPM


@dataclass
class CoverageConfig:
    """dbcan_utils cal_coverage: count reads/fragments per gene from a BAM and a GFF."""
    input: str                          # coordinate-sorted, indexed BAM
    gff: str
    output: str                         # gene read-count table (id, length, count)
    overlap_base_ratio: float = U.DEFAULT_OVERLAP_BASE_RATIO
    mapping_quality: int = U.DEFAULT_MAPPING_QUALITY
    identity: float = U.DEFAULT_IDENTITY
    threads: int = 1
    hifi: bool = False
    count_mode: str = U.DEFAULT_COUNT_MODE
    strandedness: str = U.DEFAULT_STRANDEDNESS
    include_duplicates: bool = False
    id_attr: str = U.GFF_ID_ATTR_AUTO
    feature_type: Optional[str] = None  # comma-separated GFF column-3 types; None = automatic
