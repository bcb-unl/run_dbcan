from dataclasses import dataclass
from typing import Optional

from dbcan.configs.base_config import BaseConfig
from dbcan.constants.expression_constants import (
    DEFAULT_ALPHA,
    DEFAULT_CGC_DE_RULE,
    DEFAULT_LFC_THRESHOLD,
    DEFAULT_MAX_CGC,
)


@dataclass
class ExpressionConfig(BaseConfig):
    input_dir: str = ""
    output_dir: str = "output.expression"
    samplesheet: str = ""
    reference_fasta: str = ""
    gff: str = ""
    design: str = "~condition"
    abundance: str = "TPM"
    overlap_base_ratio: float = 0.2
    mapping_quality: int = 30
    identity: float = 0.98
    threads: int = 1
    hifi: bool = False
    alpha: float = DEFAULT_ALPHA
    lfc_threshold: float = DEFAULT_LFC_THRESHOLD
    cgc_de_rule: str = DEFAULT_CGC_DE_RULE
    all_genes: bool = False
    skip_alignment: bool = False
    skip_deseq2: bool = False
    run_plots: bool = False
    run_plots_only_de: bool = False
    also_heatmap: bool = False
    max_cgc: int = DEFAULT_MAX_CGC
    force: bool = False
    db_dir: str = "db"
