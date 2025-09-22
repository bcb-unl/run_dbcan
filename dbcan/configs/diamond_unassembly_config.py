from dataclasses import dataclass
from typing import Optional

@dataclass
class DiamondAsmfreeConfig:
    # Which function to run
    function: str  # diamond_fam_abund | diamond_subfam_abund | diamond_EC_abund | diamond_substrate_abund

    # Inputs
    paf1: str
    paf2: Optional[str] = None               # Optional second PAF (paired-end)
    input: Optional[str] = None              # For downstream table-based steps (EC/substrate)
    db_dir: str = "./db"                     # Mapping tables directory
    raw_reads: Optional[str] = None          # Raw reads file (fq/fa or gz)

    # Outputs
    output: str = "asmfree_fam_abund"        # Output file path

    # Normalization: FPKM | RPM | TPM
    normalized: str = "TPM"
