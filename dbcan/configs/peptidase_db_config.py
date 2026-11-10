from dataclasses import dataclass
from typing import Optional

import dbcan.constants.peptidase_db_constants as P
from dbcan.configs.base_config import BaseConfig


@dataclass
class PeptidaseDbBuildConfig(BaseConfig):
    """Rebuild peptidase_db.dmnd from MEROPS pepunit.lib (dbcan.utils.build_peptidase_db)."""
    # MEROPS pepunit.lib (plain or .gz)
    input: str
    # DIAMOND database to write (peptidase_db.dmnd)
    output: str = P.DEFAULT_OUTPUT_DB
    # cleaned FASTA; None: next to the database, '<output stem>.faa'
    fasta: Optional[str] = None
    # optional TSV with one row per entry that was left out (accession, MEROPS id, reason, header)
    report: Optional[str] = None
    # True: keep MEROPS inhibitor units (valid I identifiers) as decoys; run_dbcan drops hits to them
    keep_inhibitors: bool = True   # inhibitor units act as decoys (see prepare_the_database.rst)
    # False: write the FASTA only (no 'diamond makedb')
    makedb: bool = True
    diamond_bin: str = P.DIAMOND_CMD
    threads: int = 1
