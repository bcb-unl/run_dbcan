from dataclasses import dataclass, field
from typing import List, Optional
import dbcan.constants.diamond_constants as D
from dbcan.configs.base_config import BaseConfig

@dataclass
class DiamondConfig(BaseConfig):
    db_dir: str
    output_dir: str
    threads: int = 1
    verbose_option: bool = False
    diamond_bin: str = D.DIAMOND_CMD
    # input_faa: relative to output_dir when not absolute
    input_faa: str = D.INPUT_PROTEIN_NAME
    e_value_threshold: float = 1e-5
    coverage_threshold: Optional[float] = None


@dataclass
class DiamondCAZyConfig(DiamondConfig):
    db_file: str = D.CAZY_DIAMOND_DB
    output_file: str = D.CAZY_DIAMOND_OUTPUT
    outfmt_fields: Optional[List[str]] = None
    column_names: List[str] = field(default_factory=lambda: list(D.CAZY_COLUMN_NAMES))
    id_column: Optional[str] = None
    label: Optional[str] = None
    e_value_threshold: float = D.DIAMOND_CAZY_EVALUE_DEFAULT
    #coverage_threshold: float = D.DIAMOND_CAZY_COVERAGE_DEFAULT

@dataclass
class DiamondTCConfig(DiamondConfig):
    db_file: str = D.TCDB_DIAMOND_DB
    output_file: str = D.TCDB_DIAMOND_OUTPUT
    outfmt_fields: Optional[List[str]] = field(default_factory=lambda: list(D.TCDB_DIAMOND_OUTFMT_FIELDS))
    column_names: List[str] = field(default_factory=lambda: list(D.TCDB_COLUMN_NAMES))
    id_column: Optional[str] = D.TCDB_ID_COLUMN
    label: Optional[str] = D.TC
    e_value_threshold: float = D.DIAMOND_TCDB_EVALUE_DEFAULT
    coverage_threshold: float = D.DIAMOND_TCDB_COVERAGE_DEFAULT

@dataclass
class DiamondSulfataseConfig(DiamondConfig):
    db_file: str = D.SULFATLAS_DIAMOND_DB
    output_file: str = D.SULFATLAS_DIAMOND_OUTPUT
    outfmt_fields: Optional[List[str]] = field(default_factory=lambda: list(D.SULFATLAS_DIAMOND_OUTFMT_FIELDS))
    column_names: List[str] = field(default_factory=lambda: list(D.SULFATLAS_COLUMN_NAMES))
    id_column: Optional[str] = D.SULFATLAS_ID_COLUMN
    label: Optional[str] = D.SULFATASE
    e_value_threshold: float = D.DIAMOND_SULFATLAS_EVALUE_DEFAULT
    coverage_threshold: float = D.DIAMOND_SULFATLAS_COVERAGE_DEFAULT

@dataclass
class DiamondPeptidaseConfig(DiamondConfig):
    db_file: str = D.PEPTIDASE_DIAMOND_DB
    output_file: str = D.PEPTIDASE_DIAMOND_OUTPUT
    outfmt_fields: Optional[List[str]] = field(default_factory=lambda: list(D.PEPTIDASE_DIAMOND_OUTFMT_FIELDS))
    column_names: List[str] = field(default_factory=lambda: list(D.PEPTIDASE_COLUMN_NAMES))
    id_column: Optional[str] = D.PEPTIDASE_ID_COLUMN
    label: Optional[str] = D.PEPTIDASE
    e_value_threshold: float = D.DIAMOND_PEPTIDASE_EVALUE_DEFAULT
    coverage_threshold: float = D.DIAMOND_PEPTIDASE_COVERAGE_DEFAULT

@dataclass
class DiamondTFConfig(DiamondConfig):
    db_file: str = D.TF_DIAMOND_DB
    output_file: str = D.TF_DIAMOND_OUTPUT
    outfmt_fields: Optional[List[str]] = field(default_factory=lambda: list(D.TF_DIAMOND_OUTFMT_FIELDS))
    column_names: List[str] = field(default_factory=lambda: list(D.TF_COLUMN_NAMES))
    id_column: Optional[str] = D.TF_ID_COLUMN
    label: Optional[str] = D.TF_DATABASE
    e_value_threshold: float = D.DIAMOND_TF_EVALUE_DEFAULT
    coverage_threshold: float = D.DIAMOND_TF_COVERAGE_DEFAULT
    prokaryotic: bool = True  # whether to use prokaryotic TF database
