"""Config for structure-based CAZyme annotation (Foldseek vs CAZyme3D)."""
from dataclasses import dataclass, field
from typing import List

from dbcan.configs.base_config import BaseConfig
import dbcan.constants.structure_search_constants as S


@dataclass
class StructureSearchConfig(BaseConfig):
    """Configuration for Foldseek search against CAZyme3D structure database."""

    db_dir: str = None
    output_dir: str = None
    threads: int = 1
    input_faa: str = S.INPUT_PROTEIN_NAME
    output_file: str = S.STRUCTURE_SEARCH_RESULT_FILE
    db_file: str = S.CAZYME3D_FOLDSEEK_DB_PREFIX
    mapping_file: str = S.CAZYME3D_MAPPING_FILE
    e_value_threshold: float = S.STRUCTURE_SEARCH_EVALUE_DEFAULT
    column_names: List[str] = field(default_factory=lambda: list(S.STRUCTURE_SEARCH_COLUMN_NAMES_OVERVIEW))
