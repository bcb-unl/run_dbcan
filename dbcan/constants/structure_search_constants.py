"""Constants for structure-based CAZyme annotation (Foldseek vs CAZyme3D)."""
import dbcan.constants.base_constants as base_constants

# Input (same as other CAZyme annotation methods)
INPUT_PROTEIN_NAME = base_constants.INPUT_PROTEIN_NAME

# CAZyme3D / Foldseek database file names (in db_dir)
CAZYME3D_FOLDSEEK_DB_PREFIX = "CAZyme3D"
CAZYME3D_MAPPING_FILE = "cazyme3d_to_cazy.tsv"

# Output file name (written under output_dir)
STRUCTURE_SEARCH_RESULT_FILE = "structure_search_results.tsv"

# Output columns for overview compatibility (same semantics as Diamond: Gene ID, CAZy ID)
STRUCTURE_SEARCH_COLUMN_NAMES_OVERVIEW = ["Gene ID", "CAZy ID"]

# Default e-value threshold for Foldseek hits
STRUCTURE_SEARCH_EVALUE_DEFAULT = 1e-3

# Foldseek binary name for availability check
FOLDSEEK_CMD = "foldseek"
