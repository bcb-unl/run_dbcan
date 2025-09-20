import dbcan.constants.base_constants as base_constants

"""Constants for SignalP annotation (DeepTMHMM removed)"""

# File names
OVERVIEW_FILE = base_constants.OVERVIEW_FILE
INPUT_PROTEIN_NAME = base_constants.INPUT_PROTEIN_NAME
TOPOLOGY_INPUT_FILE = "topology_input_recommended.faa"
SIGNALP_RESULT_TSV = "prediction_results.txt"

# Column names
SIGNALP_COL = "SignalP"
GENE_ID_COL = base_constants.GFF_GENE_ID_COL
RECOMMEND_RESULTS_COL = base_constants.GFF_RECOMMEND_RESULTS_COL

# Default values
DEFAULT_EMPTY = "-"
DEFAULT_SIGNALP_ORG = "other"
DEFAULT_MODE = "fast"
DEFAULT_FORMAT = "none"

# Output directories
SIGNALP_OUT_DIR = "signalp6_out"

# SignalP parameters
SIGNALP_ORGANISMS = "other"
SIGNALP_MODES = "fast"
OUTPUT_FORMATS = "none"

BIOLIB_SIGNALP_VARIANTS = [
    "DTU/SignalP-6.0",
    "signalp",
    "signalp6"
]