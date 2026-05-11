import dbcan.constants.base_constants as base_constants

"""Constants for optional SignalP and local DeepTMHMM annotation."""

# File names
OVERVIEW_FILE = base_constants.OVERVIEW_FILE
INPUT_PROTEIN_NAME = base_constants.INPUT_PROTEIN_NAME
TOPOLOGY_INPUT_FILE = "topology_input_recommended.faa"
SIGNALP_RESULT_TSV = "prediction_results.txt"
DEEPTMHMM_RESULT_3LINE = "predicted_topologies.3line"
DEEPTMHMM_PREDICT_SCRIPT = "predict.py"

# Column names
SIGNALP_COL = "SignalP"
DEEPTMHMM_COL = "DeepTMHMM"
GENE_ID_COL = base_constants.GFF_GENE_ID_COL
RECOMMEND_RESULTS_COL = base_constants.GFF_RECOMMEND_RESULTS_COL

# Default values
DEFAULT_EMPTY = "-"
DEFAULT_SIGNALP_ORG = "other"
DEFAULT_MODE = "fast"
DEFAULT_FORMAT = "none"
DEFAULT_DEEPTMHMM_PYTHON = "python"

# Output directories
SIGNALP_OUT_DIR = "signalp6_out"
DEEPTMHMM_OUT_DIR = "deeptmhmm_out"

# SignalP parameters (must be sequences for "x in CONSTANT" validation)
SIGNALP_ORGANISMS = ("other", "euk")
SIGNALP_MODES = ("fast",)
OUTPUT_FORMATS = ("none",)

BIOLIB_SIGNALP_VARIANTS = [
    "DTU/SignalP-6.0",
    "signalp",
    "signalp6"
]