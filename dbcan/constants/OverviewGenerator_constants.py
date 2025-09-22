import dbcan.constants.base_constants as base_constants

# Constants for Overview Generator


# File names
INPUT_PRODIGAL_GFF_NAME = base_constants.INPUT_PRODIGAL_GFF_NAME
INPUT_PROTEIN_NAME = base_constants.INPUT_PROTEIN_NAME
OVERVIEW_FILE = base_constants.OVERVIEW_FILE

# Result file names
DIAMOND_RESULT_FILE = base_constants.CAZY_DIAMOND_OUTPUT
DBCAN_SUB_RESULT_FILE = base_constants.DBCAN_SUB_HMM_RESULT_FILE
DBCAN_HMM_RESULT_FILE = base_constants.DBCAN_HMM_RESULT_FILE

# Column names and data structures
OVERVIEW_COLUMNS = ['Gene ID', 'EC#', 'dbCAN_hmm', 'dbCAN_sub', 'DIAMOND', '#ofTools', 'Recommend Results']
DIAMOND_COLUMN_NAMES_OVERVIEW = ['Gene ID', 'CAZy ID']
DBCAN_SUB_COLUMN_NAMES_OVERVIEW = ['Target Name', 'Subfam Name', 'Subfam EC', 'Target From', 'Target To', 'i-Evalue']
DBCAN_HMM_COLUMN_NAMES_OVERVIEW = ['Target Name', 'HMM Name', 'Target From', 'Target To', 'i-Evalue']

DBCAN_SUB_SUBSTRATE_COLUMN = "Substrate"

# Special fields and values
GENE_ID_FIELD = "Gene ID"
EC_FIELD = "EC#"
DBCAN_HMM_FIELD = "dbCAN_hmm"
DBCAN_SUB_FIELD = "dbCAN_sub"
DIAMOND_FIELD = "DIAMOND"
TOOLS_COUNT_FIELD = "#ofTools"
RECOMMEND_RESULTS_FIELD = "Recommend Results"
EMPTY_RESULT_PLACEHOLDER = "-"
SUBFAMILY_NAME_FIELD = "Subfam Name"
HMM_NAME_FIELD = "HMM Name"
TARGET_NAME_FIELD = "Target Name"
TARGET_FROM_FIELD = "Target From"
TARGET_TO_FIELD = "Target To"
I_EVALUE_FIELD = "i-Evalue"
CAZY_ID_FIELD = "CAZy ID"
SUBFAMILY_EC_FIELD = "Subfam EC"
SUBSTRATE_FIELD = "Substrate"

# Configuration values
OVERVIEW_OVERLAP_THRESHOLD = 0.5
MIN_TOOLS_FOR_RECOMMENDATION = 2

# Regex patterns
CAZY_ID_PATTERN = r"^(GH|GT|CBM|AA|CE|PL)"

# Separators
RESULT_SEPARATOR = "+"
EC_SEPARATOR = "|"
SUB_SEPARATOR = ";"
RANGE_SEPARATOR = "-"
