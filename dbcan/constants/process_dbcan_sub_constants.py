import re
import dbcan.constants.base_constants as base_constants
import dbcan.constants.process_utils_constants as process_utils_constants

####################################### Constants for process_dbcan_sub.py ##############################

# File paths
DBCAN_SUB_HMM_RAW_FILE = base_constants.DBCAN_SUB_HMM_RAW_FILE
SUBSTRATE_MAPPING_FILE = 'fam-substrate-mapping.tsv'

DBCAN_SUB_HMM_RESULT_FILE = base_constants.DBCAN_SUB_HMM_RESULT_FILE

# Column names for input data
# DBCAN_SUB_HMM_NAME_COLUMN = "HMM Name"
# DBCAN_SUB_TARGET_NAME_COLUMN = "Target Name"
# DBCAN_SUB_TARGET_LENGTH_COLUMN = "Target Length"
# DBCAN_SUB_IEVALUE_COLUMN = "i-Evalue"
# DBCAN_SUB_HMM_LENGTH_COLUMN = "HMM Length"
# DBCAN_SUB_HMM_FROM_COLUMN = "HMM From"
# DBCAN_SUB_HMM_TO_COLUMN = "HMM To"
# DBCAN_SUB_TARGET_FROM_COLUMN = "Target From"
# DBCAN_SUB_TARGET_TO_COLUMN = "Target To"
# DBCAN_SUB_COVERAGE_COLUMN = "Coverage"
# DBCAN_SUB_HMM_FILE_COLUMN = "HMM File Name"

DBCAN_SUB_HMM_NAME_COLUMN = process_utils_constants.HMM_NAME_COLUMN
DBCAN_SUB_TARGET_NAME_COLUMN = process_utils_constants.TARGET_NAME_COLUMN
DBCAN_SUB_TARGET_LENGTH_COLUMN = process_utils_constants.TARGET_LENGTH_COLUMN
DBCAN_SUB_IEVALUE_COLUMN = process_utils_constants.IEVALUE_COLUMN
DBCAN_SUB_HMM_LENGTH_COLUMN = process_utils_constants.HMM_LENGTH_COLUMN
DBCAN_SUB_HMM_FROM_COLUMN = process_utils_constants.HMM_FROM_COLUMN
DBCAN_SUB_HMM_TO_COLUMN = process_utils_constants.HMM_TO_COLUMN
DBCAN_SUB_TARGET_FROM_COLUMN = process_utils_constants.TARGET_FROM_COLUMN
DBCAN_SUB_TARGET_TO_COLUMN = process_utils_constants.TARGET_TO_COLUMN
DBCAN_SUB_COVERAGE_COLUMN = process_utils_constants.COVERAGE_COLUMN
DBCAN_SUB_HMM_FILE_COLUMN = process_utils_constants.HMM_FILE_COLUMN



# Column names for processed data
DBCAN_SUB_SUBFAMILY_NAME_COLUMN = "Subfam Name"
DBCAN_SUB_SUBFAMILY_COMP_COLUMN = "Subfam Composition"
DBCAN_SUB_SUBFAMILY_EC_COLUMN = "Subfam EC"
DBCAN_SUB_SUBSTRATE_COLUMN = "Substrate"

# Collection of all columns in final output
DBCAN_SUB_COLUMN_NAMES = [
    DBCAN_SUB_SUBFAMILY_NAME_COLUMN,
    DBCAN_SUB_SUBFAMILY_COMP_COLUMN,
    DBCAN_SUB_SUBFAMILY_EC_COLUMN,
    DBCAN_SUB_SUBSTRATE_COLUMN,
    DBCAN_SUB_HMM_LENGTH_COLUMN,
    DBCAN_SUB_TARGET_NAME_COLUMN,
    DBCAN_SUB_TARGET_LENGTH_COLUMN,
    DBCAN_SUB_IEVALUE_COLUMN,
    DBCAN_SUB_HMM_FROM_COLUMN,
    DBCAN_SUB_HMM_TO_COLUMN,
    DBCAN_SUB_TARGET_FROM_COLUMN,
    DBCAN_SUB_TARGET_TO_COLUMN,
    DBCAN_SUB_COVERAGE_COLUMN,
    DBCAN_SUB_HMM_FILE_COLUMN
]

# Special family prefixes
DBCAN_SUB_CBM_PREFIX = "CBM"

# File formats
DBCAN_SUB_HMM_SUFFIX = ".hmm"
DBCAN_SUB_SEPARATOR = "|"


# EC number core inside dbCAN-sub eHMM names: a.b.c.d where b/c/d may be '-' and d may be 'nX'
EC_CORE_PATTERN = re.compile(r"^\d+\.(?:\d+|-)\.(?:\d+|-)\.(?:n?\d+|-)$")

# fam-substrate-mapping.tsv columns
MAPPING_SUBSTRATE_COLUMN = "Substrate_high_level"
MAPPING_NEW_SUBSTRATE_COLUMN = "new_Substrate_high_level"
MAPPING_FAMILY_COLUMN = "Family"
MAPPING_EC_COLUMN = "EC_Number"
MAPPING_TYPE_COLUMN = "type"

# CAZy family token in an eHMM composition: GH5_5 -> GH5
FAMILY_TOKEN_PATTERN = re.compile(r"^((?:GH|GT|PL|CE|AA|CBM)\d+)(?:_\d+)?$")
# Optional per-subfamily substrate table exported from the dbCAN-sub web server
DBCANSUB_SUBSTRATE_TABLE_FILE = "dbCAN-sub-substrate.tsv"
