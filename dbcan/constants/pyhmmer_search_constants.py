

import dbcan.constants.process_utils_constants as process_utils_constants
import dbcan.constants.process_dbcan_sub_constants as process_dbcan_sub_constants
import dbcan.constants.base_constants as base_constants


DBCAN_EVALUE_DEFAULT = 1e-15
DBCAN_COVERAGE_DEFAULT = 0.35


DBCAN_SUB_EVALUE_DEFAULT = 1e-15
DBCAN_SUB_COVERAGE_DEFAULT = 0.35
TF_EVALUE_DEFAULT = 1e-4
TF_COVERAGE_DEFAULT = 0.35
STP_EVALUE_DEFAULT = 1e-4
STP_COVERAGE_DEFAULT = 0.35
PFAM_EVALUE_DEFAULT = 1e-4
PFAM_COVERAGE_DEFAULT = 0.35

HMMER_COLUMN_NAMES=process_utils_constants.HMMER_COLUMN_NAMES

# HMMER_COLUMN_NAMES = [
#     'HMM Name',
#     'HMM Length',
#     'Target Name',
#     'Target Length',
#     'i-Evalue',
#     'HMM From',
#     'HMM To',
#     'Target From',
#     'Target To',
#     'Coverage',
#     'HMM File Name'
# ]


# HMM Database files
DBCAN_HMM_FILE = "dbCAN.hmm"
DBCAN_SUB_HMM_FILE = "dbCAN-sub.hmm"
TF_HMM_FILE = "TF.hmm"
STP_HMM_FILE = "STP.hmm"
PFAM_HMM_FILE = "Pfam-A.hmm"

# Input/Output files
INPUT_PROTEIN_FILE = base_constants.INPUT_PROTEIN_NAME
#NON_CAZYME_PROTEIN_FILE = "non_CAZyme.faa"
NULL_PROTEIN_FILE = base_constants.NULL_PROTEIN_FILE
DBCAN_HMM_RESULT_FILE = base_constants.DBCAN_HMM_RESULT_FILE
DBCAN_SUB_HMM_RESULT_FILE = base_constants.DBCAN_SUB_HMM_RESULT_FILE
TF_HMM_RESULT_FILE = base_constants.TF_HMM_RESULT_FILE
STP_HMM_RESULT_FILE = base_constants.STP_HMM_RESULT_FILE
PFAM_HMM_RESULT_FILE = base_constants.PFAM_HMM_RESULT_FILE

# Mapping files
SUBSTRATE_MAPPING_FILE = process_dbcan_sub_constants.SUBSTRATE_MAPPING_FILE
DBCAN_SUB_HMM_RAW_FILE = base_constants.DBCAN_SUB_HMM_RAW_FILE

# Special case handling
GT2_FAMILY_NAME = "GT2.hmm"
GT2_PREFIX = "GT2_"


# Overlap ratio threshold for filtering overlapping hits
OVERLAP_RATIO_THRESHOLD = 0.5