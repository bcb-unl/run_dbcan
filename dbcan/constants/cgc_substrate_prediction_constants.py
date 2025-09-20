import dbcan.constants.base_constants as base_constants
##########################cgc_substrate_prediction constant###############################

# CAZYME="CAZyme"
# TC="TC"
# TF="TF"
# STP="STP"
# PUL="PUL"
# NULL="null"

CAZYME=base_constants.CAZYME
TC=base_constants.TC
TF=base_constants.TF
STP=base_constants.STP
PUL=base_constants.PUL
NULL=base_constants.NULL



CGC_RESULT_FILE = base_constants.CGC_RESULT_FILE
DBCAN_SUB_OUT_FILE = base_constants.DBCAN_SUB_HMM_RESULT_FILE
OVERVIEW_FILE = base_constants.OVERVIEW_FILE
INPUT_PROTEIN_NAME = base_constants.INPUT_PROTEIN_NAME

CGC_SUB_PREDICTION_FILE= base_constants.CGC_SUB_PREDICTION_FILE
PUL_DIAMOND_FILE = base_constants.PUL_DIAMOND_FILE
CGC_FAA_FILE = base_constants.CGC_FAA_FILE
PUL_DIAMOND_DB = base_constants.PUL_DIAMOND_DB
PUL_EXCEL_FILE = base_constants.PUL_EXCEL_FILE
CAZYME_FAA_FILE = base_constants.CAZYME_FAA_FILE
PUL_FAA_FILE = base_constants.PUL_FAA_FILE

# CGC_RESULT_FILE = 'cgc_standard_out.tsv'
# DBCAN_SUB_OUT_FILE = "dbCANsub_hmm_results.tsv"
# OVERVIEW_FILE = "overview.tsv"
# INPUT_PROTEIN_NAME = "uniInput.faa"

# CGC_SUB_PREDICTION_FILE= "substrate_prediction.tsv"
# PUL_DIAMOND_FILE = "PUL_blast.out"
# CGC_FAA_FILE = "CGC.faa"
# PUL_DIAMOND_DB = "PUL.dmnd"
# PUL_EXCEL_FILE = "dbCAN-PUL.xlsx"
# CAZYME_FAA_FILE = "CAZyme.faa"
# PUL_FAA_FILE = "PUL.faa"

DBCANPUL_TMP="dbcanpul.tmp.txt"
DBCAN_SUB_TMP="dbcan_sub.tmp.txt"

DIAMOND_PUL_EVALUE = 0.01

# ===================== Default parameter constants (extracted from configs) =====================
# Homology / dbCAN-PUL matching thresholds
DEFAULT_UNIQ_PUL_GENE_HIT_NUM = 2          # upghn
DEFAULT_UNIQ_QUERY_CGC_GENE_NUM = 2        # uqcgn
DEFAULT_CAZYME_PAIR_NUM = 1                # cpn
DEFAULT_TOTAL_PAIR_NUM = 2                 # tpn
DEFAULT_EXTRA_PAIR_TYPE = None             # e.g. "TC-TC,STP-STP"
DEFAULT_EXTRA_PAIR_TYPE_NUM = "0"          # e.g. "1,2"

DEFAULT_IDENTITY_CUTOFF = 0.0
DEFAULT_COVERAGE_CUTOFF = 0.0
DEFAULT_BITSCORE_CUTOFF = 50.0
DEFAULT_EVALUE_CUTOFF = 1e-2

# dbCAN-sub substrate inference thresholds
DEFAULT_SUB_HMM_EVALUE = 1e-2
DEFAULT_SUB_HMM_COVERAGE = 0.0
DEFAULT_SUB_DOMAIN_SUBSTRATE_CUTOFF = 2
DEFAULT_SUB_PROTEIN_SUBSTRATE_CUTOFF = 2
DEFAULT_SUBSTRATE_SCORE_THRESHOLD = 2.0

# Execution / runtime
DEFAULT_THREADS = 0               # 0 => auto detect
DEFAULT_OUTPUT_DBCAN_SUB = False
DEFAULT_OUTPUT_DBCAN_PUL = False

################################################################################################
