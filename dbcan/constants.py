import os
from pathlib import Path

TEST_ROOT = Path(__file__).parent
DATA_ROOT = os.path.join(TEST_ROOT, "_data")

INPUT_PROTEIN_NAME="uniInput.faa"
INPUT_PRODIGAL_GFF_NAME="uniInput.gff"



#################################### Constants for CGCFinder.py#####################
# CGCFinder related constants

# DataFrame column names
CGC_ANNOTATION_COLUMN = 'CGC_annotation'
PROTEIN_ID_COLUMN = 'Protein_ID'
CONTIG_ID_COLUMN = 'Contig ID'
START_COLUMN = 'start'
END_COLUMN = 'end'
STRAND_COLUMN = 'strand'
ATTRIBUTES_COLUMN = 'attributes'

GFF_COLUMNS = ['Contig ID', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']


# Gene marker columns
IS_CORE_COLUMN = 'is_core'
IS_ADDITIONAL_COLUMN = 'is_additional'
IS_SIGNATURE_COLUMN = 'is_signature'

# GFF file attribute names
CGC_ANNOTATION_ATTR = 'CGC_annotation'
PROTEIN_ID_ATTR = 'protein_id'

# Selected columns for CGC processing
CGC_SELECTED_COLUMNS = [CONTIG_ID_COLUMN, START_COLUMN, END_COLUMN, STRAND_COLUMN,
                        CGC_ANNOTATION_COLUMN, PROTEIN_ID_COLUMN]

# CGC output fields
CGC_ID_FIELD = 'CGC#'
CGC_PROTEIN_ID_FIELD = 'Protein ID'
GENE_TYPE_FIELD = 'Gene Type'
GENE_START_FIELD = 'Gene Start'
GENE_STOP_FIELD = 'Gene Stop'
GENE_STRAND_FIELD = 'Gene Strand'
GENE_ANNOTATION_FIELD = 'Gene Annotation'
NULL_GENE_TYPE = 'null'


CGC_CORE_SIG_TYPES = ['CAZyme']
CGC_DEFAULT_NULL_GENE = 2
CGC_DEFAULT_BP_DISTANCE = 15000
CGC_DEFAULT_USE_NULL_GENES = True
CGC_DEFAULT_USE_DISTANCE = False
CGC_DEFAULT_ADDITIONAL_GENES = ['TC']
CGC_GFF_FILE = 'cgc.gff'
CGC_RESULT_FILE = 'cgc_standard_out.tsv'
CGC_OUTPUT_COLUMNS = [CGC_ID_FIELD, GENE_TYPE_FIELD, CONTIG_ID_COLUMN, CGC_PROTEIN_ID_FIELD,
                    GENE_START_FIELD, GENE_STOP_FIELD, GENE_STRAND_FIELD, GENE_ANNOTATION_FIELD]


############################################################################################################


##########################cgc_substrate_prediction constant###############################

CAZYME="CAZyme"
TC="TC"
TF="TF"
STP="STP"
PUL="PUL"
NULL="null"

CGC_RESULT_FILE = CGC_RESULT_FILE
DBCAN_SUB_OUT_FILE = "dbCANsub_hmm_results.tsv"
OVERVIEW_FILE = "overview.tsv"
INPUT_PROTEIN_NAME = INPUT_PROTEIN_NAME

CGC_SUB_PREDICTION_FILE= "substrate_prediction.tsv"
PUL_DIAMOND_FILE = "PUL_blast.out"
CGC_FAA_FILE = "CGC.faa"
PUL_DIAMOND_DB = "PUL.dmnd"
PUL_EXCEL_FILE = "dbCAN-PUL.xlsx"
CAZYME_FAA_FILE = "CAZyme.faa"
PUL_FAA_FILE = "PUL.faa"



DBCANPUL_TMP="dbcanpul.tmp.txt"
DBCAN_SUB_TMP="dbcan_sub.tmp.txt"

DIAMOND_PUL_EVALUE = 0.01


################################################################################################


####################################### Constants for diamond.py ##############################

# Diamond database file names
CAZY_DIAMOND_DB = "CAZy.dmnd"
TCDB_DIAMOND_DB = "TCDB.dmnd"

# Output file names
CAZY_DIAMOND_OUTPUT = "diamond.out"
TCDB_DIAMOND_OUTPUT = "diamond.out.tc"

# Default parameters
DIAMOND_CAZY_EVALUE_DEFAULT = 1e-15
DIAMOND_TCDB_EVALUE_DEFAULT = 1e-5
DIAMOND_TCDB_COVERAGE_DEFAULT = 0.4
DIAMOND_MAX_TARGET_SEQS = "1"
DIAMOND_DEFAULT_OUTFMT = "6"

# Diamond command arguments
DIAMOND_CMD= "diamond"
DIAMOND_BLASTP_CMD = "blastp"
DIAMOND_CMD_DB = "--db"
DIAMOND_CMD_QUERY = "--query"
DIAMOND_CMD_OUT = "--out"
DIAMOND_CMD_OUTFMT = "--outfmt"
DIAMOND_CMD_EVALUE = "--evalue"
DIAMOND_CMD_MAX_TARGET = "--max-target-seqs"
DIAMOND_CMD_THREADS = "--threads"
DIAMOND_CMD_VERBOSE = "-v"
DIAMOND_CMD_QUIET = "--quiet"
DIAMOND_CMD_QUERY_COVER = "--query-cover"

# TCDB output format
TCDB_DIAMOND_OUTFMT_FIELDS = ['sseqid', 'slen', 'qseqid', 'qlen', 'evalue', 'sstart', 'send', 'qstart', 'qend', 'qcovhsp']

TCDB_ID_COLUMN = 'TCDB ID'


CAZY_COLUMN_NAMES = [
    'Gene ID',
    'CAZy ID',
    '% Identical',
    'Length',
    'Mismatches',
    'Gap Open',
    'Gene Start',
    'Gene End',
    'CAZy Start',
    'CAZy End',
    'E Value',
    'Bit Score'
]

TCDB_COLUMN_NAMES = [
    'TCDB ID',
    'TCDB Length',
    'Target ID',
    'Target Length',
    'EVALUE',
    'TCDB START',
    'TCDB END',
    'QSTART',
    'QEND',
    'COVERAGE'
]





################################################################################################

####################################### Constants for pyhmmer_search.py ##############################

HMMER_COLUMN_NAMES = [
    'HMM Name',
    'HMM Length',
    'Target Name',
    'Target Length',
    'i-Evalue',
    'HMM From',
    'HMM To',
    'Target From',
    'Target To',
    'Coverage',
    'HMM File Name'
]


# HMM Database files
DBCAN_HMM_FILE = "dbCAN.hmm"
DBCAN_SUB_HMM_FILE = "dbCAN-sub.hmm"
TF_HMM_FILE = "TF.hmm"
STP_HMM_FILE = "STP.hmm"

# Input/Output files
INPUT_PROTEIN_FILE = INPUT_PROTEIN_NAME
NON_CAZYME_PROTEIN_FILE = "non_CAZyme.faa"
DBCAN_HMM_RESULT_FILE = "dbCAN_hmm_results.tsv"
DBCAN_SUB_HMM_RESULT_FILE = "dbCANsub_hmm_results.tsv"
TF_HMM_RESULT_FILE = "TF_hmm_results.tsv"
STP_HMM_RESULT_FILE = "STP_hmm_results.tsv"

# Mapping files
SUBSTRATE_MAPPING_FILE = "fam-substrate-mapping.tsv"

# Special case handling
GT2_FAMILY_NAME = "GT2.hmm"
GT2_PREFIX = "GT2_"

################################################################################################


####################################### Constants for gff.py ##############################

# Input/Output file names
GFF_INPUT_PROTEIN_FILE = INPUT_PROTEIN_NAME
GFF_CAZYME_OVERVIEW_FILE = OVERVIEW_FILE
GFF_CGC_SIG_FILE = "total_cgc_info.tsv"
GFF_OUTPUT_FILE = CGC_GFF_FILE
GFF_TEMP_SUFFIX = ".temp"

# Column names and indices
GFF_PROTEIN_ID_COL = "protein_id"
GFF_CAZYME_COL = "CAZyme"
GFF_GENE_ID_COL = "Gene ID"
GFF_TOOLS_COUNT_COL = "#ofTools"
GFF_RECOMMEND_RESULTS_COL = "Recommend Results"
GFF_CGC_ANNOTATION_COL = CGC_ANNOTATION_COLUMN
GFF_FUNCTION_ANNOTATION_COL = "function_annotation"
GFF_TYPE_COL = "type"
GFF_CGC_SIG_COLUMNS = [0, 2, 10]
GFF_MIN_TOOL_COUNT = 2

# Annotation prefixes and defaults
GFF_CAZYME_PREFIX = CAZYME + "|"
GFF_OTHER_PREFIX = "Other|"
GFF_NULL_ANNOTATION = NULL
GFF_UNKNOWN_ANNOTATION = "unknown"
GFF_NA_PROTEIN_ID = "NA"

# GFF feature types
GFF_GENE_FEATURE = "gene"
GFF_MRNA_FEATURE = "mRNA"
GFF_CDS_FEATURE = "CDS"

# GFF format types
GFF_FORMAT_NCBI_EUK = "NCBI_euk"
GFF_FORMAT_NCBI_PROK = "NCBI_prok"
GFF_FORMAT_JGI = "JGI"
GFF_FORMAT_PRODIGAL = "prodigal"

# GFF attribute names
GFF_PROTEIN_ID_ATTR = PROTEIN_ID_ATTR
GFF_NAME_ATTR = "Name"
GFF_ID_ATTR = "ID"
GFF_JGI_PROTEIN_ID_ATTR = "proteinId"

################################################################################################

####################################### Constants for OverviewGenerator.py ##############################

# File names
OVERVIEW_FILE = OVERVIEW_FILE
NON_CAZYME_FAA_FILE = NON_CAZYME_PROTEIN_FILE

# Result file names
DIAMOND_RESULT_FILE = CAZY_DIAMOND_OUTPUT
DBCAN_SUB_RESULT_FILE = DBCAN_SUB_HMM_RESULT_FILE
DBCAN_HMM_RESULT_FILE = DBCAN_HMM_RESULT_FILE

# Column names and data structures
OVERVIEW_COLUMNS = ['Gene ID', 'EC#', 'dbCAN_hmm', 'dbCAN_sub', 'DIAMOND', '#ofTools', 'Recommend Results']
DIAMOND_COLUMN_NAMES_OVERVIEW = ['Gene ID', 'CAZy ID']
DBCAN_SUB_COLUMN_NAMES_OVERVIEW = ['Target Name', 'Subfam Name', 'Subfam EC', 'Target From', 'Target To', 'i-Evalue']
DBCAN_HMM_COLUMN_NAMES_OVERVIEW = ['Target Name', 'HMM Name', 'Target From', 'Target To', 'i-Evalue']

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

# Configuration values
OVERVIEW_OVERLAP_THRESHOLD = 0.5
MIN_TOOLS_FOR_RECOMMENDATION = 2

# Regex patterns
CAZY_ID_PATTERN = r"^(GH|GT|CBM|AA|CE|PL)"

# Separators
RESULT_SEPARATOR = "+"
EC_SEPARATOR = "|"
RANGE_SEPARATOR = "-"


####################################### Constants for plot_cgc_circle.py ##############################

# File paths and names
CGC_GFF_FILE = CGC_GFF_FILE
CGC_RESULT_FILE = CGC_RESULT_FILE
CGC_CIRCOS_DIR = "cgc_circos"
CGC_CIRCOS_PLOT_FILE = "cgc_circos_plot.svg"
CGC_CIRCOS_CONTIG_FILE_TEMPLATE = "cgc_circos_{contig_name}.svg"
DEG_FILE = "DEG.tsv"

# Feature types
CGC_FEATURE_TYPE = "gene"
CGC_ANNOTATION_ATTR = CGC_ANNOTATION_ATTR
PROTEIN_ID_ATTR = PROTEIN_ID_ATTR

# TSV column names
CGC_ID_COLUMN = CGC_ID_FIELD
CONTIG_ID_COLUMN = CONTIG_ID_COLUMN
CGC_PROTEIN_ID_FIELD = CGC_PROTEIN_ID_FIELD
GENE_START_COLUMN = GENE_START_FIELD
GENE_STOP_COLUMN = GENE_STOP_FIELD

# Circos track parameters
CGC_OUTER_TRACK_RANGE = (86.7, 87)
CGC_CAZYME_TRACK_RANGE = (35, 40)
CGC_FEATURE_TRACK_RANGE = (45, 50)
CGC_RANGE_TRACK_RANGE = (52, 55)
DEG_TRACK_RANGE = (65, 70)
DEG_LOG2FC_RANGE=(75,85)

CGC_TRACK_PADDING = 0.1
CGC_MAJOR_INTERVAL = 100000
CGC_MINOR_INTERVAL_DIVISOR = 10

# Visual properties
CGC_TRACK_BG_COLOR = "#EEEEEE"
CGC_GRID_COLOR = "black"
CGC_RANGE_COLOR = "lightblue"
CGC_RANGE_BORDER_COLOR = "black"
CGC_AXIS_COLOR = "black"
CGC_LABEL_SIZE = 10
CGC_LEGEND_POSITION = (0.5, 0.4)
CGC_LEGEND_FONT_SIZE = 20
CGC_TITLE_FONT_SIZE = 40

# Feature colors
CGC_FEATURE_LEGEND = ["CAZyme", "TC", "TF", "STP", "PEPTIDASE", "SULFATLAS"]
CGC_FEATURE_COLORS = {
    "CAZyme": "#E67E22",      # orange
    "TC": "#2ECC71",          # green
    "TF": "#9B59B6",          # purple
    "STP": "#F1C40F",         # golden yellow
    "PEPTIDASE": "#16A085",   # greenish
    "SULFATLAS": "#34495E",   # dark blue
    "default": "#95A5A6"      # light gray
}


# Plot scaling parameters
CGC_MIN_FIGURE_SIZE = 15
CGC_MAX_FIGURE_SIZE = 30
CGC_FIGURE_SIZE_SCALING_FACTOR = 0.5

# Text constants
CGC_PLOT_TITLE = "CGC Annotation Circos Plot"
CGC_CONTIG_TITLE_TEMPLATE = "CGC Annotation - {contig_name}"
CGC_LEGEND_TITLE = "Types"




####################################### Constants for process_dbcan_sub.py ##############################

# File paths
DBCAN_SUB_RESULT_FILE = DBCAN_SUB_HMM_RESULT_FILE
SUBSTRATE_MAPPING_FILE = SUBSTRATE_MAPPING_FILE



# Column names for input data
DBCAN_SUB_HMM_NAME_COLUMN = "HMM Name"
DBCAN_SUB_TARGET_NAME_COLUMN = "Target Name"
DBCAN_SUB_TARGET_LENGTH_COLUMN = "Target Length"
DBCAN_SUB_IEVALUE_COLUMN = "i-Evalue"
DBCAN_SUB_HMM_LENGTH_COLUMN = "HMM Length"
DBCAN_SUB_HMM_FROM_COLUMN = "HMM From"
DBCAN_SUB_HMM_TO_COLUMN = "HMM To"
DBCAN_SUB_TARGET_FROM_COLUMN = "Target From"
DBCAN_SUB_TARGET_TO_COLUMN = "Target To"
DBCAN_SUB_COVERAGE_COLUMN = "Coverage"
DBCAN_SUB_HMM_FILE_COLUMN = "HMM File Name"

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

# dbCAN_sub_COLUMN_NAMES = [
#     'Subfam Name',
#     'Subfam Composition',
#     'Subfam EC',
#     'Substrate',
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

# Special family prefixes
DBCAN_SUB_CBM_PREFIX = "CBM"

# File formats
DBCAN_SUB_HMM_SUFFIX = ".hmm"
DBCAN_SUB_SEPARATOR = "|"



#########constants for process_utils.py##############################
OVERLAP_RATIO_THRESHOLD = 0.5

#############constants for downloading dbCAN databases##############################


CAZY_DB_URL = "https://bcb.unl.edu/dbCAN2/download/run_dbCAN_database_total/CAZy.dmnd"
HMMER_DB_URL = "https://bcb.unl.edu/dbCAN2/download/run_dbCAN_database_total/dbCAN.hmm"
DBCAN_SUB_DB_URL = "https://bcb.unl.edu/dbCAN2/download/run_dbCAN_database_total/dbCAN_sub.hmm"
DBCAN_SUB_MAP_URL = "https://bcb.unl.edu/dbCAN2/download/run_dbCAN_database_total/fam-substrate-mapping.tsv"

TCDB_DB_URL = "https://bcb.unl.edu/dbCAN2/download/run_dbCAN_database_total/tcdb.dmnd"
TF_DB_URL = "https://bcb.unl.edu/dbCAN2/download/run_dbCAN_database_total/TF.hmm"
STP_DB_URL = "https://bcb.unl.edu/dbCAN2/download/run_dbCAN_database_total/STP.hmm"

PEPTIDASE_DB_URL = "https://bcb.unl.edu/dbCAN2/download/run_dbCAN_database_total/peptidase_db.dmnd"
SULFATLAS_DB_URL = "https://bcb.unl.edu/dbCAN2/download/run_dbCAN_database_total/sulfatlas_db.dmnd"

PUL_DB_URL = "https://bcb.unl.edu/dbCAN2/download/run_dbCAN_database_total/PUL.dmnd"
PUL_MAP_URL = "https://bcb.unl.edu/dbCAN2/download/run_dbCAN_database_total/dbCAN-PUL.xlsx"
PUL_ALL_URL ="https://bcb.unl.edu/dbCAN2/download/run_dbCAN_database_total/dbCAN-PUL.tar.gz"



FILE_PATHS = {
    'diamond': 'diamond_results.tsv',
    'dbcan_sub': 'dbCAN-sub.substrate.tsv',
    'dbcan_hmm': 'dbCAN_hmm_results.tsv'
}

COLUMN_NAMES = {
    'diamond': ['Gene ID', 'CAZy ID'],
    'dbcan_sub': ['Target Name', 'Subfam Name', 'Subfam EC', 'Target From', 'Target To', 'i-Evalue'],
    'dbcan_hmm': ['Target Name', 'HMM Name', 'Target From', 'Target To', 'i-Evalue']
}



#############constants for newly added sulfatlas and peptidase database ##############################
# Sulfatlas Constants
SULFATLAS = "SULFATLAS"
SULFATLAS_DIAMOND_DB = "sulfatlas_db.dmnd"
SULFATLAS_DIAMOND_OUTPUT = "diamond.out.sulfatlas"
SULFATLAS_DIAMOND_OUTFMT_FIELDS = ['sseqid', 'slen', 'qseqid', 'qlen', 'evalue', 'sstart', 'send', 'qstart', 'qend', 'qcovhsp']
SULFATLAS_ID_COLUMN = 'Sul ID'

SULFATLAS_COLUMN_NAMES = [
    'Sul ID',
    'Sul Length',
    'Target ID',
    'Target Length',
    'EVALUE',
    'Sul START',
    'Sul END',
    'QSTART',
    'QEND',
    'COVERAGE'
]
DIAMOND_SULFATLAS_EVALUE_DEFAULT = 1e-5
DIAMOND_SULFATLAS_COVERAGE_DEFAULT = 0.3

# Peptidase Constants
PEPTIDASE = "PEPTIDASE"
PEPTIDASE_DIAMOND_DB = "peptidase_db.dmnd"
PEPTIDASE_DIAMOND_OUTPUT = "diamond.out.peptidase"
PEPTIDASE_DIAMOND_OUTFMT_FIELDS = ['sseqid', 'slen', 'qseqid', 'qlen', 'evalue', 'sstart', 'send', 'qstart', 'qend', 'qcovhsp']
PEPTIDASE_ID_COLUMN = 'Peptidase ID'

PEPTIDASE_COLUMN_NAMES = [
    'Peptidase ID',
    'Peptidase Length',
    'Target ID',
    'Target Length',
    'EVALUE',
    'Peptidase START',
    'Peptidase END',
    'QSTART',
    'QEND',
    'COVERAGE'
]
DIAMOND_PEPTIDASE_EVALUE_DEFAULT = 1e-5
DIAMOND_PEPTIDASE_COVERAGE_DEFAULT = 0.3
