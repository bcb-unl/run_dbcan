import dbcan.constants.base_constants as base_constants

####################################### Constants for diamond.py ##############################

INPUT_PROTEIN_NAME = base_constants.INPUT_PROTEIN_NAME

# Diamond database file names
CAZY_DIAMOND_DB = "CAZy.dmnd"
TCDB_DIAMOND_DB = "TCDB.dmnd"

# Output file names
CAZY_DIAMOND_OUTPUT = base_constants.CAZY_DIAMOND_OUTPUT
TCDB_DIAMOND_OUTPUT = base_constants.TCDB_DIAMOND_OUTPUT

# Default parameters
DIAMOND_CAZY_EVALUE_DEFAULT = 1e-102
#DIAMOND_CAZY_COVERAGE_DEFAULT = 0
DIAMOND_TCDB_EVALUE_DEFAULT = 1e-4
DIAMOND_TCDB_COVERAGE_DEFAULT = 35
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
DIAMOND_CMD_VERBOSE = "--verbose"
DIAMOND_CMD_QUIET = "--quiet"
DIAMOND_CMD_QUERY_COVER = "--query-cover"

# TCDB output format
TCDB_DIAMOND_OUTFMT_FIELDS = ['sseqid', 'slen', 'qseqid', 'qlen', 'evalue', 'sstart', 'send', 'qstart', 'qend', 'qcovhsp']

TCDB_ID_COLUMN = 'TCDB ID'
TC=base_constants.TC

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

#############constants for newly added sulfatlas and peptidase database ##############################
# Sulfatlas Constants
SULFATASE = base_constants.SULFATASE
SULFATLAS_DIAMOND_DB = "sulfatlas_db.dmnd"
SULFATLAS_DIAMOND_OUTPUT = base_constants.SULFATLAS_DIAMOND_OUTPUT
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
DIAMOND_SULFATLAS_EVALUE_DEFAULT = 1e-4
DIAMOND_SULFATLAS_COVERAGE_DEFAULT = 35

# Peptidase Constants
PEPTIDASE = base_constants.PEPTIDASE
PEPTIDASE_DIAMOND_DB = "peptidase_db.dmnd"
PEPTIDASE_DIAMOND_OUTPUT = base_constants.PEPTIDASE_DIAMOND_OUTPUT
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
DIAMOND_PEPTIDASE_EVALUE_DEFAULT = 1e-4
DIAMOND_PEPTIDASE_COVERAGE_DEFAULT = 35


#############constants for newly added TF database ##############################
# TF Constants
TF = base_constants.TF
TF_DIAMOND_DB = "TF.dmnd"
TF_DIAMOND_OUTPUT = base_constants.TF_DIAMOND_OUTPUT
TF_DIAMOND_OUTFMT_FIELDS = ['sseqid', 'slen', 'qseqid', 'qlen', 'evalue', 'sstart', 'send', 'qstart', 'qend', 'qcovhsp']
TF_ID_COLUMN = 'TF ID'
TF_DATABASE = "prodoric"
TF_COLUMN_NAMES = [
    'TF ID',
    'TF Length',
    'Target ID',
    'Target Length',
    'EVALUE',
    'TF START',
    'TF END',
    'QSTART',
    'QEND',
    'COVERAGE'
]
DIAMOND_TF_EVALUE_DEFAULT = 1e-4
DIAMOND_TF_COVERAGE_DEFAULT = 35

################################################################################################
