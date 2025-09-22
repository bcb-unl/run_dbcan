import dbcan.constants.base_constants as base_constants


#################################### Constants for CGCFinder.py#####################
# CGCFinder related constants



# DataFrame column names
# CGC_ANNOTATION_COLUMN = 'CGC_annotation'
# PROTEIN_ID_COLUMN = 'Protein_ID'
# CONTIG_ID_COLUMN = 'Contig ID'
# START_COLUMN = 'start'
# END_COLUMN = 'end'
# STRAND_COLUMN = 'strand'
# ATTRIBUTES_COLUMN = 'attributes'

CGC_ANNOTATION_COLUMN = base_constants.CGC_ANNOTATION_COLUMN
PROTEIN_ID_COLUMN = base_constants.PROTEIN_ID_COLUMN
CONTIG_ID_COLUMN = base_constants.CONTIG_ID_COLUMN
START_COLUMN = base_constants.START_COLUMN
END_COLUMN = base_constants.END_COLUMN
STRAND_COLUMN = base_constants.STRAND_COLUMN
ATTRIBUTES_COLUMN = base_constants.ATTRIBUTES_COLUMN
GFF_COLUMNS = base_constants.GFF_COLUMNS

# Gene marker columns
IS_CORE_COLUMN = 'is_core'
IS_ADDITIONAL_COLUMN = 'is_additional'
IS_SIGNATURE_COLUMN = 'is_signature'

# GFF file attribute names
CGC_ANNOTATION_ATTR = base_constants.CGC_ANNOTATION_COLUMN
PROTEIN_ID_ATTR = base_constants.GFF_PROTEIN_ID_COL

# Selected columns for CGC processing
CGC_SELECTED_COLUMNS = [CONTIG_ID_COLUMN, START_COLUMN, END_COLUMN, STRAND_COLUMN,
                        CGC_ANNOTATION_COLUMN, PROTEIN_ID_COLUMN]

# CGC output fields
# CGC_ID_FIELD = 'CGC#'
# CGC_PROTEIN_ID_FIELD = 'Protein ID'
# GENE_TYPE_FIELD = 'Gene Type'
# GENE_START_FIELD = 'Gene Start'
# GENE_STOP_FIELD = 'Gene Stop'
# GENE_STRAND_FIELD = 'Gene Strand'
# GENE_ANNOTATION_FIELD = 'Gene Annotation'
# NULL_GENE_TYPE = 'null'

CGC_ID_FIELD = base_constants.CGC_ID_FIELD
CGC_PROTEIN_ID_FIELD = base_constants.CGC_PROTEIN_ID_FIELD
GENE_TYPE_FIELD = base_constants.GENE_TYPE_FIELD
GENE_START_FIELD = base_constants.GENE_START_FIELD
GENE_STOP_FIELD = base_constants.GENE_STOP_FIELD
GENE_STRAND_FIELD = base_constants.GENE_STRAND_FIELD
GENE_ANNOTATION_FIELD = base_constants.GENE_ANNOTATION_FIELD
NULL_GENE_TYPE = base_constants.NULL


CGC_CORE_SIG_TYPES = ['CAZyme']
CGC_DEFAULT_NULL_GENE = 2
CGC_DEFAULT_BP_DISTANCE = 15000
CGC_DEFAULT_USE_NULL_GENES = True
CGC_DEFAULT_USE_DISTANCE = False
CGC_DEFAULT_ADDITIONAL_GENES = ['TC']

# NEW: defaults for CGC extension behavior
CGC_DEFAULT_EXTEND_MODE = 'none'       # 'none' | 'bp' | 'gene'
CGC_DEFAULT_EXTEND_BP = 0
CGC_DEFAULT_EXTEND_GENE_COUNT = 0


CGC_DEFAULT_MIN_CORE_CAZYME = 1
CGC_DEFAULT_MIN_CLUSTER_GENES = 2


CGC_GFF_FILE = base_constants.CGC_GFF_FILE
CGC_RESULT_FILE = base_constants.CGC_RESULT_FILE
CGC_OUTPUT_COLUMNS = [CGC_ID_FIELD, GENE_TYPE_FIELD, CONTIG_ID_COLUMN, CGC_PROTEIN_ID_FIELD,
                    GENE_START_FIELD, GENE_STOP_FIELD, GENE_STRAND_FIELD, GENE_ANNOTATION_FIELD]

PRIORITY = {'CAZyme': 0, 'TC': 1, 'TF': 2, 'STP': 3, 'SULFATLAS': 4, 'PEPTIDASE': 5}