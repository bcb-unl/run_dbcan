"""Constants for run_dbcan expression analysis."""

COUNT_MATRIX_FILE = "count_matrix.tsv"
SAMPLE_METADATA_FILE = "sample_metadata.tsv"
GENE_DEG_FILE = "gene_deg.tsv"
DEG_FILE = "DEG.tsv"
CGC_DEG_FILE = "cgc_deg.tsv"
CGC_GENE_DEG_MAP_FILE = "cgc_gene_deg_map.tsv"
CGC_EXPRESSION_MATRIX_FILE = "cgc_expression_matrix.tsv"
CGC_GENE_COUNTS_FILE = "cgc_gene_counts.tsv"

PLOTS_DIR = "plots"
CGC_EXPRESSION_PLOTS_DIR = "CGC_expression"
DEG_VOLCANO_PDF = "DEG_volcano.pdf"
DEG_CGC_BAR_PDF = "DEG_CGC_bar.pdf"

DEFAULT_ALPHA = 0.05
DEFAULT_LFC_THRESHOLD = 1.0
DEFAULT_MAX_CGC = 500
DEFAULT_CGC_DE_RULE = "any"

CGC_DE_RULES = ("any", "majority", "all")
