"""Constants for run_dbcan expression analysis."""

COUNT_MATRIX_FILE = "count_matrix.tsv"
GENE_NORM_COUNTS_FILE = "gene_norm_counts.tsv"
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

# CGC expression / heatmap y-axis metrics (--plot-metric / dbcan_plot --metric)
# Default log2_* scales are recommended for the expression track panel.
PLOT_METRIC_CHOICES = (
    "log2_tpm",
    "log2_fpkm",
    "log2_rpkm",
    "log2_rpm",
    "log2_cpm",
    "log2_norm",
    "log2fc",
    "gene_zscore",
    "tpm",
    "fpkm",
    "rpkm",
    "rpm",
    "counts",
)
DEFAULT_PLOT_METRIC = "log2_tpm"
