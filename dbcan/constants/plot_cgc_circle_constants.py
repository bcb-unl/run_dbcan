import dbcan.constants.base_constants as base_constants

# File paths and names
CGC_GFF_FILE = base_constants.CGC_GFF_FILE
CGC_RESULT_FILE = base_constants.CGC_RESULT_FILE
CGC_CIRCOS_DIR = "cgc_circos"
CGC_CIRCOS_PLOT_FILE = "cgc_circos_plot.svg"
CGC_CIRCOS_CONTIG_FILE_TEMPLATE = "cgc_circos_{contig_name}.svg"
DEG_FILE = "DEG.tsv"

# Feature types
#CGC_FEATURE_TYPE = ["gene", "CDS"]
CGC_ANNOTATION_ATTR = base_constants.GFF_CGC_ANNOTATION_COL
PROTEIN_ID_ATTR = base_constants.GFF_PROTEIN_ID_COL

# TSV column names
CGC_ID_COLUMN = base_constants.CGC_ID_FIELD
CONTIG_ID_COLUMN = base_constants.CONTIG_ID_COLUMN
CGC_PROTEIN_ID_FIELD = base_constants.CGC_PROTEIN_ID_FIELD
GENE_START_COLUMN = base_constants.GENE_START_FIELD
GENE_STOP_COLUMN = base_constants.GENE_STOP_FIELD

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
# CGC_FEATURE_LEGEND = ["CAZyme", "TC", "TF", "STP", "PEPTIDASE", "SULFATLAS"]
# CGC_FEATURE_COLORS = {
#     "CAZyme": "#E67E22",      # orange
#     "TC": "#2ECC71",          # green
#     "TF": "#9B59B6",          # purple
#     "STP": "#F1C40F",         # golden yellow
#     "PEPTIDASE": "#16A085",   # greenish
#     "SULFATLAS": "#34495E",   # dark blue
#     "default": "#95A5A6"      # light gray
# }
CGC_FEATURE_LEGEND=base_constants.CGC_FEATURE_LEGEND
CGC_FEATURE_COLORS=base_constants.CGC_FEATURE_COLORS

# Plot scaling parameters
CGC_MIN_FIGURE_SIZE = 15
CGC_MAX_FIGURE_SIZE = 30
CGC_FIGURE_SIZE_SCALING_FACTOR = 0.5

# Text constants
CGC_PLOT_TITLE = "CGC Annotation Circos Plot"
CGC_CONTIG_TITLE_TEMPLATE = "CGC Annotation - {contig_name}"
CGC_LEGEND_TITLE = "Types"



