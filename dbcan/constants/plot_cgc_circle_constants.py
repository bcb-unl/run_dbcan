import dbcan.constants.base_constants as base_constants

# File paths and names
CGC_GFF_FILE = base_constants.CGC_GFF_FILE
CGC_RESULT_FILE = base_constants.CGC_RESULT_FILE
CGC_SUB_PREDICTION_FILE = base_constants.CGC_SUB_PREDICTION_FILE
CGC_CIRCOS_DIR = "cgc_circos"
CGC_CIRCOS_FILE_TEMPLATE = "cgc_circos_plot.{ext}"
CGC_OVERVIEW_FILE_TEMPLATE = "cgc_overview.{ext}"
CGC_PLOT_FORMATS = ("png", "pdf")
DEG_FILE = "DEG.tsv"

# --layout of run_dbcan cgc_circle_plot: auto = circos for one contig, overview otherwise
CGC_LAYOUT_CHOICES = ("auto", "circos", "overview")
CGC_DEFAULT_LAYOUT = "auto"

# Circos rings (radius 0-100); the CGC ring is outermost so its labels sit outside the circle
CGC_AXIS_TRACK_RANGE = (99.0, 99.3)
CGC_RANGE_TRACK_RANGE = (90, 95)
CGC_FEATURE_TRACK_RANGE = (82, 88)
CGC_CAZYME_TRACK_RANGE = (74, 80)
DEG_TRACK_RANGE = (66, 72)
DEG_LOG2FC_RANGE = (50, 64)
CGC_LABEL_MIN_R = 101
CGC_LABEL_MAX_R = 108

CGC_TRACK_BG_COLOR = "#F2F2F2"
CGC_RANGE_COLOR = base_constants.CGC_FEATURE_COLORS["Other"]
CGC_RANGE_PREDICTED_COLOR = base_constants.CGC_FEATURE_COLORS["CAZyme"]
CGC_LABEL_SIZE = 7
CGC_TICK_LABEL_SIZE = 7
CGC_FIGURE_SIZE = 9.0
DEG_UP_COLOR = "#FF0000"
DEG_DOWN_COLOR = "#4169E1"
DEG_LINE_COLOR = "#C0392B"

# Overview (multi-contig) layout
OVERVIEW_MAX_CONTIGS = 60
OVERVIEW_CONTIG_COLOR = "#E4E4E4"
OVERVIEW_CGC_COLOR = base_constants.CGC_FEATURE_COLORS["Other"]

# Feature colours: the shared palette
CGC_FEATURE_LEGEND = base_constants.CGC_FEATURE_LEGEND
CGC_FEATURE_COLORS = base_constants.CGC_FEATURE_COLORS
