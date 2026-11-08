import dbcan.constants.base_constants as base_constants
import dbcan.constants.plots_constants as plots_constants

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
# Several contigs (forced circos): contig names inside the innermost ring, radial
CGC_CONTIG_LABEL_R = 62          # radius of the centre of the (radial) contig names, without DEG rings
CGC_CONTIG_LABEL_R_DEG = 38      # same with the DEG rings
CGC_CONTIG_LABEL_SIZE = 6.5
CGC_CONTIG_LABEL_MIN_DEG = 2.5   # narrower sectors are not named (their name would overlap)

CGC_TRACK_BG_COLOR = "#F2F2F2"
CGC_RANGE_COLOR = plots_constants.CGC_EXTENT_COLOR
CGC_RANGE_PREDICTED_COLOR = plots_constants.CGC_EXTENT_PREDICTED_COLOR
CGC_LABEL_SIZE = 7
CGC_TICK_LABEL_SIZE = 7
CGC_FIGURE_SIZE = 9.0            # minimum figure width (inches): room for the legend
CGC_CIRCLE_DIAMETER_IN = 6.2     # diameter of the outer ring (r = 100), inches
CGC_HEADER_IN = 0.7              # title + subtitle
CGC_MARGIN_IN = 0.25             # around the circle and its labels
CGC_LEGEND_ROW_IN = 0.24         # one legend row
DEG_UP_COLOR = "#FF0000"
DEG_DOWN_COLOR = "#4169E1"
DEG_LINE_COLOR = "#C0392B"

# Overview (multi-contig) layout
OVERVIEW_MAX_CONTIGS = 60
OVERVIEW_CONTIG_COLOR = "#ECECEC"
OVERVIEW_CGC_COLOR = plots_constants.CGC_EXTENT_COLOR
OVERVIEW_CGC_PREDICTED_COLOR = plots_constants.CGC_EXTENT_PREDICTED_COLOR

# Feature colours: the shared palette
CGC_FEATURE_LEGEND = base_constants.CGC_FEATURE_LEGEND
CGC_FEATURE_COLORS = base_constants.CGC_FEATURE_COLORS
CGC_CIRCOS_MAX_ID_LABELS = 60   # more CGCs than this: only CGCs with a substrate are labelled
