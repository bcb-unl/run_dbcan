import dbcan.constants.base_constants as base_constants


CGC_RESULT_FILE = base_constants.CGC_RESULT_FILE
CGC_SUB_PREDICTION_FILE = base_constants.CGC_SUB_PREDICTION_FILE
CGC_GFF_FILE = base_constants.CGC_GFF_FILE
PUL_DIAMOND_FILE = base_constants.PUL_DIAMOND_FILE

GFF_PROTEIN_ID_COL = base_constants.GFF_PROTEIN_ID_COL
GFF_CGC_ANNOTATION_COL = base_constants.GFF_CGC_ANNOTATION_COL
CGC_ID_FIELD = base_constants.CGC_ID_FIELD

# ---------------------------------------------------------------------------
# Gene-type colours: single source in base_constants (see CGC_FEATURE_COLORS there).
# ---------------------------------------------------------------------------
CGC_FEATURE_LEGEND = base_constants.CGC_FEATURE_LEGEND
CGC_FEATURE_COLORS = base_constants.CGC_FEATURE_COLORS
CGC_GENE_TYPE_PRIORITY = base_constants.CGC_GENE_TYPE_PRIORITY
CGC_GENE_TYPE_ALIASES = base_constants.CGC_GENE_TYPE_ALIASES
OTHER_GENE_TYPE = base_constants.CGC_OTHER_GENE_TYPE   # palette key for genes without a signature hit (null) and RNAs

# ---------------------------------------------------------------------------
# Shared look of the linear figures (CGC map, coverage, synteny, overview)
# ---------------------------------------------------------------------------
PLOT_FONT_FAMILY = "DejaVu Sans"
PLOT_BASE_FONT_SIZE = 8
PLOT_TITLE_FONT_SIZE = 10.5
PLOT_SUBTITLE_FONT_SIZE = 8.5
PLOT_LABEL_FONT_SIZE = 7.5
PLOT_SMALL_FONT_SIZE = 6.8
PLOT_INK = "#1F1F1F"            # main text
PLOT_INK_MUTED = "#555555"      # secondary text
PLOT_BACKBONE_COLOR = "#9A9A9A"  # line joining the genes of a track
PLOT_LEADER_COLOR = "#B5B5B5"   # leader line from a lifted label to its gene
PLOT_GENE_GAP_COLOR = "white"   # thin outline that separates abutting arrows
PLOT_DEFAULT_DPI = 200
PLOT_DEFAULT_FORMATS = ("pdf",)
PLOT_SUPPORTED_FORMATS = ("pdf", "png", "svg")
# matplotlib backends that cannot open a window (dbcan_plot heatmap_plot --show-fig)
PLOT_NON_INTERACTIVE_BACKENDS = ("agg", "pdf", "ps", "svg", "pgf", "cairo", "template")

# CGC extent blocks (circos CGC ring, multi-contig overview). Deliberately NOT gene-type
# colours, so that "CGC with substrate prediction" cannot be mistaken for a CAZyme gene.
CGC_EXTENT_COLOR = "#ABB2B9"            # CGC without a substrate prediction (cool grey)
CGC_EXTENT_PREDICTED_COLOR = "#34495E"  # CGC with a substrate prediction (slate)

# Layout geometry of the linear figures, in inches unless noted ("U" = inches per y data unit)
CGC_MAP_UNIT_IN = 0.27          # inches per y data unit of the gene axis
CGC_MAP_GENE_H = 1.0            # gene arrow height (data units)
CGC_MAP_LEFT_IN = 0.45
CGC_MAP_RIGHT_IN = 0.4
CGC_MAP_TITLE_IN = 0.62         # title + substrate line
CGC_MAP_LEGEND_IN = 0.42
CGC_MAP_DEPTH_PANEL_IN = 0.85   # height of one read-depth panel
CGC_MAP_DEPTH_GAP_IN = 0.14
CGC_MAP_RULER_ONLY_DEPTH = 0.75  # room under the ruler when the gene track has no scale bar (data units)
SYNTENY_UNIT_IN = 0.27
SYNTENY_GENE_H = 1.0
SYNTENY_CGC_TRACK_Y = 3.2       # y of the CGC track (data units); the PUL track is at y = 0
SYNTENY_TITLE_IN = 0.62
SYNTENY_LEGEND_IN = 0.42
SYNTENY_RIGHT_IN = 0.35
OVERVIEW_ROW_IN = 0.46          # height of one contig row
OVERVIEW_BAR_IN = 0.09          # contig bar thickness
OVERVIEW_CGC_IN = 0.2           # CGC block height
OVERVIEW_TITLE_IN = 0.62
OVERVIEW_FOOT_IN = 0.5
OVERVIEW_FIG_W = 11.0
OVERVIEW_RIGHT_IN = 0.75        # room for the length label after the longest bar

# read-depth panels (dbcan_plot CGC_coverage_plot / CGC_synteny_coverage_plot);
# the original dbCAN coverage colours (red line over a light red fill)
COVERAGE_LINE_COLOR = "red"
COVERAGE_LINE_ALPHA = 0.6
COVERAGE_FILL_COLOR = "red"
COVERAGE_FILL_ALPHA = 0.3
COVERAGE_GENE_SPAN_COLOR = "#F2F2F2"
COVERAGE_YLABEL = "read depth"

# synteny links: % identity on a single-hue grey ramp
SYNTENY_IDENTITY_RAMP = ("#E6E6E6", "#BDBDBD", "#6E6E6E")
SYNTENY_IDENTITY_MIN = 20.0
SYNTENY_IDENTITY_MAX = 100.0
# % identity printed on each homology link, between the two tracks
SYNTENY_LINK_LABEL_FORMAT = "{:.0f}%"
SYNTENY_LINK_LABEL_FONT_SIZE = 8.0
# Link ribbons are narrower than the genes so crossing links stay visible; a gene with several
# links gives each one slot of its width
SYNTENY_LINK_WIDTH_FRACTION = 0.35       # ribbon width at a gene, as a fraction of the gene width
SYNTENY_LINK_SLOT_FILL = 0.8             # with several links on one gene: fraction of a slot used
SYNTENY_LINK_EDGE_COLOR = "white"        # hairline that separates crossing links
SYNTENY_LINK_EDGE_LW = 0.5
# where a % tag sits along its ribbon (0 = CGC track, 1 = PUL track); later positions are tried when tags collide
SYNTENY_LINK_LABEL_POSITIONS = (0.5, 0.32, 0.68, 0.2, 0.8)
SYNTENY_PDF_DIR = "synteny_pdf"
SYNTENY_FILE_TEMPLATE = "{cgc}_{pul}-syntenic.{ext}"
DBCAN_PUL_DIR = "dbCAN-PUL"
DBCAN_PUL_GFF = "cgc.gff"

# dbcan_plot output names
CGC_MAP_FILE_TEMPLATE = "{cgc}.cgc.{ext}"
CGC_COVERAGE_FILE_TEMPLATE = "{cgc}.cgc-coverage.{ext}"
CGC_SYNTENY_COVERAGE_FILE_TEMPLATE = "{cgc}_{pul}-syntenic-cov.{ext}"
HEATMAP_FILE = "heatmap.pdf"
HEATMAP_CLUSTER_FILE = "heatmap_cluster.pdf"
BAR_PLOT_FILE = "bar_plot.pdf"

# matplotlib >= 3.6 renamed the seaborn styles; old names are mapped to the new ones
PLOT_STYLE_CHOICES = ("ggplot", "seaborn-v0_8", "seaborn-v0_8-poster", "seaborn-v0_8-whitegrid",
                      "default", "classic", "bmh", "fivethirtyeight")
PLOT_STYLE_ALIASES = {"seaborn": "seaborn-v0_8", "seaborn-poster": "seaborn-v0_8-poster",
                      "seaborn-whitegrid": "seaborn-v0_8-whitegrid"}
