import dbcan.constants.base_constants as base_constants

CGC_RESULT_FILE = base_constants.CGC_RESULT_FILE
CGC_SUB_PREDICTION_FILE = base_constants.CGC_SUB_PREDICTION_FILE
DBCAN_SUB_RESULT_FILE = base_constants.DBCAN_SUB_HMM_RESULT_FILE
OVERVIEW_FILE = base_constants.OVERVIEW_FILE
CGC_CAZYME_TYPE = base_constants.CAZYME   # cgc_standard_out.tsv Gene Type values
CGC_NULL_TYPE = base_constants.NULL

# ---------------------------------------------------------------------------
# dbcan_utils abundance outputs (written into -o/--output, default: current directory)
# ---------------------------------------------------------------------------
FAM_ABUND_FILE = "fam_abund.out"
SUBFAM_ABUND_FILE = "subfam_abund.out"
EC_ABUND_FILE = "EC_abund.out"
FAM_SUBSTRATE_ABUND_FILE = "fam_substrate_abund.out"
CGC_ABUND_FILE = "CGC_abund.out"
CGC_SUBSTRATE_HOMOLOGY_FILE = "CGC_substrate_PUL_homology.out"
CGC_SUBSTRATE_MAJORITY_FILE = "CGC_substrate_majority_voting.out"

ABUND_FUNCTIONS = ("fam_abund", "fam_substrate_abund", "CGC_abund", "CGC_substrate_abund")

# Normalisation methods. RPKM is accepted as an alias of FPKM (same formula).
ABUNDANCE_METHODS = ("FPKM", "RPKM", "RPM", "TPM")
DEFAULT_ABUNDANCE_METHOD = "RPM"

# ---------------------------------------------------------------------------
# cal_coverage (read counting)
# ---------------------------------------------------------------------------
COUNT_MODE_FRAGMENTS = "fragments"
COUNT_MODE_READS = "reads"
COUNT_MODES = (COUNT_MODE_FRAGMENTS, COUNT_MODE_READS)
DEFAULT_COUNT_MODE = COUNT_MODE_FRAGMENTS

STRANDEDNESS_NONE = "none"
STRANDEDNESS_FORWARD = "forward"   # read 1 (or a single-end read) on the gene strand
STRANDEDNESS_REVERSE = "reverse"   # read 1 antisense to the gene (dUTP / TruSeq stranded)
STRANDEDNESS_CHOICES = (STRANDEDNESS_NONE, STRANDEDNESS_FORWARD, STRANDEDNESS_REVERSE)
DEFAULT_STRANDEDNESS = STRANDEDNESS_NONE

DEFAULT_OVERLAP_BASE_RATIO = 0.2
DEFAULT_MAPPING_QUALITY = 30
DEFAULT_IDENTITY = 0.98

# GFF column 9 attributes tried (in this order) for the gene ID when --id-attr is 'auto'.
# protein_id: NCBI GFF CDS rows and run_dbcan's cgc.gff; proteinId: JGI; ID: Prodigal/pyrodigal.
GFF_ID_ATTR_AUTO = "auto"
GFF_ID_ATTR_PRIMARY = (base_constants.GFF_PROTEIN_ID_ATTR_NCBI, base_constants.GFF_JGI_PROTEIN_ID_ATTR)
GFF_ID_ATTR_FALLBACK = (base_constants.GFF_ID_ATTR_PRODIGAL, base_constants.GFF_NAME_ATTR, "locus_tag")
GFF_CDS_FEATURE = base_constants.GFF_CDS_FEATURE
# Feature types that never describe a gene (whole sequences / comments)
GFF_SKIP_FEATURE_TYPES = frozenset({
    "region", "chromosome", "contig", "scaffold", "supercontig", "source", "databank_entry", "remark",
})
# Protein IDs run_dbcan writes for non-coding / unresolved features (cgc.gff, cgc_standard_out.tsv)
PLACEHOLDER_GENE_IDS = frozenset({
    base_constants.GFF_NA_PROTEIN_ID, base_constants.GFF_UNKNOWN_ANNOTATION, "", ".", "-",
})
