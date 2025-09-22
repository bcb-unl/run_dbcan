from pathlib import Path

# root path of the package
TEST_ROOT = Path(__file__).parent

# common input files
INPUT_PROTEIN_NAME = "uniInput.faa"
INPUT_PRODIGAL_GFF_NAME = "uniInput.gff"
NULL_PROTEIN_FILE = "null_proteins.faa"

# common output files (shared across modules)
OVERVIEW_FILE = "overview.tsv"
CGC_GFF_FILE = "cgc.gff"
CGC_RESULT_FILE = "cgc_standard_out.tsv"


# module specific output files
DBCAN_HMM_RESULT_FILE = "dbCAN_hmm_results.tsv"
DBCAN_SUB_HMM_RAW_FILE = "dbCANsub_hmm_raw.tsv"
DBCAN_SUB_HMM_RESULT_FILE = "dbCANsub_hmm_results.tsv"
TF_HMM_RESULT_FILE = "TF_hmm_results.tsv"
STP_HMM_RESULT_FILE = "STP_hmm_results.tsv"
PFAM_HMM_RESULT_FILE = "Pfam_hmm_results.tsv"
CAZY_DIAMOND_OUTPUT = "diamond.out"
TCDB_DIAMOND_OUTPUT = "diamond.out.tc"
SULFATLAS_DIAMOND_OUTPUT = "diamond.out.sulfatase"
PEPTIDASE_DIAMOND_OUTPUT = "diamond.out.peptidase"
TF_DIAMOND_OUTPUT = "diamond.out.tf"



#substrate output
CGC_SUB_PREDICTION_FILE = "substrate_prediction.tsv"
PUL_DIAMOND_FILE = "PUL_blast.out"
CGC_FAA_FILE = "CGC.faa"
PUL_DIAMOND_DB = "PUL.dmnd"
PUL_EXCEL_FILE = "dbCAN-PUL.xlsx"
CAZYME_FAA_FILE = "CAZyme.faa"
PUL_FAA_FILE = "PUL.faa"


# gene annotation related constants
CONTIG_ID_COLUMN = "Contig ID"
PROTEIN_ID_COLUMN = "Protein_ID"
CGC_ANNOTATION_COLUMN = "CGC_annotation"
START_COLUMN = "start"
END_COLUMN = "end"
STRAND_COLUMN = "strand"
ATTRIBUTES_COLUMN = "attributes"
GFF_COLUMNS = ["Contig ID", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]




# general CGC fields
CGC_ID_FIELD = "CGC#"
CGC_PROTEIN_ID_FIELD = "Protein ID"
GENE_TYPE_FIELD = "Gene Type"
GENE_START_FIELD = "Gene Start"
GENE_STOP_FIELD = "Gene Stop"
GENE_STRAND_FIELD = "Gene Strand"
GENE_ANNOTATION_FIELD = "Gene Annotation"

#plot for CGCs
CGC_FEATURE_LEGEND = ['CAZyme', 'TC', 'TF', 'STP', 'Peptidase', 'Sulfatase']
CGC_FEATURE_COLORS = {
    "CAZyme": "#E67E22",      # orange
    "TC": "#2ECC71",          # green
    "TF": "#9B59B6",          # purple
    "STP": "#F1C40F",         # golden yellow
    "Peptidase": "#16A085",   # greenish
    "Sulfatase": "#010E1B",   # dark blue
    "Other": "#95A5A6"      # light gray
}
GENE_LABEL_COLOR = list(CGC_FEATURE_COLORS.values())
GENE_LABELS = list(CGC_FEATURE_COLORS.keys())

# common label for different feature types
CAZYME = "CAZyme"
TC = "TC"
TF = "TF"
STP = "STP"
PUL = "PUL"
NULL = "null"
PEPTIDASE = "Peptidase"
SULFATASE = "Sulfatase"



# Input/Output file names

GFF_CGC_SIG_FILE = "total_cgc_info.tsv"
GFF_TEMP_SUFFIX = ".temp"

# Column names and indices
GFF_PROTEIN_ID_COL = "protein_id"
GFF_CAZYME_COL = "CAZyme"
GFF_GENE_ID_COL = "Gene ID"
GFF_TOOLS_COUNT_COL = "#ofTools"
GFF_RECOMMEND_RESULTS_COL = "Recommend Results"
GFF_CGC_ANNOTATION_COL = 'CGC_annotation'
GFF_FUNCTION_ANNOTATION_COL = "function_annotation"
GFF_TYPE_COL = "type"
GFF_CGC_SIG_COLUMNS = [0, 2, 10]
GFF_MIN_TOOL_COUNT = 2

# Annotation prefixes and defaults
GFF_CAZYME_PREFIX = 'CAZyme|'
GFF_OTHER_PREFIX = "Other|"
GFF_NULL_ANNOTATION = 'null'
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
GFF_PROTEIN_ID_ATTR_NCBI = "protein_id"
GFF_NAME_ATTR = "Name"
GFF_ID_ATTR_PRODIGAL = "ID"
GFF_JGI_PROTEIN_ID_ATTR = "proteinId"

################################################################################################




__all__ = [name for name in globals() if name.isupper()]
