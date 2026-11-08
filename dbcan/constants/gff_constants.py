import dbcan.constants.base_constants as base_constants


# Input/Output file names
GFF_INPUT_PROTEIN_FILE = base_constants.INPUT_PROTEIN_NAME
GFF_CAZYME_OVERVIEW_FILE = base_constants.OVERVIEW_FILE
GFF_CGC_SIG_FILE = base_constants.GFF_CGC_SIG_FILE 
GFF_OUTPUT_FILE = base_constants.CGC_GFF_FILE
GFF_TEMP_SUFFIX = base_constants.GFF_TEMP_SUFFIX

# Column names and indices
GFF_PROTEIN_ID_COL = base_constants.GFF_PROTEIN_ID_COL
GFF_CAZYME_COL = base_constants.GFF_CAZYME_COL
GFF_GENE_ID_COL = base_constants.GFF_GENE_ID_COL
GFF_TOOLS_COUNT_COL = base_constants.GFF_TOOLS_COUNT_COL
GFF_RECOMMEND_RESULTS_COL = base_constants.GFF_RECOMMEND_RESULTS_COL
GFF_CGC_ANNOTATION_COL = base_constants.GFF_CGC_ANNOTATION_COL
GFF_FUNCTION_ANNOTATION_COL = base_constants.GFF_FUNCTION_ANNOTATION_COL
GFF_TYPE_COL = base_constants.GFF_TYPE_COL
GFF_CGC_SIG_COLUMNS = base_constants.GFF_CGC_SIG_COLUMNS
GFF_MIN_TOOL_COUNT = base_constants.GFF_MIN_TOOL_COUNT



# Annotation prefixes and defaults
GFF_CAZYME_PREFIX = base_constants.GFF_CAZYME_PREFIX
GFF_OTHER_PREFIX = base_constants.GFF_OTHER_PREFIX
GFF_NULL_ANNOTATION = base_constants.GFF_NULL_ANNOTATION
GFF_UNKNOWN_ANNOTATION = base_constants.GFF_UNKNOWN_ANNOTATION
GFF_NA_PROTEIN_ID = base_constants.GFF_NA_PROTEIN_ID

# GFF feature types
GFF_GENE_FEATURE = base_constants.GFF_GENE_FEATURE
GFF_MRNA_FEATURE = base_constants.GFF_MRNA_FEATURE
GFF_CDS_FEATURE = base_constants.GFF_CDS_FEATURE

# GFF format types
GFF_FORMAT_NCBI_EUK = base_constants.GFF_FORMAT_NCBI_EUK
GFF_FORMAT_NCBI_PROK = base_constants.GFF_FORMAT_NCBI_PROK
GFF_FORMAT_JGI = base_constants.GFF_FORMAT_JGI
GFF_FORMAT_PRODIGAL = base_constants.GFF_FORMAT_PRODIGAL
GFF_FORMAT_FUNANNOTATE = base_constants.GFF_FORMAT_FUNANNOTATE
GFF_FORMAT_BRAKER = base_constants.GFF_FORMAT_BRAKER
GFF_FORMAT_GFF3 = base_constants.GFF_FORMAT_GFF3
GFF_FORMAT_GTF = base_constants.GFF_FORMAT_GTF
GFF_LEGACY_FORMATS = base_constants.GFF_LEGACY_FORMATS
GFF_RESOLVER_FORMATS = base_constants.GFF_RESOLVER_FORMATS
GFF_FORMAT_CHOICES = base_constants.GFF_FORMAT_CHOICES

# gene -> protein resolution report
GFF_RESOLUTION_FILE = base_constants.GFF_RESOLUTION_FILE
GFF_MIN_MATCH_DEFAULT = base_constants.GFF_MIN_MATCH_DEFAULT
GFF_STATUS_MATCHED = "matched"
GFF_STATUS_UNMATCHED = "unmatched"
GFF_STATUS_NONCODING = "noncoding"
GFF_RESOLUTION_COLUMNS = [
    "gene_id", "contig", "start", "end", "strand", "protein_id", "matched_by", "status",
    "n_transcripts", "isoform_proteins",
]

# GFF3 feature types used by the resolver
GFF_RESOLVER_CDS_TYPES = ("CDS", "cds")
GFF_RESOLVER_GENE_TYPES = (
    "gene", "protein_coding_gene", "pseudogene", "ncRNA_gene", "tRNA_gene", "rRNA_gene",
)
GFF_RESOLVER_TRANSCRIPT_TYPES = (
    "mRNA", "transcript", "tRNA", "rRNA", "ncRNA", "lnc_RNA", "lncRNA", "snRNA", "snoRNA",
    "pseudogenic_transcript", "primary_transcript", "miRNA", "pre_miRNA",
)
# Pseudogenes: kept as non-coding null rows (Other|pseudogene) unless their CDS names a protein
GFF_RESOLVER_PSEUDO_TYPES = ("pseudogene", "pseudogenic_transcript")
GFF_PSEUDOGENE_FEATURE = "pseudogene"
# gffread-style transcripts without a gene parent: the gene is named by these attributes
GFF_RESOLVER_GENE_KEY_ATTRS = ("geneID", "gene_id", "gene")
# Identifiers matched first (exact FASTA ID), before Name/locus_tag/gene-level keys and aliases
GFF_RESOLVER_STRONG_KEY_ATTRS = (
    "protein_id", "proteinId", "ID", "transcript_id", "orig_protein_id", "orig_transcript_id",
)
# Rows below the transcript level that never identify a gene or protein; not kept in memory.
GFF_RESOLVER_IGNORED_TYPES = (
    "exon", "intron", "five_prime_UTR", "three_prime_UTR", "5'-UTR", "3'-UTR", "UTR",
    "start_codon", "stop_codon", "region", "chromosome", "contig", "scaffold", "supercontig",
    "remark", "source", "tss", "tts", "transcription_start_site", "transcription_end_site",
)
# Prefixes / suffixes that annotation pipelines add to the protein-bearing identifier
# (NCBI 'cds-'/'rna-', Ensembl 'CDS:'/'transcript:', funannotate '.cds', BRAKER '.CDS1', gffread '-P1' ...).
GFF_RESOLVER_ID_PREFIXES = (
    "cds-", "cds:", "CDS:", "rna-", "transcript:", "mRNA:", "gene-", "gene:", "mrna-",
)
GFF_RESOLVER_ID_SUFFIX_PATTERNS = (
    r"\.cds\d*$", r"-CDS\d*$", r"\.CDS\.?\d*$", r":cds$", r"_cds$", r"-P\d*$",
)

# GFF attribute names
GFF_PROTEIN_ID_ATTR_NCBI = base_constants.GFF_PROTEIN_ID_ATTR_NCBI
GFF_NAME_ATTR = base_constants.GFF_NAME_ATTR
GFF_ID_ATTR_PRODIGAL = base_constants.GFF_ID_ATTR_PRODIGAL
GFF_JGI_PROTEIN_ID_ATTR = base_constants.GFF_JGI_PROTEIN_ID_ATTR

################################################################################################
