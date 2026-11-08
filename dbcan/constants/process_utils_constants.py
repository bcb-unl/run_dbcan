
#########constants for process_utils.py##############################
OVERLAP_RATIO_THRESHOLD = 0.5

HMMER_COLUMN_NAMES = [
        'HMM Name', 'HMM Length', 'Target Name', 'Target Length', 'i-Evalue',
        'HMM From', 'HMM To', 'Target From', 'Target To', 'Coverage', 'HMM File Name'
    ]

# Column names for input hmm data
HMM_NAME_COLUMN = "HMM Name"
TARGET_NAME_COLUMN = "Target Name"
TARGET_LENGTH_COLUMN = "Target Length"
IEVALUE_COLUMN = "i-Evalue"
HMM_LENGTH_COLUMN = "HMM Length"
HMM_FROM_COLUMN = "HMM From"
HMM_TO_COLUMN = "HMM To"
TARGET_FROM_COLUMN = "Target From"
TARGET_TO_COLUMN = "Target To"
COVERAGE_COLUMN = "Coverage"
HMM_FILE_COLUMN = "HMM File Name"
# Internal-only column, produced when matched-region overlap resolution is active.
# It carries the match-state target intervals of a domain as a comma-separated
# "start:end" string (1-based inclusive), e.g. "21:252,340:414". It is used for
# domain resolution and mirrored into ALIGNMENT_DETAILS_COLUMNS, but is never part
# of a public result table: HMMER_COLUMN_NAMES is the public schema in every mode.
MATCHED_INTERVALS_COLUMN = "Matched Intervals"
ANNOTATE_NAME_COLUMN = "Annotate Name"
ANNOTATE_LENGTH_COLUMN = "Annotate Length"
ANNOTATE_FROM_COLUMN = "Annotate From"
ANNOTATE_TO_COLUMN = "Annotate To"

# Optional diagnostics sidecar written next to each HMM result table
# (e.g. dbCAN_hmm_results.tsv -> dbCAN_hmm_alignment_details.tsv). It keeps the
# resolved domain identity — alignment span *and* match-state intervals — so
# downstream steps can consume resolved domains instead of re-inferring the
# architecture from the continuous span.
ALIGNMENT_DETAILS_COLUMNS = [
    HMM_NAME_COLUMN,
    TARGET_NAME_COLUMN,
    TARGET_LENGTH_COLUMN,
    TARGET_FROM_COLUMN,
    TARGET_TO_COLUMN,
    MATCHED_INTERVALS_COLUMN,
    IEVALUE_COLUMN,
    COVERAGE_COLUMN,
    HMM_FILE_COLUMN,
]


CGC_SIG_RESULT_COLUMN = [
            'Annotate Name', 'Annotate Length', 'Target Name', 'Target Length',
            'i-Evalue', 'Annotate From', 'Annotate To', 'Target From', 'Target To',
            'Coverage', 'Annotate File Name'
        ]
# CGC signature sources (process_cgc_sig_results). DIAMOND sources search full-length
# protein databases and give protein-level labels; HMM sources give Pfam-style domains.
# total_cgc_info.tsv column naming the search a signature hit came from
SIGNATURE_SOURCE_COLUMN = 'Source'
DIAMOND_SIGNATURE_SOURCES = frozenset({'TC', 'TF_prok', 'Sulfatase', 'Peptidase'})
HMM_SIGNATURE_SOURCES = frozenset({'STP', 'TF_fungi'})

# MEROPS identifiers in peptidase_db.dmnd (pepunit.lib): '<family>.<id>', e.g. S09.940,
# C56.UNW. Families I* are peptidase inhibitors, not peptidases.
MEROPS_PEPTIDASE_ID_PATTERN = r"^[A-HJ-Z]\d+[A-Z]?\.[A-Za-z0-9]+$"
MEROPS_INHIBITOR_ID_PATTERN = r"^I\d+[A-Z]?\.[A-Za-z0-9]+$"

# Signature sources whose hits are domains (coordinates comparable with dbCAN domains):
# STP / TF HMM domains and MEROPS peptidase units (pepunit.lib).
DOMAIN_LEVEL_SIGNATURE_SOURCES = frozenset({'STP', 'TF_fungi', 'Peptidase'})
# Full-length protein databases (TCDB, prokaryotic TF DIAMOND, SulfAtlas): a hit labels the whole
# protein and its position cannot be compared with a domain. On a protein that is a CGC CAZyme these
# labels are dropped from cgc.gff (they stay in total_cgc_info.tsv): the CAZyme annotation wins.
FULL_LENGTH_SIGNATURE_SOURCES = frozenset({'TC', 'TF_prok', 'Sulfatase'})
# On a protein that is NOT a CGC CAZyme, only the best of its full-length labels is written to
# cgc.gff (the others stay in total_cgc_info.tsv). Best = highest DIAMOND bit score; ties go to
# the higher min(query coverage, subject coverage), then the lower E-value, then this order.
FULL_LENGTH_SOURCE_PRIORITY = ('TC', 'TF_prok', 'Sulfatase')
# total_cgc_info.tsv column (last) with the DIAMOND bit score of a hit; empty for HMM sources.
SIGNATURE_BITSCORE_COLUMN = 'Bitscore'
# A domain-level signature on a CAZyme domain (overlap > this fraction of the shorter
# interval) describes the same region and is not an extra CGC signature.
SAME_REGION_OVERLAP_THRESHOLD = 0.5
