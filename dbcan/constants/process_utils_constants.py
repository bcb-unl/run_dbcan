
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