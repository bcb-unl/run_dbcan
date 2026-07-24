
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
# Optional extra column, emitted only when the opt-in matched-region overlap mode
# is active. It carries the insertion-free matched target intervals of a domain as
# a comma-separated "start:end" string (1-based inclusive), e.g. "21:252,340:414".
MATCHED_INTERVALS_COLUMN = "Matched Intervals"
ANNOTATE_FROM_COLUMN = "Annotate From"
ANNOTATE_TO_COLUMN = "Annotate To"


CGC_SIG_RESULT_COLUMN = [
            'Annotate Name', 'Annotate Length', 'Target Name', 'Target Length',
            'i-Evalue', 'Annotate From', 'Annotate To', 'Target From', 'Target To',
            'Coverage', 'Annotate File Name'
        ]