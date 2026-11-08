
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
ANNOTATE_NAME_COLUMN = "Annotate Name"
ANNOTATE_FROM_COLUMN = "Annotate From"
ANNOTATE_TO_COLUMN = "Annotate To"


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
# A domain-level signature on a CAZyme domain (overlap > this fraction of the shorter
# interval) describes the same region and is not an extra CGC signature.
SAME_REGION_OVERLAP_THRESHOLD = 0.5
