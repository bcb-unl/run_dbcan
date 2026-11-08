"""Constants for rebuilding peptidase_db.dmnd from MEROPS pepunit.lib (maintainer tool).

pepunit.lib (https://ftp.ebi.ac.uk/pub/databases/merops/current_release/pepunit.lib) holds
the peptidase units *and* the inhibitor units of MEROPS, one FASTA record per unit:

    >MER0000001 - chymotrypsin B (Homo sapiens) [S01.152]#S01A#{peptidase unit: 34-263}~source CTRB_HUMAN~
    >MER0004386 - alpha-1-antitrypsin (Homo sapiens) [I04.001]#I04#{inhibitor unit: 25-418}~source A1AT_HUMAN~

The MEROPS identifier is the *last* bracketed token before ``#<family>#``; names and
organisms can contain brackets of their own (``[misleading]``, ``[Ser->Ala replacement]``,
``([Clostridium] irregulare)``), so the first bracket is not the identifier.
"""
import dbcan.constants.diamond_constants as D
import dbcan.constants.process_utils_constants as PU

# Release files
MEROPS_RELEASE_URL = "https://ftp.ebi.ac.uk/pub/databases/merops/current_release"
MEROPS_PEPUNIT_FILE = "pepunit.lib"
MEROPS_PEPUNIT_URL = f"{MEROPS_RELEASE_URL}/{MEROPS_PEPUNIT_FILE}"

# Header parsing
MEROPS_ACCESSION_PATTERN = r"^MER\d+$"
# '[<MEROPS id>]#<family/subfamily>#{<peptidase|inhibitor> unit: <ranges>}'; the last match
# in the header is used.
MEROPS_ID_FIELD_PATTERN = r"\[([^\[\]]+)\]#([^#]*)#\{\s*([A-Za-z]+) unit\s*:\s*([^}]*)\}"
# Fallback when the '#...#{... unit' tail is missing: the last bracketed token.
MEROPS_BRACKET_PATTERN = r"\[([^\[\]]+)\]"
MEROPS_SOURCE_PATTERN = r"~source\s+([^~]*)~"
UNIT_TYPE_PEPTIDASE = "peptidase"
UNIT_TYPE_INHIBITOR = "inhibitor"

# Same identifier rules as the Peptidase signature step (process_utils._drop_non_peptidase_units).
MEROPS_PEPTIDASE_ID_PATTERN = PU.MEROPS_PEPTIDASE_ID_PATTERN
MEROPS_INHIBITOR_ID_PATTERN = PU.MEROPS_INHIBITOR_ID_PATTERN

# Output FASTA ID: '<MER accession>|<MEROPS id>'. DIAMOND reports it as sseqid and
# PeptidaseDiamondProcessor keeps the part after '|' as 'Peptidase ID'.
OUTPUT_ID_SEPARATOR = "|"

# Reasons an entry is left out (keys of the summary and of the --report table).
DROP_INHIBITOR = "inhibitor_unit"
DROP_INVALID_ID = "invalid_merops_id"
DROP_INVALID_ACCESSION = "invalid_accession"
DROP_EMPTY_SEQUENCE = "empty_sequence"
DROP_DUPLICATE = "duplicate_id_and_sequence"
DROP_DUPLICATE_ACCESSION = "duplicate_accession"
DROP_REASONS = (
    DROP_INHIBITOR,
    DROP_INVALID_ID,
    DROP_INVALID_ACCESSION,
    DROP_EMPTY_SEQUENCE,
    DROP_DUPLICATE,
    DROP_DUPLICATE_ACCESSION,
)
REPORT_COLUMNS = ("accession", "merops_id", "reason", "header")

# Defaults
DEFAULT_OUTPUT_DB = D.PEPTIDASE_DIAMOND_DB
FASTA_SUFFIX = ".faa"
DIAMOND_CMD = D.DIAMOND_CMD
DIAMOND_MAKEDB = "makedb"
FASTA_LINE_WIDTH = 60
