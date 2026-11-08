"""peptidase_db.dmnd builder (dbcan.utils.build_peptidase_db) on a synthetic pepunit.lib.

Headers follow the real MEROPS pepunit.lib format; the MEROPS identifier is the last
bracket before '#<family>#{... unit: ...}', not the first one ('[misleading]',
'([Clostridium] ...)'), which is how the previous peptidase_db.dmnd got IDs such as
'misleading' and 'Clostridium'.
"""
import gzip
import re
import shutil

import pandas as pd
import pytest
from click.testing import CliRunner

import dbcan.constants.peptidase_db_constants as P
import dbcan.constants.process_utils_constants as PU
from dbcan.annotation.diamond import PeptidaseDiamondProcessor
from dbcan.configs.peptidase_db_config import PeptidaseDbBuildConfig
from dbcan.utils.build_peptidase_db import (
    PeptidaseDbBuildError, build_peptidase_db, cli, parse_header, write_clean_fasta,
)

SEQ_A = "IVNGEDAVPGSWPWQVSLQDKTGFHFCGGSLISEDWVVTAAHCGVRTSDVVVAGEFDQGS"
SEQ_B = "MKVLLTGAGGFIGSHLVDRLLAEGHEVTVLDNLSTGRRENLPAGVRFVEGDIRDAAAL"
SEQ_C = "QWERTYIPASDFGHKLCVNMQWERTYIPASDFGHKLCVNM"

RECORDS = [
    # kept
    (">MER0000001 - chymotrypsin B (Homo sapiens) [S01.152]#S01A#{peptidase unit: 34-263}~source CTRB_HUMAN~", SEQ_A),
    # inhibitor unit
    (">MER0004386 - alpha-1-antitrypsin (Homo sapiens) [I04.001]#I04#{inhibitor unit: 25-418}~source A1AT_HUMAN~", SEQ_B),
    # bracket in the name: ID is C19.972, not 'misleading'
    (">MER0003996 - ubiquitin-specific endopeptidase 39 [misleading] (Mus musculus) [C19.972]#C19#{peptidase unit: 187-549}~source EDK98966~", SEQ_C),
    # bracket in the organism: ID is M15.UPB, not 'Clostridium'
    (">MER0333466 - subfamily M15B unassigned peptidases ([Clostridium] lituseburense) [M15.UPB]#M15B#{peptidase unit: 33-222}~source WP_092727540~", SEQ_B),
    # two brackets before the ID
    (">MER0037280 - EGF-like module containing mucin-like hormone receptor-like 1 ({Mus musculus}) [Ser->Ala replacement] (Mus musculus) [P02.951]#P02A#{peptidase unit: 591-637}~source EMR1_MOUSE~", SEQ_A[:30]),
    # malformed identifier (lower-case family)
    (">MER0697988 - subfamily A1A non-peptidase homologues (Gossypium raimondii) [a01.UNa]#A01A#{peptidase unit: 100-442}~source A0A0D2QDX5~", SEQ_C),
    # exact duplicate of MER0000001 (same MEROPS id, same sequence)
    (">MER0000003 - chymotrypsin B (Rattus norvegicus) [S01.152]#S01A#{peptidase unit: 34-263}~source CTRB_RAT~", SEQ_A.lower()),
    # same MEROPS id, different sequence: kept
    (">MER0000004 - chymotrypsin B (Bos taurus) [S01.152]#S01A#{peptidase unit: 16-245}~source CTRB_BOVIN~", SEQ_B),
    # same sequence, different MEROPS id: kept
    (">MER0000006 - chymotrypsin B ({Gadus}-type) (Gadus morhua) [S01.437]#S01A#{peptidase unit: 34-263}~source P47796~", SEQ_A),
    # no identifier at all
    (">MER0099999 - broken header without identifier", SEQ_B),
    # not a MEROPS accession
    (">XYZ0000001 - something (Homo sapiens) [S01.001]#S01A#{peptidase unit: 1-10}~source X~", SEQ_B),
    # empty sequence
    (">MER0000010 - empty (Homo sapiens) [S01.001]#S01A#{peptidase unit: 1-10}~source X~", ""),
]

EXPECTED_IDS = [
    "MER0000001|S01.152", "MER0003996|C19.972", "MER0333466|M15.UPB",
    "MER0037280|P02.951", "MER0000004|S01.152", "MER0000006|S01.437",
]


def _write_lib(path, records=RECORDS, gz=False):
    text = "".join(f"{h}\n{s[:25]}\n{s[25:]}\n" if s else f"{h}\n" for h, s in records)
    if gz:
        with gzip.open(path, "wt") as fh:
            fh.write(text)
    else:
        path.write_text(text)
    return path


def _read_fasta(path):
    ids, seqs = [], []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            ids.append(line[1:])
            seqs.append("")
        else:
            seqs[-1] += line
    return ids, seqs


@pytest.mark.parametrize("header, merops_id, unit", [
    (RECORDS[0][0], "S01.152", "peptidase"),
    (RECORDS[1][0], "I04.001", "inhibitor"),
    (RECORDS[2][0], "C19.972", "peptidase"),
    (RECORDS[3][0], "M15.UPB", "peptidase"),
    (RECORDS[4][0], "P02.951", "peptidase"),
    (">MER0000002 - name (Homo sapiens) [S01.001] trailing text without family tail", "S01.001", None),
    (RECORDS[9][0], None, None),
])
def test_parse_header_takes_identifier_before_family_field(header, merops_id, unit):
    entry = parse_header(header)
    assert entry.merops_id == merops_id
    assert entry.unit_type == unit
    assert entry.accession == header[1:].split()[0]


def test_clean_fasta_keeps_peptidase_units_only(tmp_path):
    lib = _write_lib(tmp_path / "pepunit.lib")
    out = tmp_path / "pep.faa"
    report = tmp_path / "dropped.tsv"
    summary = write_clean_fasta(lib, out, report)

    ids, seqs = _read_fasta(out)
    assert ids == EXPECTED_IDS
    assert seqs[0] == SEQ_A and all(s == s.upper() for s in seqs)
    assert summary.total == len(RECORDS) and summary.kept == len(EXPECTED_IDS)
    assert dict(summary.dropped) == {
        P.DROP_INHIBITOR: 1, P.DROP_INVALID_ID: 2, P.DROP_DUPLICATE: 1,
        P.DROP_INVALID_ACCESSION: 1, P.DROP_EMPTY_SEQUENCE: 1,
    }
    assert summary.families == {"S01": 3, "C19": 1, "M15": 1, "P02": 1}

    rep = pd.read_csv(report, sep="\t", dtype=str, keep_default_na=False)
    assert list(rep.columns) == list(P.REPORT_COLUMNS)
    assert dict(zip(rep["accession"], rep["reason"])) == {
        "MER0004386": P.DROP_INHIBITOR, "MER0697988": P.DROP_INVALID_ID, "MER0000003": P.DROP_DUPLICATE,
        "MER0099999": P.DROP_INVALID_ID, "XYZ0000001": P.DROP_INVALID_ACCESSION,
        "MER0000010": P.DROP_EMPTY_SEQUENCE,
    }


def test_output_ids_are_what_the_peptidase_step_reports(tmp_path):
    """sseqid '<MER>|<id>' -> 'Peptidase ID' <id>, accepted by the Peptidase signature filter."""
    out = tmp_path / "pep.faa"
    write_clean_fasta(_write_lib(tmp_path / "pepunit.lib"), out)
    ids, _ = _read_fasta(out)
    df = pd.DataFrame({"Peptidase ID": ids})
    PeptidaseDiamondProcessor._postprocess_ids(None, df, "Peptidase ID")
    assert list(df["Peptidase ID"]) == [i.split("|")[1] for i in EXPECTED_IDS]
    assert df["Peptidase ID"].str.match(PU.MEROPS_PEPTIDASE_ID_PATTERN).all()


def test_keep_inhibitors_as_decoys(tmp_path):
    out = tmp_path / "pep.faa"
    summary = write_clean_fasta(_write_lib(tmp_path / "pepunit.lib"), out, keep_inhibitors=True)
    ids, _ = _read_fasta(out)
    assert "MER0004386|I04.001" in ids
    assert summary.inhibitors_kept == 1 and P.DROP_INHIBITOR not in summary.dropped
    inh = [i.split("|")[1] for i in ids if re.match(PU.MEROPS_INHIBITOR_ID_PATTERN, i.split("|")[1])]
    assert inh == ["I04.001"]


def test_gzip_input(tmp_path):
    out = tmp_path / "pep.faa"
    write_clean_fasta(_write_lib(tmp_path / "pepunit.lib.gz", gz=True), out)
    assert _read_fasta(out)[0] == EXPECTED_IDS


def test_no_peptidase_units_is_an_error(tmp_path):
    lib = _write_lib(tmp_path / "pepunit.lib", records=[RECORDS[1], RECORDS[5]])
    cfg = PeptidaseDbBuildConfig(input=str(lib), output=str(tmp_path / "x.dmnd"), makedb=False)
    with pytest.raises(PeptidaseDbBuildError):
        build_peptidase_db(cfg)


def test_cli_fasta_only(tmp_path):
    lib = _write_lib(tmp_path / "pepunit.lib")
    res = CliRunner().invoke(cli, ["-i", str(lib), "-o", str(tmp_path / "db" / "peptidase_db.dmnd"), "--no-makedb",
                                   "--drop-inhibitors"])
    assert res.exit_code == 0, res.output
    assert _read_fasta(tmp_path / "db" / "peptidase_db.faa")[0] == EXPECTED_IDS
    assert not (tmp_path / "db" / "peptidase_db.dmnd").exists()


def test_cli_keeps_inhibitor_units_by_default(tmp_path):
    """Inhibitor units are decoys: with DIAMOND --max-target-seqs 1 they stop an inhibitor-like protein from
    taking a weaker peptidase label (run_dbcan never reports them as Peptidase)."""
    lib = _write_lib(tmp_path / "pepunit.lib")
    res = CliRunner().invoke(cli, ["-i", str(lib), "-o", str(tmp_path / "db" / "peptidase_db.dmnd"), "--no-makedb"])
    assert res.exit_code == 0, res.output
    ids = _read_fasta(tmp_path / "db" / "peptidase_db.faa")[0]
    assert set(EXPECTED_IDS) < set(ids)
    assert any(i.split("|")[1].startswith("I") for i in ids)


def test_cli_missing_diamond_exits_nonzero(tmp_path):
    lib = _write_lib(tmp_path / "pepunit.lib")
    res = CliRunner().invoke(cli, ["-i", str(lib), "-o", str(tmp_path / "peptidase_db.dmnd"),
                                   "--diamond_bin", "no-such-diamond-binary"])
    assert res.exit_code == 1


def test_cli_no_peptidase_units_exits_nonzero(tmp_path):
    lib = _write_lib(tmp_path / "pepunit.lib", records=[RECORDS[1]])
    res = CliRunner().invoke(cli, ["-i", str(lib), "-o", str(tmp_path / "peptidase_db.dmnd"), "--no-makedb"])
    assert res.exit_code == 1


@pytest.mark.skipif(shutil.which("diamond") is None, reason="DIAMOND not installed")
def test_makedb_writes_dmnd(tmp_path):
    lib = _write_lib(tmp_path / "pepunit.lib")
    out = tmp_path / "peptidase_db.dmnd"
    build_peptidase_db(PeptidaseDbBuildConfig(input=str(lib), output=str(out)))
    assert out.is_file() and out.stat().st_size > 0
