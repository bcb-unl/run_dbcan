"""One full-length label per non-CAZyme protein in cgc.gff, chosen by DIAMOND bit score.

TCDB, the prokaryotic TF DIAMOND database and SulfAtlas hold whole proteins: a hit names the
whole protein and cannot be placed at domain level. On a CGC CAZyme such labels are dropped
(test_cazyme_wins_over_full_length_labels.py). On any other protein with labels from more than
one of these databases only the best one is used: highest bit score (independent of database
size, unlike the E-value), then higher min(query, subject) coverage, then lower E-value, then
TC > TF_prok > Sulfatase. Domain-level labels (STP, fungal TF HMM, MEROPS units) are not touched.
"""
import subprocess
from types import SimpleNamespace
from unittest.mock import patch

import pandas as pd
import pytest

from dbcan.annotation.diamond import (PeptidaseDiamondProcessor, SulfatlasDiamondProcessor,
                                      TCDBDiamondProcessor, TFDiamondProcessor)
from dbcan.configs.base_config import GFFConfig
from dbcan.configs.diamond_config import (DiamondPeptidaseConfig, DiamondSulfataseConfig,
                                          DiamondTCConfig, DiamondTFConfig)
from dbcan.IO.gff import get_gff_processor
from dbcan.IO.table import read_tsv
from dbcan.process.process_utils import _best_hit_per_protein, process_cgc_sig_results

# --------------------------------------------------------------------------- DIAMOND searches

SEARCHES = [
    (TCDBDiamondProcessor, DiamondTCConfig),
    (TFDiamondProcessor, DiamondTFConfig),
    (SulfatlasDiamondProcessor, DiamondSulfataseConfig),
    (PeptidaseDiamondProcessor, DiamondPeptidaseConfig),
]


def _setup(tmp_path, config_cls):
    (tmp_path / "uniInput.faa").write_text(">p1\nMKV\n")
    cfg = config_cls(db_dir=str(tmp_path), output_dir=str(tmp_path))
    (tmp_path / cfg.db_file).write_text("")
    return cfg


@pytest.mark.parametrize("processor_cls,config_cls", SEARCHES)
def test_diamond_command_reports_bitscore_last(tmp_path, processor_cls, config_cls):
    cfg = _setup(tmp_path, config_cls)
    seen = {}

    def fake_run(cmd, **kwargs):
        seen["cmd"] = cmd
        return subprocess.CompletedProcess(cmd, 0, "", "")

    with patch("dbcan.annotation.diamond.subprocess.run", side_effect=fake_run), \
         patch.object(processor_cls, "format_results", lambda self: None):
        processor_cls(cfg).run()
    cmd = seen["cmd"]
    fields = cmd[cmd.index("--outfmt") + 1: cmd.index("--evalue")]
    assert fields[0] == "6" and fields[-1] == "bitscore"
    # the first ten fields (and thus the old column positions) are unchanged
    assert fields[1:4] == ["sseqid", "slen", "qseqid"] and len(fields) == 12


@pytest.mark.parametrize("processor_cls,config_cls", SEARCHES)
def test_formatted_output_keeps_old_columns_and_puts_bitscore_after_database(tmp_path, processor_cls, config_cls):
    cfg = _setup(tmp_path, config_cls)
    (tmp_path / cfg.output_file).write_text(
        "x|y|z|w\t300\tp1\t310\t1e-50\t1\t300\t5\t305\t97.0\t180.5\n")
    processor_cls(cfg).format_results()
    df = read_tsv(tmp_path / cfg.output_file)
    assert list(df.columns)[:10] == list(cfg.column_names)[:10]
    assert list(df.columns)[10:] == ["Database", "BITSCORE"]
    assert df["BITSCORE"].tolist() == [180.5]


@pytest.mark.parametrize("processor_cls,config_cls", SEARCHES)
def test_header_only_output_has_bitscore_last(tmp_path, processor_cls, config_cls):
    cfg = _setup(tmp_path, config_cls)
    processor_cls(cfg).format_results()
    header = (tmp_path / cfg.output_file).read_text().rstrip("\n").split("\t")
    assert header[-2:] == ["Database", "BITSCORE"]


# --------------------------------------------------------------------------- best hit within one database

def _hits(rows, with_bits=True):
    cols = ["Annotate Name", "Target Name", "i-Evalue", "Coverage"] + (["Bitscore"] if with_bits else [])
    return pd.DataFrame(rows, columns=cols)


def test_best_hit_within_database_uses_bitscore():
    # both E-values are reported as 0: only the bit score separates the hits
    df = _hits([("A", "p1", 0.0, 100.0, 900.0), ("B", "p1", 0.0, 100.0, 1200.0),
                ("C", "p2", 1e-30, 90.0, 150.0), ("D", "p2", 1e-50, 60.0, 200.0)])
    best = dict(zip(*_best_hit_per_protein(df)[["Target Name", "Annotate Name"]].T.values))
    assert best == {"p1": "B", "p2": "D"}


def test_best_hit_within_database_without_bitscore_uses_evalue_then_coverage():
    df = _hits([("A", "p1", 1e-10, 100.0), ("B", "p1", 1e-20, 50.0),
                ("C", "p2", 1e-5, 50.0), ("D", "p2", 1e-5, 80.0)], with_bits=False)
    best = dict(zip(*_best_hit_per_protein(df)[["Target Name", "Annotate Name"]].T.values))
    assert best == {"p1": "B", "p2": "D"}


# --------------------------------------------------------------------------- total_cgc_info.tsv

TC_HEADER = "TCDB ID\tTCDB Length\tTarget ID\tTarget Length\tEVALUE\tTCDB START\tTCDB END\tQSTART\tQEND\tCOVERAGE\tDatabase"
SUL_HEADER = "Sul ID\tSul Length\tTarget ID\tTarget Length\tEVALUE\tSul START\tSul END\tQSTART\tQEND\tCOVERAGE\tDatabase"
HMM_HEADER = ("HMM Name\tHMM Length\tTarget Name\tTarget Length\ti-Evalue\tHMM From\tHMM To\t"
              "Target From\tTarget To\tCoverage\tHMM File Name\n")


def _run_sig(d):
    cfg = SimpleNamespace(output_dir=str(d), prokaryotic=True, fungi=False)
    process_cgc_sig_results(cfg, cfg, None, cfg, cfg, cfg)
    return read_tsv(d / "total_cgc_info.tsv")


def test_total_cgc_info_carries_bitscore_as_last_column(tmp_path):
    (tmp_path / "diamond.out.tc").write_text(
        TC_HEADER + "\tBITSCORE\n"
        "2.A.1.1.1\t400\tp1\t410\t1e-80\t1\t400\t5\t405\t97.6\tTC\t300.4\n"
        "2.A.1.1.2\t400\tp1\t410\t1e-60\t1\t400\t5\t405\t97.6\tTC\t250.0\n")
    # a file written before the bit score existed: still read, no bit score
    (tmp_path / "diamond.out.sulfatase").write_text(
        SUL_HEADER + "\nS1_4\t500\tp1\t410\t1e-70\t1\t480\t1\t400\t97.0\tSulfatase\n")
    (tmp_path / "STP_hmm_results.tsv").write_text(
        HMM_HEADER + "LacI\t46\tp1\t410\t1e-20\t1\t46\t3\t48\t1.0\tSTP\n")
    df = _run_sig(tmp_path)
    assert list(df.columns)[11:] == ["Type", "Source", "Bitscore"]
    by = {(r["Source"], r["Annotate Name"]): r["Bitscore"] for _, r in df.iterrows()}
    assert set(by) == {("TC", "2.A.1.1.1"), ("Sulfatase", "S1_4"), ("STP", "LacI")}
    assert by[("TC", "2.A.1.1.1")] == pytest.approx(300.4)
    assert pd.isna(by[("Sulfatase", "S1_4")]) and pd.isna(by[("STP", "LacI")])
    # Annotate File Name is still the Database label (column 10)
    assert df.loc[df["Source"] == "TC", "Annotate File Name"].tolist() == ["TC"]


# --------------------------------------------------------------------------- cgc.gff labels

OVERVIEW = (
    "Gene ID\tEC#\tdbCAN_hmm\tdbCAN_sub\tDIAMOND\t#ofTools\tRecommend Results\tSubstrate\n"
    "cz\t-\tGH13_31(1-400)\tGH13_e1(1-400)\tGH13_31\t3\tGH13_31\t-\n"
)
HEAD = ("Annotate Name\tAnnotate Length\tTarget Name\tTarget Length\ti-Evalue\tAnnotate From\tAnnotate To\t"
        "Target From\tTarget To\tCoverage\tAnnotate File Name\tType\tSource\tBitscore\n")
TYPE = {"TC": "TC", "TF_prok": "TF", "TF_fungi": "TF", "Sulfatase": "Sulfatase", "STP": "STP",
        "Peptidase": "Peptidase"}


def _row(name, prot, source, bits="", evalue="1e-30", slen=400, sfrom=1, sto=400, qlen=400, qfrom=1, qto=400):
    typ = TYPE[source]
    return (f"{name}\t{slen}\t{prot}\t{qlen}\t{evalue}\t{sfrom}\t{sto}\t{qfrom}\t{qto}\t90.0\t{typ}\t{typ}\t"
            f"{source}\t{bits}\n")


def _labels(tmp_path, rows, head=HEAD):
    (tmp_path / "overview.tsv").write_text(OVERVIEW)
    (tmp_path / "total_cgc_info.tsv").write_text(head + "".join(rows))
    proc = get_gff_processor(GFFConfig(output_dir=str(tmp_path), input_gff=str(tmp_path / "x.gff"),
                                       gff_type="prodigal"))
    return {k: set(v["CGC_annotation"].split("+")) for k, v in proc.load_cgc_type().items()}


def test_non_cazyme_keeps_the_full_length_label_with_highest_bitscore(tmp_path):
    ann = _labels(tmp_path, [
        # p1: TC vs sulfatase vs TF; the E-values come from databases of different size and are
        # not used: the sulfatase hit has the worse E-value but the higher bit score
        _row("2.A.1.1.1", "p1", "TC", bits="310.0", evalue="1e-90"),
        _row("S1_4", "p1", "Sulfatase", bits="420.0", evalue="1e-80"),
        _row("LACI_ECOLI", "p1", "TF_prok", bits="120.0", evalue="1e-95"),
        # p2: TF beats TC
        _row("3.A.1.1.1", "p2", "TC", bits="90.0"),
        _row("CCPA_BACSU", "p2", "TF_prok", bits="95.5"),
        # p3: a single full-length label is never touched
        _row("1.A.1.1.1", "p3", "TC", bits="50.0"),
    ])
    assert ann["p1"] == {"Sulfatase|S1_4"}
    assert ann["p2"] == {"TF|CCPA_BACSU"}
    assert ann["p3"] == {"TC|1.A.1.1.1"}


def test_ties_go_to_coverage_then_evalue_then_tc_tf_sulfatase(tmp_path):
    ann = _labels(tmp_path, [
        # t1: equal bit scores; the TC alignment covers only 50 % of the TCDB protein
        _row("2.A.1.1.1", "t1", "TC", bits="200.0", slen=800, sfrom=1, sto=400),
        _row("S1_4", "t1", "Sulfatase", bits="200.0", slen=400),
        # t2: query coverage counts too (min of query and subject coverage)
        _row("2.A.1.1.2", "t2", "TC", bits="200.0", qlen=400, qfrom=1, qto=400),
        _row("S1_8", "t2", "Sulfatase", bits="200.0", qlen=400, qfrom=101, qto=400),
        # t3: equal bit score and coverage, lower E-value wins
        _row("2.A.1.1.3", "t3", "TC", bits="200.0", evalue="1e-40"),
        _row("CCPA_BACSU", "t3", "TF_prok", bits="200.0", evalue="1e-45"),
        # t4: everything equal -> TC > TF > Sulfatase
        _row("S1_1", "t4", "Sulfatase", bits="200.0"),
        _row("NHAR_ECOLI", "t4", "TF_prok", bits="200.0"),
        _row("2.A.1.1.4", "t4", "TC", bits="200.0"),
        # t5: everything equal, TF before Sulfatase
        _row("S1_2", "t5", "Sulfatase", bits="200.0"),
        _row("ARAC_ECOLI", "t5", "TF_prok", bits="200.0"),
    ])
    assert ann["t1"] == {"Sulfatase|S1_4"}
    assert ann["t2"] == {"TC|2.A.1.1.2"}
    assert ann["t3"] == {"TF|CCPA_BACSU"}
    assert ann["t4"] == {"TC|2.A.1.1.4"}
    assert ann["t5"] == {"TF|ARAC_ECOLI"}


def test_domain_level_labels_are_kept_next_to_the_best_full_length_label(tmp_path):
    ann = _labels(tmp_path, [
        _row("2.A.1.1.1", "p1", "TC", bits="300.0"),
        _row("S1_4", "p1", "Sulfatase", bits="100.0"),
        _row("LacI", "p1", "STP", qfrom=3, qto=48),
        _row("Peripla_BP_3", "p1", "STP", qfrom=158, qto=315),
        _row("S09.940", "p1", "Peptidase", bits="80.0"),
        _row("Zn_clus", "p1", "TF_fungi", qfrom=10, qto=50),
    ])
    assert ann["p1"] == {"TC|2.A.1.1.1", "STP|LacI", "STP|Peripla_BP_3", "Peptidase|S09.940", "TF|Zn_clus"}


def test_cazyme_rule_unchanged(tmp_path):
    ann = _labels(tmp_path, [
        _row("8.A.9.1.1", "cz", "TC", bits="500.0"),
        _row("S1_7", "cz", "Sulfatase", bits="600.0"),
        _row("LacI", "cz", "STP"),
        _row("M38.UNW", "cz", "Peptidase", bits="70.0"),
    ])
    assert ann["cz"] == {"CAZyme|GH13_31", "STP|LacI", "Peptidase|M38.UNW"}


def test_old_table_without_bitscore_keeps_all_labels(tmp_path):
    old_head = HEAD.replace("\tBitscore", "")
    rows = [_row("2.A.1.1.1", "p1", "TC").rsplit("\t", 1)[0] + "\n",
            _row("S1_4", "p1", "Sulfatase").rsplit("\t", 1)[0] + "\n"]
    ann = _labels(tmp_path, rows, head=old_head)
    assert ann["p1"] == {"TC|2.A.1.1.1", "Sulfatase|S1_4"}


def test_protein_with_a_label_without_bitscore_keeps_all_labels(tmp_path):
    ann = _labels(tmp_path, [
        _row("2.A.1.1.1", "p1", "TC", bits="300.0"),
        _row("S1_4", "p1", "Sulfatase", bits=""),      # e.g. an old diamond.out.sulfatase
        _row("2.A.1.1.2", "p2", "TC", bits="300.0"),
        _row("S1_5", "p2", "Sulfatase", bits="100.0"),
    ])
    assert ann["p1"] == {"TC|2.A.1.1.1", "Sulfatase|S1_4"}
    assert ann["p2"] == {"TC|2.A.1.1.2"}
