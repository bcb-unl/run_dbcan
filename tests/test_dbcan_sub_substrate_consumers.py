"""Consumers of the per-EC Substrate column of dbCANsub_hmm_results.tsv.

The column has one entry per EC token joined with ';' ('-' = EC without substrate) and
',' inside an entry that lists several substrates. overview.tsv and the CGC dbCAN-sub
scoring must read it as a list of substrates, like the web server (run_dbcan 4), whose
column was 'substrate, substrate'.
"""
from types import SimpleNamespace

import pytest

from dbcan.annotation.cgc_substrate_prediction import (
    dbCAN_substrate_prediction,
    dbSub_record,
    split_dbsub_substrates,
)
from dbcan.configs.base_config import OverviewGeneratorConfig
from dbcan.configs.cgc_substrate_config import CGCSubstrateConfig
from dbcan.IO.OverviewGenerator import OverviewGenerator
from dbcan.IO.table import read_tsv

SUB_HEADER = (
    "Subfam Name\tSubfam Composition\tSubfam EC\tSubstrate\tHMM Length\tTarget Name\tTarget Length\t"
    "i-Evalue\tHMM From\tHMM To\tTarget From\tTarget To\tCoverage\tHMM File Name\n"
)


def _sub_row(subfam, ec, substrate, gene):
    return f"{subfam}\t-\t{ec}\t{substrate}\t300\t{gene}\t400\t1e-50\t1\t290\t10\t300\t0.96\tdbCAN-sub\n"


@pytest.mark.parametrize(
    "field, expected",
    [
        ("-", []),
        ("-;-;-", []),
        ("lignin", ["lignin"]),
        ("xylan;-", ["xylan"]),
        ("alpha-glucan;alpha-glucan", ["alpha-glucan"]),
        ("glucosylglycerate;alpha-glucan", ["glucosylglycerate", "alpha-glucan"]),
        ("chitin,peptidoglycan;-", ["chitin", "peptidoglycan"]),
        ("beta-glucuronan, alginate;alginate", ["beta-glucuronan", "alginate"]),
    ],
)
def test_split_dbsub_substrates(field, expected):
    assert split_dbsub_substrates(field) == expected


def test_record_without_real_substrate_is_empty():
    rec = dbSub_record(_sub_row("GT8_e147", "2.4.1.44:3;2.4.1.58:1", "-;-", "g1").rstrip("\n").split("\t"))
    assert rec.Substrate == ""


def test_overview_substrate_drops_placeholder_entries(tmp_path):
    (tmp_path / "dbCANsub_hmm_results.tsv").write_text(
        SUB_HEADER + _sub_row("CE15_e1", "3.1.1.117:1;3.1.1.-:1", "lignin;-", "g1")
    )
    OverviewGenerator(OverviewGeneratorConfig(output_dir=str(tmp_path))).run()
    ov = read_tsv(tmp_path / "overview.tsv").set_index("Gene ID")
    assert ov.loc["g1", "Substrate"] == "lignin"


def test_cgc_scoring_splits_per_ec_entries(tmp_path):
    """Two CAZymes of one CGC: 'alpha-glucan;alpha-glucan' and 'alpha-glucan' -> alpha-glucan 2.0."""
    (tmp_path / "dbCANsub_hmm_results.tsv").write_text(
        SUB_HEADER
        + _sub_row("GH13_e62", "3.2.1.1:1;3.2.1.-:1", "alpha-glucan;alpha-glucan", "g1")
        + _sub_row("GH13_e383", "2.4.1.352:8;2.4.1.-:1", "glucosylglycerate;alpha-glucan", "g2")
        + _sub_row("GT8_e147", "2.4.1.44:3;2.4.1.58:1", "-;-", "g3")
    )
    (tmp_path / "cgc_standard_out.tsv").write_text(
        "CGC#\tGene Type\tContig ID\tProtein ID\tGene Start\tGene Stop\tGene Strand\tGene Annotation\n"
        "CGC1\tCAZyme\tc1\tg1\t1\t100\t+\tCAZyme|GH13\n"
        "CGC1\tCAZyme\tc1\tg2\t200\t300\t+\tCAZyme|GH13\n"
        "CGC1\tCAZyme\tc1\tg3\t400\t500\t+\tCAZyme|GT8\n"
    )
    sp = dbCAN_substrate_prediction(CGCSubstrateConfig(output_dir=str(tmp_path), db_dir=str(tmp_path)))
    sp.dbcan_sub_subfamily_substrate_prediction()
    cgc = next(iter(sp.dbcan_sub_CGC2scores))
    assert sp.dbcan_sub_CGC2scores[cgc] == {"alpha-glucan": 2.0, "glucosylglycerate": 1.0}
    assert sp.dbcan_sub_CGC2substrates[cgc] == "alpha-glucan"
    # g3 has no real substrate and must not count towards the cut-offs
    assert sp.cgcid2CAZyme_domain_substrate_num[cgc] == 2
