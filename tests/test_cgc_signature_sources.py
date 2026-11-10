"""CGC signature merge (total_cgc_info.tsv): no overlap filtering across sources.

TCDB, PRODORIC (TF), SulfAtlas and MEROPS are full-length protein databases searched
with DIAMOND: a hit is a protein-level label and its query interval is an alignment
span, not a domain. The old code put all sources in one table and removed
"overlapping" rows by comparing DIAMOND spans with Pfam domain coordinates and DIAMOND
E-values with HMMER i-Evalues. Rows below are taken from a real run on a
Bacteroidetes PUL (PUL0230_GCMeta_00102082_C2), where that removed four STP domains.
"""
from types import SimpleNamespace

import pandas as pd

from dbcan.IO.table import read_tsv
from dbcan.process.process_utils import process_cgc_sig_results

TC_HEADER = "TCDB ID\tTCDB Length\tTarget ID\tTarget Length\tEVALUE\tTCDB START\tTCDB END\tQSTART\tQEND\tCOVERAGE\tDatabase\n"
TF_HEADER = "TF ID\tTF Length\tTarget ID\tTarget Length\tEVALUE\tTF START\tTF END\tQSTART\tQEND\tCOVERAGE\tDatabase\n"
HMM_HEADER = ("HMM Name\tHMM Length\tTarget Name\tTarget Length\ti-Evalue\tHMM From\tHMM To\t"
              "Target From\tTarget To\tCoverage\tHMM File Name\n")


def _write_inputs(d, extra_tc=""):
    (d / "diamond.out.tc").write_text(
        TC_HEADER
        + "3.A.1.1.27\t415\tp12\t408\t2.62e-85\t1\t411\t1\t405\t99.3\tTC\n"
        + "3.A.1.1.20\t368\tp13\t368\t5.77e-234\t1\t368\t1\t368\t100.0\tTC\n"
        + extra_tc
    )
    (d / "diamond.out.tf").write_text(
        TF_HEADER + "CCPA_BACSU\t334\tp20\t316\t5.27e-34\t5\t333\t3\t316\t99.4\tprodoric\n"
    )
    (d / "STP_hmm_results.tsv").write_text(
        HMM_HEADER
        + "SBP_bac_1\t315\tp12\t408\t4.85e-28\t9\t314\t54\t321\t0.97\tSTP\n"
        + "TOBE\t64\tp13\t368\t2.39e-07\t6\t64\t304\t357\t0.92\tSTP\n"
        + "LacI\t46\tp20\t316\t4.62e-22\t1\t46\t3\t48\t1.0\tSTP\n"
        + "Peripla_BP_3\t164\tp20\t316\t2.24e-33\t2\t164\t158\t315\t0.99\tSTP\n"
    )


def _run(d):
    cfg = SimpleNamespace(output_dir=str(d), prokaryotic=True, fungi=False)
    process_cgc_sig_results(cfg, cfg, None, cfg, cfg, cfg)
    return read_tsv(d / "total_cgc_info.tsv")


def test_stp_domains_are_not_removed_by_full_length_diamond_hits(tmp_path):
    _write_inputs(tmp_path)
    df = _run(tmp_path)
    got = set(zip(df["Target Name"], df["Annotate Name"], df["Type"]))
    assert got == {
        ("p12", "3.A.1.1.27", "TC"), ("p12", "SBP_bac_1", "STP"),
        ("p13", "3.A.1.1.20", "TC"), ("p13", "TOBE", "STP"),
        ("p20", "CCPA_BACSU", "TF"), ("p20", "LacI", "STP"), ("p20", "Peripla_BP_3", "STP"),
    }


def test_diamond_source_keeps_one_best_hit_per_protein(tmp_path):
    # a second, weaker TC hit on p12 (e.g. another HSP) is a duplicate protein-level label
    _write_inputs(tmp_path, extra_tc="3.A.1.1.99\t500\tp12\t408\t1e-20\t1\t300\t10\t300\t70.0\tTC\n")
    df = _run(tmp_path)
    tc = df[(df["Type"] == "TC") & (df["Target Name"] == "p12")]
    assert tc["Annotate Name"].tolist() == ["3.A.1.1.27"]


def test_gff_annotation_lists_all_sources(tmp_path):
    """The merged labels reach cgc.gff; the gene type of a TF+STP protein is TF."""
    from dbcan.annotation.CGCFinder import CGCFinder
    from dbcan.configs.base_config import GFFConfig
    from dbcan.IO.gff import get_gff_processor

    _write_inputs(tmp_path)
    _run(tmp_path)
    (tmp_path / "overview.tsv").write_text(
        "Gene ID\tEC#\tdbCAN_hmm\tdbCAN_sub\tDIAMOND\t#ofTools\tRecommend Results\n"
    )
    (tmp_path / "in.gff").write_text(
        "c1\tx\tCDS\t100\t1300\t.\t-\t0\tID=p12\n"
        "c1\tx\tCDS\t1400\t2500\t.\t-\t0\tID=p13\n"
        "c1\tx\tCDS\t2600\t3550\t.\t-\t0\tID=p20\n"
    )
    cfg = GFFConfig(output_dir=str(tmp_path), input_gff=str(tmp_path / "in.gff"), gff_type="prodigal")
    assert get_gff_processor(cfg).process_gff()
    gff = pd.read_csv(tmp_path / "cgc.gff", sep="\t", header=None, dtype=str)[8].tolist()
    ann = {a.split(";")[0].split("=")[1]: a.split("CGC_annotation=")[1] for a in gff}
    assert set(ann["p12"].split("+")) == {"TC|3.A.1.1.27", "STP|SBP_bac_1"}
    assert set(ann["p13"].split("+")) == {"TC|3.A.1.1.20", "STP|TOBE"}
    assert set(ann["p20"].split("+")) == {"TF|CCPA_BACSU", "STP|LacI", "STP|Peripla_BP_3"}
    assert CGCFinder.get_gene_type(ann["p20"]) == "TF"
    assert CGCFinder.get_gene_type("CAZyme|GH13_31+TC|8.A.9.1.1") == "CAZyme"


def test_tf_route_bacteria_diamond_fungi_hmm(tmp_path):
    """Default: TF DIAMOND only. --fungi: TF HMM only. Explicit --prokaryotic --fungi: both."""
    from unittest.mock import patch

    from click.testing import CliRunner

    from dbcan.main import cli

    gff = tmp_path / "in.gff"
    gff.write_text("c1\tx\tCDS\t1\t90\t.\t+\t0\tID=p1\n")
    base = ["gff_process", "--output_dir", str(tmp_path), "--db_dir", str(tmp_path),
            "--input_gff", str(gff), "--gff_type", "prodigal"]
    expected = {(): (True, False), ("--fungi",): (False, True),
                ("--fungi", "--prokaryotic"): (True, True), ("--no-prokaryotic",): (False, False)}
    for extra, (tf_diamond, tf_hmm) in expected.items():
        seen = {}

        def fake(tc, tfd, tf, stp, sul, pep, gffcfg):
            seen.update(diamond=tfd.prokaryotic, hmm=tf.fungi)

        with patch("dbcan.core.run_dbCAN_CGCFinder_preprocess", side_effect=fake):
            result = CliRunner().invoke(cli, base + list(extra))
        assert result.exit_code == 0, result.output
        assert (seen["diamond"], seen["hmm"]) == (tf_diamond, tf_hmm), extra


def test_peptidase_signature_excludes_inhibitor_units(tmp_path):
    """MEROPS pepunit.lib also holds inhibitor units (I families) and a few malformed IDs."""
    _write_inputs(tmp_path)
    (tmp_path / "diamond.out.peptidase").write_text(
        "Peptidase ID\tPeptidase Length\tTarget ID\tTarget Length\tEVALUE\tPeptidase START\tPeptidase END\t"
        "QSTART\tQEND\tCOVERAGE\tDatabase\n"
        "S09.940\t270\tpA\t277\t4.7e-205\t1\t270\t4\t274\t100\tPeptidase\n"
        "I39.UPW\t161\tpB\t1037\t1.4e-60\t1\t161\t73\t228\t98.8\tPeptidase\n"
        "I04.UPW\t160\tpC\t674\t2.6e-21\t1\t160\t388\t544\t94.0\tPeptidase\n"
        "Bacteroides\t300\tpD\t500\t1e-30\t1\t300\t10\t300\t95.0\tPeptidase\n"
        "C56.UNW\t180\tpE\t200\t1e-40\t1\t180\t5\t185\t99.0\tPeptidase\n"
    )
    df = _run(tmp_path)
    pep = df[df["Type"] == "Peptidase"]
    assert sorted(zip(pep["Target Name"], pep["Annotate Name"])) == [("pA", "S09.940"), ("pE", "C56.UNW")]
