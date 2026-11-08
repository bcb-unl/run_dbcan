"""run_dbcan expression without DESeq2: counting and abundance go through dbcan_utils cal_coverage / *_abund."""
import random

import pandas as pd
import pysam

from dbcan.configs.expression_config import ExpressionConfig
from dbcan.expression.pipeline import run_expression_pipeline

L = 3000
random.seed(7)
REF = "".join(random.choice("ACGT") for _ in range(L))
HEADER = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "c1", "LN": L}]})
GENES = [("gA", 101, 1000, "+"), ("gB", 1101, 1600, "+"), ("gD", 2001, 2500, "-")]


def _bam(path, fragments):
    """``fragments``: list of (name, pos1, pos2) read pairs (100M each)."""
    recs = []
    for name, p1, p2 in fragments:
        for pos, mate, flag in ((p1, p2, 0x1 | 0x2 | 0x20 | 0x40), (p2, p1, 0x1 | 0x2 | 0x10 | 0x80)):
            a = pysam.AlignedSegment(HEADER)
            a.query_name, a.reference_id, a.reference_start = name, 0, pos - 1
            a.cigarstring, a.mapping_quality, a.flag = "100M", 60, flag
            a.next_reference_id, a.next_reference_start = 0, mate - 1
            a.query_sequence = REF[pos - 1:pos - 1 + 100]
            a.query_qualities = pysam.qualitystring_to_array("I" * 100)
            a.set_tags([("NM", 0), ("MD", "100")])
            recs.append(a)
    recs.sort(key=lambda r: r.reference_start)
    with pysam.AlignmentFile(str(path), "wb", header=HEADER) as out:
        for r in recs:
            out.write(r)
    pysam.index(str(path))
    return str(path)


def _dbcan_dir(tmp_path):
    d = tmp_path / "sample.dbCAN"
    d.mkdir()
    (d / "overview.tsv").write_text(
        "Gene ID\tEC#\tdbCAN_hmm\tdbCAN_sub\tDIAMOND\t#ofTools\tRecommend Results\tSubstrate\n"
        "gA\t-\tGH13(1-290)\tGH13_e1(1-290)\tGH13\t3\tGH13\t-\n"
        "gD\t-\tGT2(1-160)\tGT2_e1(1-160)\tGT2\t3\tGT2\t-\n")
    (d / "cgc_standard_out.tsv").write_text(
        "CGC#\tGene Type\tContig ID\tProtein ID\tGene Start\tGene Stop\tGene Strand\tGene Annotation\n"
        "CGC1\tCAZyme\tc1\tgA\t101\t1000\t+\tCAZyme|GH13\n"
        "CGC1\tTC\tc1\tgB\t1101\t1600\t+\tTC|1.A.1\n"
        "CGC1\tCAZyme\tc1\tgD\t2001\t2500\t-\tCAZyme|GT2\n")
    with open(d / "cgc.gff", "w") as fh:
        for gid, s, e, st in GENES:
            fh.write(f"c1\t.\tgene\t{s}\t{e}\t.\t{st}\t.\tprotein_id={gid};CGC_annotation=null\n")
    return d


def test_counts_are_fragments_and_tables_are_written(tmp_path):
    d = _dbcan_dir(tmp_path)
    b1 = _bam(tmp_path / "s1.bam", [("f1", 201, 401), ("f2", 301, 501), ("f3", 1201, 1301)])
    b2 = _bam(tmp_path / "s2.bam", [("f1", 2101, 2301)])
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(f"sample_id\tcondition\tbam\nS1\tcontrol\t{b1}\nS2\ttreat\t{b2}\n")
    out = tmp_path / "expr"
    cfg = ExpressionConfig(input_dir=str(d), output_dir=str(out), samplesheet=str(sheet), skip_alignment=True,
                           skip_deseq2=True, abundance="TPM", run_plots=True, also_heatmap=True)
    run_expression_pipeline(cfg)

    counts = pd.read_csv(out / "count_matrix.tsv", sep="\t", index_col=0)
    # both mates of a pair count once (cal_coverage fragments mode)
    assert counts.loc["S1", "gA"] == 2 and counts.loc["S1", "gB"] == 1 and counts.loc["S1", "gD"] == 0
    assert counts.loc["S2", "gD"] == 1
    long = pd.read_csv(out / "cgc_gene_counts.tsv", sep="\t")
    assert set(long["cgcid"]) == {"c1|CGC1"} and len(long) == 6
    # per-sample dbcan_utils abundance tables
    assert (out / "S1" / "fam_abund.out").is_file()
    assert (out / "S1" / "CGC_abund.out").is_file()
    assert (out / "cgc_expression_matrix.tsv").is_file()
    # CGC expression figures (gene track helpers in dbcan.plot.expression_track)
    figs = out / "plots" / "CGC_expression"
    assert (figs / "c1_CGC1.expression.pdf").is_file()
    assert (out / "plots" / "CGC_expression_heatmap.pdf").is_file()


def test_dbcan_plot_expression_commands_fail_without_inputs(tmp_path):
    from click.testing import CliRunner
    from dbcan.utils.plots import cli

    d = _dbcan_dir(tmp_path)
    expr = tmp_path / "expr"
    expr.mkdir()
    for fn in ("DEG_volcano", "DEG_CGC_bar", "CGC_expression_heatmap", "CGC_expression_plot"):
        res = CliRunner().invoke(cli, [fn, "-i", str(d), "--expression-dir", str(expr)])
        assert res.exit_code == 1, (fn, res.output)
    res = CliRunner().invoke(cli, ["DEG_volcano", "-i", str(d)])
    assert res.exit_code == 1 and "--expression-dir is required" in res.output
