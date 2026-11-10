"""dbcan_utils: gene read counting (cal_coverage) and family / substrate / CGC abundance."""
import os
import random

import pysam
import pytest
from click.testing import CliRunner

from dbcan.configs.utils_config import AbundanceConfig, CoverageConfig
from dbcan.utils.utils import (
    AbundanceError,
    OverView,
    Read_cgc_standard_out,
    Read_dbcansub_out,
    alignment_identity,
    cal_coverage,
    cli,
    parse_ec_field,
    read_gene_regions,
    run_abundance,
    split_substrates,
)

L = 3000
random.seed(1)
REF = "".join(random.choice("ACGT") for _ in range(L))

# 1-based inclusive
GENES = [
    ("gA", 101, 1000, "+"),
    ("gB", 1001, 1600, "+"),
    ("gC", 1500, 2200, "-"),   # overlaps gB (1500-1600) on the other strand
    ("gD", 2401, 2500, "+"),
]


# ---------------------------------------------------------------------------
# parsing
# ---------------------------------------------------------------------------
def test_parse_ec_field_splits_domains_and_tokens():
    assert parse_ec_field("3.2.1.1:3;3.2.1.20:13|3.2.1.1;3.2.1.20;3.2.1.133") == ["3.2.1.1", "3.2.1.20", "3.2.1.133"]
    assert parse_ec_field("3.2.1.37:20;3.2.1.55:3;3.2.1.-|3.2.1.37:35;3.2.1.55:6;3.2.1.72") == [
        "3.2.1.37", "3.2.1.55", "3.2.1.-", "3.2.1.72"]
    assert parse_ec_field("-|-|-") == []
    assert parse_ec_field("-") == []
    assert parse_ec_field("2.4.1.352:8;2.4.1.-") == ["2.4.1.352", "2.4.1.-"]


def test_overview_preds_have_no_joined_ec_keys():
    ov = OverView("p1", "3.2.1.1:3;3.2.1.20:13|3.2.1.1;3.2.1.20", "CBM34(5-111)+GH13_21(188-486)",
                  "CBM34_e8(5-112)+GH13_e88(176-491)", "GH13_21+CBM34")
    assert ov.ECs == ["3.2.1.1", "3.2.1.20"]
    assert not any(";" in p or ":" in p or "|" in p for p in ov.preds)
    assert set(ov.preds) == {"CBM34", "GH13_21", "CBM34_e8", "GH13_e88", "3.2.1.1", "3.2.1.20"}


def test_split_substrates():
    assert split_substrates("xylan;xylan;arabinan") == ["xylan", "arabinan"]
    assert split_substrates("-;-;-") == []
    assert split_substrates("chitin,peptidoglycan;-") == ["chitin", "peptidoglycan"]
    assert split_substrates("beta-glucan, xylan") == ["beta-glucan", "xylan"]
    # a comma between digits is part of a chemical name
    assert split_substrates("alpha-1,3-glucan;-") == ["alpha-1,3-glucan"]
    assert split_substrates("") == [] and split_substrates("-") == []


def test_read_dbcansub_keeps_every_row_of_a_protein(tmp_path):
    f = tmp_path / "dbCANsub_hmm_results.tsv"
    f.write_text(
        "Subfam Name\tSubfam Composition\tSubfam EC\tSubstrate\tHMM Length\tTarget Name\tTarget Length\ti-Evalue\n"
        "GH43_e285\tx\t3.2.1.37:20;3.2.1.55:3;3.2.1.-\txylan;xylan;arabinan\t292\tP1\t500\t1e-50\n"
        "CBM91_e64\tx\t3.2.1.37:35\t-;-;-\t205\tP1\t500\t1e-30\n"
        "GH13_e88\tx\t3.2.1.1\talpha-glucan,starch\t317\tP2\t500\t1e-30\n"
    )
    recs = Read_dbcansub_out(str(f))
    assert len(recs["P1"].hits) == 2
    assert recs["P1"].substrates == ["xylan", "arabinan"]
    assert recs["P2"].substrates == ["alpha-glucan", "starch"]


def test_cgc_standard_placeholder_rows_are_skipped(tmp_path):
    f = tmp_path / "cgc_standard_out.tsv"
    f.write_text(
        "CGC#\tGene Type\tContig ID\tProtein ID\tGene Start\tGene Stop\tGene Strand\tGene Annotation\n"
        "CGC1\tCAZyme\tc1\tgA\t101\t1000\t+\tCAZyme|GH13\n"
        "CGC1\tOther\tc1\tNA\t1001\t1080\t+\tOther|tRNA\n"
        "CGC1\tTC\tc1\tgB\t1001\t1600\t+\tTC|1.A.1\n"
        "CGC2\tOther\tc1\tNA\t3001\t3080\t+\tOther|tRNA\n"
        "CGC2\tCAZyme\tc1\tgD\t2401\t2500\t+\tCAZyme|GT2\n"
    )
    seq2rec, cgc2recs = Read_cgc_standard_out(str(f))
    assert "NA" not in seq2rec
    assert [r.seqid for r in cgc2recs["c1|CGC1"]] == ["gA", "gB"]
    assert [r.seqid for r in cgc2recs["c1|CGC2"]] == ["gD"]


# ---------------------------------------------------------------------------
# GFF gene IDs
# ---------------------------------------------------------------------------
def _write_gff(path, style):
    with open(path, "w") as fh:
        fh.write("##gff-version 3\n")
        if style == "ncbi":
            fh.write(f"c1\tRefSeq\tregion\t1\t{L}\t.\t+\t.\tID=c1:1..{L}\n")
        for gid, s, e, st in GENES:
            if style == "prodigal":
                fh.write(f"c1\tsim\tCDS\t{s}\t{e}\t.\t{st}\t0\tID={gid};partial=00\n")
            elif style == "cgcgff":
                fh.write(f"c1\t.\tgene\t{s}\t{e}\t.\t{st}\t.\tprotein_id={gid};CGC_annotation=null\n")
            elif style == "ncbi":
                fh.write(f"c1\tRefSeq\tgene\t{s}\t{e}\t.\t{st}\t.\tID=gene-{gid};Name=b{gid};locus_tag=b{gid}\n")
                fh.write(f"c1\tRefSeq\tCDS\t{s}\t{e}\t.\t{st}\t0\tID=cds-{gid};Parent=gene-{gid};"
                         f"Name={gid};protein_id={gid};locus_tag=b{gid}\n")
        if style == "cgcgff":
            # non-coding features get the placeholder protein ID
            fh.write("c1\t.\tgene\t2600\t2680\t.\t+\t.\tprotein_id=NA;CGC_annotation=Other|tRNA\n")
            fh.write("c1\t.\tgene\t2700\t2780\t.\t+\t.\tprotein_id=NA;CGC_annotation=Other|tRNA\n")
        if style == "ncbi":
            fh.write("c1\tRefSeq\tgene\t2600\t2680\t.\t+\t.\tID=gene-t1;Name=t1\n")
            fh.write("c1\tRefSeq\ttRNA\t2600\t2680\t.\t+\t.\tID=rna-t1;Parent=gene-t1\n")


@pytest.mark.parametrize("style", ["prodigal", "cgcgff", "ncbi"])
def test_gene_ids_from_any_gff_style(tmp_path, style):
    gff = tmp_path / f"{style}.gff"
    _write_gff(gff, style)
    genes = read_gene_regions(str(gff))
    assert sorted(g.seqid for g in genes) == ["gA", "gB", "gC", "gD"]
    gA = next(g for g in genes if g.seqid == "gA")
    assert gA.length == 900 and gA.intervals == [("c1", 100, 1000)]


def test_multi_exon_cds_is_one_gene(tmp_path):
    gff = tmp_path / "euk.gff"
    gff.write_text(
        "c1\tsrc\tgene\t100\t900\t.\t+\t.\tID=g1\n"
        "c1\tsrc\tmRNA\t100\t900\t.\t+\t.\tID=m1;Parent=g1\n"
        "c1\tsrc\tCDS\t100\t300\t.\t+\t0\tID=cds1;Parent=m1;protein_id=XP_1.1\n"
        "c1\tsrc\tCDS\t501\t900\t.\t+\t0\tID=cds1;Parent=m1;protein_id=XP_1.1\n"
    )
    (g,) = read_gene_regions(str(gff))
    assert g.seqid == "XP_1.1" and g.length == 201 + 400
    # funannotate/BRAKER style: proteins are named after the mRNA (CDS Parent)
    (g2,) = read_gene_regions(str(gff), id_attr="Parent", feature_type="CDS")
    assert g2.seqid == "m1" and g2.length == 601


# ---------------------------------------------------------------------------
# BAM fixtures
# ---------------------------------------------------------------------------
HEADER = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "c1", "LN": L}]})


class BamBuilder:
    def __init__(self):
        self.recs = []

    def add(self, name, pos1, cigar="100M", flag=0, mapq=60, md=True, nm=True, mism=0, mate=None):
        a = pysam.AlignedSegment(HEADER)
        a.query_name = name
        a.reference_id = 0
        a.reference_start = pos1 - 1
        a.cigarstring = cigar
        a.mapping_quality = mapq
        a.flag = flag
        q, rpos, mlen = [], pos1 - 1, 0
        mism_pos = [10 + i * 7 for i in range(mism)]
        for op, n in a.cigartuples:
            if op == 4:
                q.append("N" * n)
            elif op in (0, 7, 8):
                seg = list(REF[rpos:rpos + n])
                for p in mism_pos:
                    seg[p] = "A" if seg[p] != "A" else "C"
                q.append("".join(seg))
                rpos += n
                mlen += n
        a.query_sequence = "".join(q)
        a.query_qualities = pysam.qualitystring_to_array("I" * len(a.query_sequence))
        if mate is not None:
            a.next_reference_id = 0
            a.next_reference_start = mate - 1
        tags = []
        if nm:
            tags.append(("NM", mism))
        if md:
            if mism:
                parts, last = [], 0
                for p in mism_pos:
                    parts += [str(p - last), REF[pos1 - 1 + p]]
                    last = p + 1
                parts.append(str(mlen - last))
                tags.append(("MD", "".join(parts)))
            else:
                tags.append(("MD", str(mlen)))
        a.set_tags(tags)
        self.recs.append(a)
        return a

    def pair(self, name, pos1, pos2, **kw):
        self.add(name, pos1, flag=0x1 | 0x2 | 0x20 | 0x40, mate=pos2, **kw)
        self.add(name, pos2, flag=0x1 | 0x2 | 0x10 | 0x80, mate=pos1, **kw)

    def write(self, path):
        self.recs.sort(key=lambda r: r.reference_start)
        with pysam.AlignmentFile(str(path), "wb", header=HEADER) as out:
            for r in self.recs:
                out.write(r)
        pysam.index(str(path))
        return str(path)


def _counts(path):
    out = {}
    for line in open(path):
        gid, length, count = line.rstrip("\n").split("\t")
        out[gid] = (int(length), int(count))
    return out


def _run(tmp_path, bam, gff_style="prodigal", **kw):
    gff = tmp_path / f"{gff_style}.gff"
    _write_gff(gff, gff_style)
    out = tmp_path / "sub" / "depth.txt"
    cal_coverage(CoverageConfig(input=bam, gff=str(gff), output=str(out), **kw))
    return _counts(out)


@pytest.fixture
def edge_bam(tmp_path):
    b = BamBuilder()
    b.pair("pairA", 201, 401)                         # 1 fragment in gA
    b.add("noMD_noNM", 501, md=False, nm=False)       # no identity tags: counted, filter skipped
    b.add("NMonly", 521, md=False)                    # NM only (HISAT2 / STAR style)
    b.add("secondary", 601, flag=0x100)
    b.add("qcfail", 651, flag=0x200)
    b.add("suppl", 701, flag=0x800)
    b.add("dup", 751, flag=0x400)
    b.add("lowmapq", 771, mapq=0)
    b.add("mm2", 301, mism=2)                         # identity 0.98 -> counted
    b.add("mm3", 311, mism=3)                         # 0.97 -> rejected
    b.add("antisense_overlap", 1501)                  # gB and gC overlap: counted for both (unstranded)
    b.add("gA_first_base", 2, cigar="100M")           # 1-based 2-101: 1 base on gA -> below ratio
    b.add("gD_start20", 2321)                         # 20 bases in gD (ratio 0.20) -> counted
    b.add("gD_end20", 2481)                           # 20 bases in gD -> counted
    return b.write(tmp_path / "edge.bam")


def test_edge_cases_exact_counts(tmp_path, edge_bam):
    c = _run(tmp_path, edge_bam)
    # gA: pairA (1 fragment) + noMD_noNM + NMonly + mm2
    assert c["gA"] == (900, 4)
    assert c["gB"] == (600, 1)
    assert c["gC"] == (701, 1)
    assert c["gD"] == (100, 2)


def test_count_mode_reads_counts_both_mates(tmp_path, edge_bam):
    c = _run(tmp_path, edge_bam, count_mode="reads")
    assert c["gA"][1] == 5


def test_duplicates_can_be_included(tmp_path, edge_bam):
    assert _run(tmp_path, edge_bam, include_duplicates=True)["gA"][1] == 5


def test_mapq_threshold_option(tmp_path, edge_bam):
    assert _run(tmp_path, edge_bam, mapping_quality=0)["gA"][1] == 5  # lowmapq kept


def test_identity_zero_disables_filter(tmp_path, edge_bam):
    assert _run(tmp_path, edge_bam, identity=0)["gA"][1] == 5  # mm3 kept


@pytest.mark.parametrize("style", ["cgcgff", "ncbi"])
def test_same_counts_with_cgc_gff_and_ncbi_gff(tmp_path, edge_bam, style):
    assert _run(tmp_path, edge_bam, gff_style=style) == _run(tmp_path, edge_bam, gff_style="prodigal")


def test_gene_start_base_is_counted(tmp_path):
    b = BamBuilder()
    # 1-based 51-150 covers gA 101-150: 50 bases (ratio .5); the old code started one base late
    b.add("start", 51)
    # read ending exactly on gA's first base (1-based 2-101): 1 base, counted with ratio 0.01
    b.add("first_base_only", 2)
    bam = b.write(tmp_path / "s.bam")
    assert _run(tmp_path, bam, overlap_base_ratio=0.01)["gA"][1] == 2


def test_stranded_counting(tmp_path):
    b = BamBuilder()
    # dUTP library: read 1 antisense, read 2 sense to the transcript.
    # gB is '+': read1 reverse, read2 forward ; gC is '-': read1 forward, read2 reverse
    b.add("gB_frag", 1100, flag=0x1 | 0x2 | 0x10 | 0x40, mate=1200)
    b.add("gB_frag", 1200, flag=0x1 | 0x2 | 0x20 | 0x80, mate=1100)
    b.add("gC_frag", 1900, flag=0x1 | 0x2 | 0x20 | 0x40, mate=2000)
    b.add("gC_frag", 2000, flag=0x1 | 0x2 | 0x10 | 0x80, mate=1900)
    # fragment inside the gB/gC overlap from a gC transcript
    b.add("ov_gC", 1501, flag=0x1 | 0x2 | 0x20 | 0x40, mate=1501)
    b.add("ov_gC", 1501, flag=0x1 | 0x2 | 0x10 | 0x80, mate=1501)
    bam = b.write(tmp_path / "st.bam")
    un = _run(tmp_path, bam)
    assert un["gB"][1] == 2 and un["gC"][1] == 2
    rev = _run(tmp_path, bam, strandedness="reverse")
    assert rev["gB"][1] == 1 and rev["gC"][1] == 2
    fwd = _run(tmp_path, bam, strandedness="forward")
    assert fwd["gB"][1] == 1 and fwd["gC"][1] == 0  # each read 1 checked: all are antisense


def test_multiprocess_matches_single(tmp_path, edge_bam):
    assert _run(tmp_path, edge_bam, threads=2) == _run(tmp_path, edge_bam, threads=1)


def test_alignment_identity_sources():
    b = BamBuilder()
    r = b.add("x", 301, mism=2, md=False)
    assert alignment_identity(r) == pytest.approx(0.98)
    r = b.add("y", 301, mism=2, nm=False)
    assert alignment_identity(r) == pytest.approx(0.98)
    r = b.add("z", 301, md=False, nm=False)
    assert alignment_identity(r) is None
    r = b.add("eqx", 301, cigar="50=1X49=", md=False, nm=False)
    assert alignment_identity(r) == pytest.approx(0.99)


def test_cal_coverage_contig_mismatch_fails(tmp_path, edge_bam):
    gff = tmp_path / "other.gff"
    gff.write_text("contigX\tsim\tCDS\t1\t100\t.\t+\t0\tID=g1\n")
    with pytest.raises(AbundanceError, match="None of the GFF contigs"):
        cal_coverage(CoverageConfig(input=edge_bam, gff=str(gff), output=str(tmp_path / "o.txt")))


# ---------------------------------------------------------------------------
# abundance
# ---------------------------------------------------------------------------
def _dbcan_dir(tmp_path):
    d = tmp_path / "dbcan"
    d.mkdir()
    (d / "overview.tsv").write_text(
        "Gene ID\tEC#\tdbCAN_hmm\tdbCAN_sub\tDIAMOND\t#ofTools\tRecommend Results\tSubstrate\n"
        "gA\t3.2.1.1:3;3.2.1.20:13|3.2.1.1;3.2.1.133\tCBM34(5-111)+GH13_21(188-486)\tCBM34_e8(5-112)+GH13_e88(176-491)\tGH13_21\t3\tx\t-\n"
        "gD\t-\tGT2(1-90)\tGT2_e1(1-90)\tGT2\t3\tGT2\t-\n"
    )
    (d / "dbCANsub_hmm_results.tsv").write_text(
        "Subfam Name\tSubfam Composition\tSubfam EC\tSubstrate\tHMM Length\tTarget Name\tTarget Length\ti-Evalue\n"
        "CBM34_e8\tx\t3.2.1.1:3\t-;-\t107\tgA\t300\t1e-30\n"
        "GH13_e88\tx\t3.2.1.1;3.2.1.20\talpha-glucan;starch\t317\tgA\t300\t1e-30\n"
        "GT2_e1\tx\t-\t-\t90\tgD\t33\t1e-30\n"
    )
    (d / "cgc_standard_out.tsv").write_text(
        "CGC#\tGene Type\tContig ID\tProtein ID\tGene Start\tGene Stop\tGene Strand\tGene Annotation\n"
        "CGC1\tCAZyme\tc1\tgA\t101\t1000\t+\tCAZyme|GH13_21\n"
        "CGC1\tOther\tc1\tNA\t1001\t1080\t+\tOther|tRNA\n"
        "CGC1\tTC\tc1\tgB\t1001\t1600\t+\tTC|1.A.1\n"
        "CGC2\tOther\tc1\tNA\t2600\t2680\t+\tOther|tRNA\n"
        "CGC2\tCAZyme\tc1\tgD\t2401\t2500\t+\tCAZyme|GT2\n"
    )
    (d / "substrate_prediction.tsv").write_text(
        "#cgcid\tPULID\tdbCAN-PUL substrate\tbitscore\tsignature pairs\tdbCAN-sub substrate\tdbCAN-sub substrate score\n"
        "c1|CGC1\tPUL1\tstarch\t100\tCAZyme-CAZyme\talpha-glucan,starch\t2.0\n"
        "c1|CGC2\t\t\t\t\t\t\n"
    )
    return d


def _depth(tmp_path, rows):
    p = tmp_path / "depth.txt"
    p.write_text("".join(f"{g}\t{length}\t{c}\n" for g, length, c in rows))
    return str(p)


def _table(path):
    lines = open(path).read().splitlines()
    return {l.split("\t")[0]: l.split("\t")[1:] for l in lines[1:]}


def test_fam_abund_ec_and_output_dir(tmp_path):
    d = _dbcan_dir(tmp_path)
    depth = _depth(tmp_path, [("gA", 1000, 300), ("gB", 500, 100), ("gC", 1000, 0), ("gD", 100, 100)])
    out = tmp_path / "abund_out"
    run_abundance(AbundanceConfig(input_dir=str(d), bedtools_depth=depth, output_dir=str(out), abundance="RPM"), "fam_abund")
    ec = _table(out / "EC_abund.out")
    assert set(ec) == {"3.2.1.1", "3.2.1.20", "3.2.1.133"}
    assert float(ec["3.2.1.1"][0]) == pytest.approx(300 / 500 * 1e6)
    fam = _table(out / "fam_abund.out")
    assert float(fam["GT2"][0]) == pytest.approx(100 / 500 * 1e6)
    assert "GT2_e1" in _table(out / "subfam_abund.out")
    assert not (tmp_path / "fam_abund.out").exists()


def test_tpm_and_fpkm(tmp_path):
    d = _dbcan_dir(tmp_path)
    depth = _depth(tmp_path, [("gA", 1000, 300), ("gB", 500, 100), ("gC", 1000, 0), ("gD", 100, 100)])
    denom = 300 / 1000 + 100 / 500 + 100 / 100
    for method, expect in (("TPM", (100 / 100) / denom * 1e6), ("FPKM", 100 / (500 / 1e6) / 0.1),
                           ("RPKM", 100 / (500 / 1e6) / 0.1)):
        out = tmp_path / method
        run_abundance(AbundanceConfig(str(d), depth, str(out), method), "fam_abund")
        assert float(_table(out / "fam_abund.out")["GT2"][0]) == pytest.approx(expect, rel=1e-3)


def test_fam_substrate_abund_uses_all_rows(tmp_path):
    d = _dbcan_dir(tmp_path)
    depth = _depth(tmp_path, [("gA", 1000, 300), ("gB", 500, 100), ("gD", 100, 100)])
    run_abundance(AbundanceConfig(str(d), depth, str(tmp_path / "o"), "RPM"), "fam_substrate_abund")
    t = _table(tmp_path / "o" / "fam_substrate_abund.out")
    assert set(t) == {"alpha-glucan", "starch"}
    assert t["starch"][1] == "gA"


def test_cgc_abund_skips_placeholders_and_stays_aligned(tmp_path):
    d = _dbcan_dir(tmp_path)
    depth = _depth(tmp_path, [("gA", 1000, 300), ("gB", 500, 100), ("gD", 100, 100)])
    run_abundance(AbundanceConfig(str(d), depth, str(tmp_path / "o"), "RPM"), "CGC_abund")
    t = _table(tmp_path / "o" / "CGC_abund.out")
    abund, seqids, seq_abund, fams = t["c1|CGC1"]
    assert seqids == "gA;gB" and fams == "GH13_21;TC|1.A.1"
    assert [float(x) for x in seq_abund.split(";")] == pytest.approx([600000.0, 200000.0])
    assert float(abund) == pytest.approx(400000.0)
    assert t["c1|CGC2"][1] == "gD"


def test_cgc_substrate_abund(tmp_path):
    d = _dbcan_dir(tmp_path)
    depth = _depth(tmp_path, [("gA", 1000, 300), ("gB", 500, 100), ("gD", 100, 100)])
    run_abundance(AbundanceConfig(str(d), depth, str(tmp_path / "o"), "RPM"), "CGC_substrate_abund")
    homo = _table(tmp_path / "o" / "CGC_substrate_PUL_homology.out")
    major = _table(tmp_path / "o" / "CGC_substrate_majority_voting.out")
    assert set(homo) == {"starch"} and homo["starch"][1] == "c1|CGC1"
    assert set(major) == {"alpha-glucan", "starch"}


def test_id_mismatch_fails(tmp_path):
    d = _dbcan_dir(tmp_path)
    depth = _depth(tmp_path, [("k141_1_1", 1000, 300), ("k141_1_2", 500, 100)])
    with pytest.raises(AbundanceError, match="None of the"):
        run_abundance(AbundanceConfig(str(d), depth, str(tmp_path / "o"), "RPM"), "CGC_abund")


def test_cli_exit_codes_and_output_dir(tmp_path, edge_bam):
    d = _dbcan_dir(tmp_path)
    _write_gff(tmp_path / "cgc.gff", "cgcgff")
    runner = CliRunner()
    depth = tmp_path / "cov" / "s.depth.txt"
    r = runner.invoke(cli, ["cal_coverage", "-i", edge_bam, "-g", str(tmp_path / "cgc.gff"), "-o", str(depth)])
    assert r.exit_code == 0, r.output
    assert _counts(depth)["gA"] == (900, 4)
    outdir = tmp_path / "abund"
    for fn in ("fam_abund", "fam_substrate_abund", "CGC_abund", "CGC_substrate_abund"):
        r = runner.invoke(cli, [fn, "-i", str(d), "-bt", str(depth), "-a", "rpkm", "-o", str(outdir)])
        assert r.exit_code == 0, (fn, r.output)
    assert sorted(os.listdir(outdir)) == sorted([
        "fam_abund.out", "subfam_abund.out", "EC_abund.out", "fam_substrate_abund.out", "CGC_abund.out",
        "CGC_substrate_PUL_homology.out", "CGC_substrate_majority_voting.out"])
    # failures exit non-zero
    bad = _depth(tmp_path, [("x", 10, 5)])
    assert runner.invoke(cli, ["fam_abund", "-i", str(d), "-bt", bad, "-o", str(outdir)]).exit_code != 0
    assert runner.invoke(cli, ["fam_abund", "-i", str(d), "-o", str(outdir)]).exit_code != 0
    assert runner.invoke(cli, ["cal_coverage", "-i", edge_bam, "-o", str(depth)]).exit_code != 0


@pytest.mark.parametrize("threads", [1, 2])
def test_rejected_alignments_are_counted_once(tmp_path, caplog, threads):
    """A rejected alignment that overlaps two genes (gB/gC) or two exons of one gene is tallied once."""
    import logging
    b = BamBuilder()
    b.add("mapq_two_genes", 1501, mapq=0)            # gB and gC
    b.add("mm_two_genes", 1505, mism=5)              # identity 0.95 < 0.98, gB and gC
    b.add("mapq_two_exons", 451, mapq=0)             # both exons of g2 (401-500, 551-700)
    b.add("dup_one_gene", 2401, flag=0x400)          # gD only
    b.add("ok", 2401)
    bam = b.write(tmp_path / "rej.bam")
    gff = tmp_path / "x.gff"
    rows = [(gid, s, e, st) for gid, s, e, st in GENES if gid != "gA"]
    with open(gff, "w") as fh:
        for gid, s, e, st in rows:
            fh.write(f"c1\tsim\tCDS\t{s}\t{e}\t.\t{st}\t0\tID={gid}\n")
        fh.write("c1\tsim\tCDS\t401\t500\t.\t+\t0\tID=g2\n")
        fh.write("c1\tsim\tCDS\t551\t700\t.\t+\t0\tID=g2\n")
    caplog.set_level(logging.INFO, logger="dbcan.utils.utils")
    out = tmp_path / "d.txt"
    cal_coverage(CoverageConfig(input=bam, gff=str(gff), output=str(out), threads=threads, overlap_base_ratio=0.05))
    assert _counts(out)["gD"] == (100, 1)
    line = next(r.getMessage() for r in caplog.records if r.getMessage().startswith("Alignments rejected"))
    assert line == "Alignments rejected (each alignment once): duplicate=1, low_identity=1, low_mapq=2"


def test_cgc_abund_fams_of_multi_role_and_signature_members(tmp_path):
    d = _dbcan_dir(tmp_path)
    (d / "cgc_standard_out.tsv").write_text(
        "CGC#\tGene Type\tContig ID\tProtein ID\tGene Start\tGene Stop\tGene Strand\tGene Annotation\n"
        "CGC1\tCAZyme\tc1\tgA\t101\t1000\t+\tCAZyme|GH32_e132+Sulfatase|S1_7\n"
        "CGC1\tSTP\tc1\tgB\t1001\t1600\t+\tSTP|HTH_1+STP|LysR_substrate\n"
        "CGC1\tnull\tc1\tgC\t1601\t2000\t+\tnull\n"
        "CGC1\tTC\tc1\tgD\t2401\t2500\t+\tTC|3.A.1.1.27+STP|SBP_bac_1\n"
    )
    depth = _depth(tmp_path, [("gA", 1000, 300), ("gB", 500, 100), ("gC", 1000, 0), ("gD", 100, 100)])
    run_abundance(AbundanceConfig(str(d), depth, str(tmp_path / "o"), "RPM"), "CGC_abund")
    _abund, seqids, seq_abund, fams = _table(tmp_path / "o" / "CGC_abund.out")["c1|CGC1"]
    assert seqids.split(";") == ["gA", "gB", "gC", "gD"]
    assert fams.split(";") == ["GH32_e132", "STP|HTH_1+LysR_substrate", "null", "TC|3.A.1.1.27"]
    assert len(seq_abund.split(";")) == 4


def test_cgc_substrate_abund_splits_pul_substrates(tmp_path):
    d = _dbcan_dir(tmp_path)
    (d / "substrate_prediction.tsv").write_text(
        "#cgcid\tPULID\tdbCAN-PUL substrate\tbitscore\tsignature pairs\tdbCAN-sub substrate\tdbCAN-sub substrate score\n"
        "c1|CGC1\tPUL1\tarabinan,galactan\t100\tCAZyme-CAZyme\tstarch\t2.0\n"
        "c1|CGC2\tPUL2\talpha-1,4-glucan;galactan\t100\tCAZyme-CAZyme\t\t\n"
    )
    depth = _depth(tmp_path, [("gA", 1000, 300), ("gB", 500, 100), ("gD", 100, 100)])
    run_abundance(AbundanceConfig(str(d), depth, str(tmp_path / "o"), "RPM"), "CGC_substrate_abund")
    homo = _table(tmp_path / "o" / "CGC_substrate_PUL_homology.out")
    assert set(homo) == {"arabinan", "galactan", "alpha-1,4-glucan"}
    assert homo["galactan"][1] == "c1|CGC1;c1|CGC2"


def test_missing_input_dir_names_the_path(tmp_path):
    depth = _depth(tmp_path, [("gA", 1000, 300)])
    missing = tmp_path / "no_such_run"
    r = CliRunner().invoke(cli, ["fam_abund", "-i", str(missing), "-bt", depth, "-o", str(tmp_path / "o")])
    assert r.exit_code == 2
    out = " ".join(r.output.replace("\u2502", " ").split())
    assert str(missing) in out and "/overview.tsv" not in out
