"""CGC calling with multi-label genes.

1. Same region: a domain-level signature (MEROPS peptidase unit, STP/TF HMM domain) lying on
   the gene's own CAZyme domain is one domain classified twice (e.g. CE1 vs S09 on an
   alpha/beta-hydrolase) and is not an extra signature.
2. One gene, one role: the core CAZyme(s) and each required additional category must be
   provided by different genes, so one multi-label gene cannot satisfy the CGC rule alone.
Annotations are taken from real E. coli K-12 hits (NP_414889.1 CE1 21-271 / S09.940 4-274).
"""
import pandas as pd
import pytest

from dbcan.annotation.CGCFinder import CGCFinder
from dbcan.configs.cgcfinder_config import CGCFinderConfig

OVERVIEW_HEADER = "Gene ID\tEC#\tdbCAN_hmm\tdbCAN_sub\tDIAMOND\t#ofTools\tRecommend Results\tSubstrate\n"
SIG_HEADER = ("Annotate Name\tAnnotate Length\tTarget Name\tTarget Length\ti-Evalue\tAnnotate From\tAnnotate To\t"
              "Target From\tTarget To\tCoverage\tAnnotate File Name\tType\tSource\n")


def _setup(tmp_path, genes, overview_rows, sig_rows):
    """genes: list of (protein_id, CGC_annotation); consecutive genes 1 kb apart on contig c1."""
    gff = "".join(
        f"c1\t.\tCDS\t{1000 * i + 1}\t{1000 * i + 900}\t.\t+\t.\tprotein_id={pid};CGC_annotation={ann}\n"
        for i, (pid, ann) in enumerate(genes)
    )
    (tmp_path / "cgc.gff").write_text(gff)
    (tmp_path / "overview.tsv").write_text(OVERVIEW_HEADER + "".join(overview_rows))
    (tmp_path / "total_cgc_info.tsv").write_text(SIG_HEADER + "".join(sig_rows))


def _cgcs(tmp_path, **cfg):
    finder = CGCFinder(CGCFinderConfig(output_dir=str(tmp_path), **cfg))
    assert finder.run()
    out = pd.read_csv(tmp_path / "cgc_standard_out.tsv", sep="\t", dtype=str)
    return out.groupby("CGC#")["Protein ID"].apply(list).tolist()


CE1_S09 = ("pX", "CAZyme|CE1_e5+Peptidase|S09.940")
CE1_ROW = "pX\t-\tCE1(21-271)\tCE1_e5(20-270)\tCE1\t3\tCE1_e5\t-\n"
S09_SAME_REGION = "S09.940\t270\tpX\t277\t4.7e-205\t1\t270\t4\t274\t100\tPeptidase\tPeptidase\tPeptidase\n"
GH5_ROW = "pY\t-\tGH5(30-300)\tGH5_e1(30-300)\tGH5\t3\tGH5_e1\t-\n"


def test_peptidase_on_the_cazyme_domain_is_not_an_extra_signature(tmp_path):
    _setup(tmp_path, [CE1_S09, ("pY", "CAZyme|GH5_e1")], [CE1_ROW, GH5_ROW], [S09_SAME_REGION])
    # old behaviour: pX counted as CAZyme and as Peptidase -> CGC of two CAZymes "with a peptidase"
    assert _cgcs(tmp_path, additional_genes=["Peptidase"]) == []


def test_peptidase_outside_the_cazyme_domain_still_counts(tmp_path):
    # a 600 aa protein: CE1 domain 21-271, peptidase unit 320-590 (a real two-domain protein)
    s09_elsewhere = "S09.940\t270\tpX\t600\t1e-80\t1\t270\t320\t590\t100\tPeptidase\tPeptidase\tPeptidase\n"
    _setup(tmp_path, [CE1_S09, ("pY", "CAZyme|GH5_e1")], [CE1_ROW, GH5_ROW], [s09_elsewhere])
    assert _cgcs(tmp_path, additional_genes=["Peptidase"]) == [["pX", "pY"]]


def test_full_length_tc_label_is_never_treated_as_same_region(tmp_path):
    """TC comes from full-length TCDB proteins: no coordinate comparison, it stays a TC label."""
    fusion = ("pX", "CAZyme|GH13_31+TC|8.A.9.1.1")
    tc_row = "8.A.9.1.1\t677\tpX\t553\t2.1e-60\t109\t596\t4\t516\t92.8\tTC\tTC\tTC\n"
    ov = ["pX\t-\tGH13_31(37-376)\tGH13_e119(27-378)\tGH13_31\t3\tGH13_31\t-\n", GH5_ROW]
    _setup(tmp_path, [fusion, ("pY", "CAZyme|GH5_e1")], ov, [tc_row])
    # pY is the core CAZyme and pX provides TC: two genes, two roles
    assert _cgcs(tmp_path, additional_genes=["TC"]) == [["pX", "pY"]]


def test_one_gene_cannot_be_both_the_cazyme_and_the_required_tc(tmp_path):
    fusion = ("pX", "CAZyme|GT2_e1+TC|4.D.3.1.1")
    tc_row = "4.D.3.1.1\t870\tpX\t870\t0.0\t1\t870\t1\t870\t99\tTC\tTC\tTC\n"
    ov = ["pX\t-\tGT2(140-330)\tGT2_e1(140-330)\tGT2\t3\tGT2_e1\t-\n"]
    tf_row = "CCPA_BACSU\t334\tpZ\t316\t5e-34\t5\t333\t3\t316\t99\tprodoric\tTF\tTF_prok\n"
    _setup(tmp_path, [fusion, ("pZ", "TF|CCPA_BACSU")], ov, [tc_row, tf_row])
    # 'all' of TC + TF: pX would have to be the CAZyme and the TC at the same time -> no CGC
    assert _cgcs(tmp_path, additional_genes=["TC", "TF"], additional_logic="all") == []
    # 'any' 1 category: pX is the CAZyme, pZ the TF -> valid
    assert _cgcs(tmp_path, additional_genes=["TC", "TF"], additional_logic="any") == [["pX", "pZ"]]


def test_same_region_needs_every_hit_of_the_category_on_the_domain(tmp_path):
    """Two STP domains: one on the CAZyme domain, one elsewhere -> STP still counts."""
    gene = ("pX", "CAZyme|CE1_e5+STP|A+STP|B")
    stp_on = "A\t100\tpX\t600\t1e-20\t1\t100\t30\t130\t0.9\tSTP\tSTP\tSTP\n"
    stp_off = "B\t100\tpX\t600\t1e-20\t1\t100\t400\t500\t0.9\tSTP\tSTP\tSTP\n"
    _setup(tmp_path, [gene, ("pY", "CAZyme|GH5_e1")], [CE1_ROW, GH5_ROW], [stp_on, stp_off])
    assert _cgcs(tmp_path, additional_genes=["STP"]) == [["pX", "pY"]]


@pytest.mark.parametrize("annotation, expected", [
    ("CAZyme|CBM34_e15|GH13_20", {"CAZyme"}),
    ("STP|TPR_1+TF|Response_reg", {"STP", "TF"}),
    ("STP|TCS_like", {"STP"}),          # substring 'TC' must not count as a transporter
    ("null", set()),
])
def test_annotation_categories_are_parsed_exactly(annotation, expected):
    assert CGCFinder.annotation_categories(annotation) == expected
