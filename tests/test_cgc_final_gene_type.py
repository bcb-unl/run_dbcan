"""CGC roles come from the final label (Gene Type) of each gene only.

Every gene of cgc.gff fills exactly one role: its Gene Type, chosen by the priority
CAZyme > TC > TF > STP > Sulfatase > Peptidase. Its other labels stay in the annotation text
but cannot fill a role, so e.g. a GH13 gene with an STP domain is a CAZyme and never the STP.
"""
import pandas as pd
import pytest

from dbcan.annotation.CGCFinder import CGCFinder
from dbcan.configs.cgcfinder_config import CGCFinderConfig
import dbcan.constants.cgcfinder_constants as C
from dbcan.IO.cgc_ids import gene_type, parse_roles
from dbcan.plot.cgc_viz import primary_gene_type

GH31 = ("pA", "CAZyme|GH31_e3")
GH13_STP = ("pB", "CAZyme|GH13_e1056+STP|SBP_bac_3")
STP = ("pC", "STP|SBP_bac_3")
TC_STP = ("pD", "TC|3.A.1.14.3+STP|Peripla_BP_2")
CE1_S09 = ("pE", "CAZyme|CE1_e5+Peptidase|S09.940")
TF = ("pF", "TF|CCPA_BACSU")


def _cgcs(tmp_path, genes, **cfg):
    """genes: (protein_id, CGC_annotation) in contig order, 1 kb apart; returns the CGCs as lists."""
    gff = "".join(
        f"c1\t.\tCDS\t{1000 * i + 1}\t{1000 * i + 900}\t.\t+\t.\tprotein_id={pid};CGC_annotation={ann}\n"
        for i, (pid, ann) in enumerate(genes)
    )
    (tmp_path / "cgc.gff").write_text(gff)
    finder = CGCFinder(CGCFinderConfig(output_dir=str(tmp_path), **cfg))
    assert finder.run()
    out = pd.read_csv(tmp_path / "cgc_standard_out.tsv", sep="\t", dtype=str)
    return out.groupby("CGC#")["Protein ID"].apply(list).tolist()


def test_stp_label_of_a_cazyme_gene_does_not_fill_the_stp_role(tmp_path):
    assert _cgcs(tmp_path, [GH13_STP, GH31], additional_genes=["STP"]) == []


def test_real_stp_gene_fills_the_stp_role(tmp_path):
    assert _cgcs(tmp_path, [STP, GH31], additional_genes=["STP"]) == [["pC", "pA"]]


def test_tc_plus_stp_gene_counts_only_as_tc(tmp_path):
    assert _cgcs(tmp_path, [TC_STP, GH31], additional_genes=["STP"]) == []
    assert _cgcs(tmp_path, [TC_STP, GH31], additional_genes=["TC"]) == [["pD", "pA"]]
    # TC and STP both required: the TC+STP gene gives only TC
    assert _cgcs(tmp_path, [TC_STP, GH31], additional_genes=["TC", "STP"]) == []
    assert _cgcs(tmp_path, [TC_STP, STP, GH31], additional_genes=["TC", "STP"]) == [["pD", "pC", "pA"]]


def test_cazyme_plus_peptidase_gene_counts_only_as_cazyme(tmp_path):
    assert _cgcs(tmp_path, [CE1_S09, GH31], additional_genes=["Peptidase"]) == []
    assert _cgcs(tmp_path, [CE1_S09, ("pP", "Peptidase|S09.940")], additional_genes=["Peptidase"]) == [["pE", "pP"]]


def test_default_rule_is_one_cazyme_and_one_tc(tmp_path):
    # CLI default --additional_genes TC; STP and TF genes are then null genes
    tc = dict(additional_genes=C.CGC_DEFAULT_ADDITIONAL_GENES)
    genes = [GH31, STP, TF, ("pT", "TC|2.A.1.1.1")]
    assert _cgcs(tmp_path, genes, **tc) == [["pA", "pC", "pF", "pT"]]
    assert _cgcs(tmp_path, [GH31, STP, TF, ("pN", "null"), ("pT", "TC|2.A.1.1.1")], **tc) == []
    assert _cgcs(tmp_path, [GH13_STP, GH31], **tc) == []


def test_additional_logic_any_and_min_categories(tmp_path):
    genes = [TC_STP, GH31, TF]
    cfg = dict(additional_genes=["TC", "TF", "STP"], additional_logic="any")
    assert _cgcs(tmp_path, genes, **cfg) == [["pD", "pA", "pF"]]
    # TC + TF are two categories; the STP label of pD is not a third
    assert _cgcs(tmp_path, genes, additional_min_categories=2, **cfg) == [["pD", "pA", "pF"]]
    assert _cgcs(tmp_path, genes, additional_min_categories=3, **cfg) == []
    assert _cgcs(tmp_path, [GH13_STP, GH31], **cfg) == []


def test_cazyme_pairs(tmp_path):
    # --additional_genes CAZyme alone: >= 2 CAZymes; the TC gene is a null gene
    assert _cgcs(tmp_path, [GH13_STP, GH31], additional_genes=["CAZyme"]) == [["pB", "pA"]]
    assert _cgcs(tmp_path, [GH31, ("pT", "TC|2.A.1.1.1")], additional_genes=["CAZyme"]) == []
    # CAZyme next to another type asks for one CAZyme more than --min_core_cazyme
    assert _cgcs(tmp_path, [GH31, ("pT", "TC|2.A.1.1.1")], additional_genes=["CAZyme", "TC"]) == []
    assert _cgcs(tmp_path, [GH31, GH13_STP, ("pT", "TC|2.A.1.1.1")],
                 additional_genes=["CAZyme", "TC"]) == [["pA", "pB", "pT"]]


def test_custom_type_in_additional_genes(tmp_path):
    assert _cgcs(tmp_path, [GH31, ("pR", "Other|tRNA")], additional_genes=["Other"]) == [["pA", "pR"]]
    assert _cgcs(tmp_path, [GH31, ("pR", "Other|tRNA")], additional_genes=["TC"]) == []


@pytest.mark.parametrize("annotation, expected", [
    ("CAZyme|GH13_e1056+STP|SBP_bac_3", "CAZyme"),
    ("TC|3.A.1.14.3+STP|Peripla_BP_2", "TC"),
    ("STP|TPR_1+TF|Response_reg", "TF"),
    ("STP|TCS_like", "STP"),              # substring 'TC' is not a transporter
    ("Peptidase|M38+Sulfatase|S1_8", "Sulfatase"),
    ("PEPTIDASE|M38+SULFATLAS|S1_8", "Sulfatase"),
    ("Other|tRNA", "Other"),
    ("null", "null"),
    ("", "null"),
])
def test_gene_type_is_the_highest_priority_label(annotation, expected):
    assert gene_type(annotation) == expected
    assert CGCFinder.get_gene_type(annotation) == expected
    if parse_roles(annotation):  # plots and substrate prediction use the same role
        assert parse_roles(annotation)[0][0] == expected == primary_gene_type(annotation)


def test_gene_type_column_is_the_role(tmp_path):
    _cgcs(tmp_path, [TC_STP, GH13_STP, CE1_S09], additional_genes=["TC"])
    out = pd.read_csv(tmp_path / "cgc_standard_out.tsv", sep="\t", dtype=str)
    assert out["Gene Type"].tolist() == ["TC", "CAZyme", "CAZyme"]
    assert out["Gene Annotation"].tolist()[0] == "TC|3.A.1.14.3+STP|Peripla_BP_2"


def test_additional_genes_spelling_is_normalized(tmp_path):
    # stp / SULFATLAS are matched like STP / Sulfatase (the Gene Type spelling)
    assert _cgcs(tmp_path, [STP, GH31], additional_genes=["stp"]) == [["pC", "pA"]]
    sulf = ("pS", "SULFATLAS|S1_8")
    assert _cgcs(tmp_path, [GH31, sulf], additional_genes=["sulfatlas"]) == [["pA", "pS"]]
    assert _cgcs(tmp_path, [GH31, sulf], additional_genes=["Sulfatase", "SULFATASE"]) == [["pA", "pS"]]
