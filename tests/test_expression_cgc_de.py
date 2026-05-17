import pandas as pd

from dbcan.expression.cgc_de import _cgc_de_rule_passes, build_cgc_pseudobulk_matrix


def test_cgc_de_rules():
    assert _cgc_de_rule_passes(1, 5, "any")
    assert not _cgc_de_rule_passes(0, 5, "any")
    assert _cgc_de_rule_passes(3, 5, "majority")
    assert not _cgc_de_rule_passes(2, 5, "majority")
    assert _cgc_de_rule_passes(5, 5, "all")
    assert not _cgc_de_rule_passes(4, 5, "all")


def test_pseudobulk_matrix():
    count_matrix = pd.DataFrame(
        {"g1": [10, 20], "g2": [5, 15]},
        index=["S1", "S2"],
    )

    class Rec:
        def __init__(self, seqid):
            self.seqid = seqid

    cgcid2records = {"cgcA": [Rec("g1"), Rec("g2")], "cgcB": [Rec("g2")]}
    out = build_cgc_pseudobulk_matrix(count_matrix, cgcid2records)
    assert out.loc["S1", "cgcA"] == 15
    assert out.loc["S2", "cgcB"] == 15
