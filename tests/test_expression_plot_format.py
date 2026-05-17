import pandas as pd

from dbcan.plot.expression_plots import format_gene_label


def test_format_cazyme_annotation():
    assert format_gene_label("CAZyme", "CAZyme|GH24_e112", "g1") == "GH24"


def test_format_tc_annotation():
    label = format_gene_label("TC", "TC|1.A.115.1.3", "g2")
    assert label.startswith("TC:")
    assert "1.A.115" in label


def test_format_null_from_overview():
    row = pd.Series({
        "Recommend Results": "GH18_e325",
        "dbCAN_sub": "-",
        "DIAMOND": "-",
    })
    assert format_gene_label("null", "null", "g3", row) == "GH18"


def test_format_pfam_from_overview():
    row = pd.Series({
        "Recommend Results": "-",
        "dbCAN_sub": "-",
        "DIAMOND": "PF01234",
    })
    assert format_gene_label("null", None, "g4", row) == "PF01234"
