import pandas as pd

from dbcan.plot.expression_plots import build_cgc_gene_matrix_with_labels
from dbcan.utils.plots import CGC, cgc_standard_line


def _minimal_cgc():
    genes = []
    for i, (s, e, pid) in enumerate([(100, 200, "g1"), (300, 400, "g2")]):
        line = ["CGC1", "CAZyme", "scaf", pid, s, e, "+", f"CAZyme|GH{i+1}"]
        g = cgc_standard_line(line)
        g.annotation = line[-1]
        genes.append(g)
    return CGC(genes)


def test_heatmap_matrix_condition_rows(tmp_path):
    cgc = _minimal_cgc()
    counts = pd.DataFrame([
        {"cgcid": "scaf|CGC1", "protein_id": "g1", "sample_id": "S1", "count": 10, "tpm": 1.0},
        {"cgcid": "scaf|CGC1", "protein_id": "g1", "sample_id": "S2", "count": 20, "tpm": 2.0},
        {"cgcid": "scaf|CGC1", "protein_id": "g2", "sample_id": "S1", "count": 5, "tpm": 0.5},
        {"cgcid": "scaf|CGC1", "protein_id": "g2", "sample_id": "S2", "count": 8, "tpm": 0.8},
    ])
    cgc.ID = "scaf|CGC1"
    meta = pd.DataFrame(
        {"condition": ["control", "treat"]},
        index=["S1", "S2"],
    )
    expr = tmp_path / "expr"
    expr.mkdir()
    meta.reset_index().rename(columns={"index": "sample_id"}).to_csv(
        expr / "sample_metadata.tsv", sep="\t", index=False
    )
    mat, labels = build_cgc_gene_matrix_with_labels(
        cgc, counts, meta, "log2_tpm", expr, pd.DataFrame(), str(tmp_path), row_by="condition"
    )
    assert mat.shape[0] == 2
    assert mat.shape[1] == 2
    assert len(labels) == 2
