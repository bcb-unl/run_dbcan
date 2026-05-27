import pytest
import pandas as pd

pydeseq2 = pytest.importorskip("pydeseq2")

from dbcan.configs.expression_config import ExpressionConfig
from dbcan.expression.deseq2_runner import run_deseq2
from dbcan.expression.samplesheet import Samplesheet, SampleRecord


@pytest.mark.expression
def test_run_deseq2_synthetic(tmp_path):
    counts = pd.DataFrame(
        {
            "g1": [100, 120, 110, 10, 12, 11],
            "g2": [50, 55, 52, 48, 51, 49],
            "g3": [5, 6, 5, 5, 6, 5],
        },
        index=["C1", "C2", "C3", "T1", "T2", "T3"],
    )
    sheet = Samplesheet(
        samples=[
            SampleRecord(f"C{i}", "control", bam="x") for i in range(1, 4)
        ]
        + [SampleRecord(f"T{i}", "treat", bam="x") for i in range(1, 4)],
        path=tmp_path / "s.tsv",
    )
    cfg = ExpressionConfig(output_dir=str(tmp_path), design="~condition", alpha=0.05, lfc_threshold=0.5)
    gene_deg = run_deseq2(cfg, sheet, counts)
    assert "protein_id" in gene_deg.columns
    assert "is_deg" in gene_deg.columns
    assert (tmp_path / "gene_deg.tsv").is_file()
    assert (tmp_path / "DEG.tsv").is_file()
