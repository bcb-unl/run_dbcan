import pytest

from dbcan.expression.deseq2_runner import validate_deseq2_design
from dbcan.expression.samplesheet import SampleRecord, Samplesheet


def _sheet(*rows):
    samples = [
        SampleRecord(sid, cond, bam="x")
        for sid, cond in rows
    ]
    return Samplesheet(samples=samples, path="s.tsv")


def test_validate_rejects_one_sample_per_condition():
    sheet = _sheet(
        ("S1", "control"),
        ("S2", "test1"),
        ("S3", "test2"),
        ("S4", "test3"),
    )
    with pytest.raises(ValueError, match="4 sample"):
        validate_deseq2_design(sheet, "~condition")


def test_validate_accepts_replicated_conditions():
    sheet = _sheet(
        ("C1", "control"),
        ("C2", "control"),
        ("T1", "treatment"),
        ("T2", "treatment"),
    )
    validate_deseq2_design(sheet, "~condition")
