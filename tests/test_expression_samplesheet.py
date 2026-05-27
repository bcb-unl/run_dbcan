import pytest
from pathlib import Path
import tempfile

from dbcan.expression.samplesheet import load_samplesheet


def test_load_samplesheet_bam_only():
    content = "sample_id\tcondition\tbam\nS1\tcontrol\t/path/a.bam\nS2\ttreat\t/path/b.bam\n"
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        f.write(content)
        path = f.name
    sheet = load_samplesheet(path)
    assert len(sheet.samples) == 2
    assert sheet.samples[0].sample_id == "S1"
    Path(path).unlink()


def test_samplesheet_requires_alignment_input():
    content = "sample_id\tcondition\nS1\tcontrol\nS2\ttreat\n"
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        f.write(content)
        path = f.name
    with pytest.raises(ValueError, match="bam or r1"):
        load_samplesheet(path)
    Path(path).unlink()


def test_samplesheet_header_with_trailing_spaces():
    content = "sample_id\tcondition   \tbam\nS1\tcontrol\ta.bam\nS2\ttreat\tb2.bam\n"
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        f.write(content)
        path = f.name
    sheet = load_samplesheet(path)
    assert len(sheet.samples) == 2
    assert sheet.samples[0].condition == "control"
    Path(path).unlink()


def test_duplicate_sample_id():
    content = "sample_id\tcondition\tbam\nS1\tc\ta.bam\nS1\tt\ta2.bam\n"
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        f.write(content)
        path = f.name
    with pytest.raises(ValueError, match="Duplicate"):
        load_samplesheet(path)
    Path(path).unlink()
