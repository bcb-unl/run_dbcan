"""cgc_finder: a crash while clustering fails the step; a run that finds no CGC still succeeds."""
import pandas as pd
import pytest

from dbcan.annotation.CGCFinder import CGCFinder
from dbcan.configs.cgcfinder_config import CGCFinderConfig
from dbcan.core import DbcanStepError, run_dbCAN_CGCFinder


def _setup(tmp_path, layout):
    """layout: string, one char per gene 1 kb apart on contig c1: C = CAZyme, T = TC, . = null."""
    ann = {"C": "CAZyme|GH5_e1", "T": "TC|2.A.1.1.1", ".": "null"}
    gff = "".join(
        f"c1\t.\tCDS\t{1000 * i + 1}\t{1000 * i + 900}\t.\t+\t.\tprotein_id=g{i};CGC_annotation={ann[c]}\n"
        for i, c in enumerate(layout)
    )
    (tmp_path / "cgc.gff").write_text(gff)


def _config(tmp_path):
    return CGCFinderConfig(output_dir=str(tmp_path), additional_genes=["TC"])


def _crash(*args, **kwargs):
    raise RuntimeError("boom")


def test_clustering_crash_fails_without_writing_results(tmp_path, monkeypatch):
    _setup(tmp_path, "CT...CT...")
    monkeypatch.setattr(CGCFinder, "_resolve_extensions", _crash)
    assert CGCFinder(_config(tmp_path)).run() is False
    assert not (tmp_path / "cgc_standard_out.tsv").exists()
    assert not (tmp_path / "cgc_standard_out_summary.tsv").exists()


def test_clustering_crash_raises_step_error(tmp_path, monkeypatch):
    _setup(tmp_path, "CT...CT...")
    monkeypatch.setattr(CGCFinder, "validate_cluster", _crash)
    with pytest.raises(DbcanStepError):
        run_dbCAN_CGCFinder(_config(tmp_path))
    assert not (tmp_path / "cgc_standard_out.tsv").exists()


def test_no_cgc_writes_empty_tables_and_succeeds(tmp_path):
    _setup(tmp_path, "C...T.....")  # CAZyme and TC 3 null genes apart: no CGC
    run_dbCAN_CGCFinder(_config(tmp_path))
    out = pd.read_csv(tmp_path / "cgc_standard_out.tsv", sep="\t")
    summary = pd.read_csv(tmp_path / "cgc_standard_out_summary.tsv", sep="\t")
    assert out.empty and len(out.columns) > 0
    assert summary.empty and "CGC#" in summary.columns


def test_cgcs_found_still_written(tmp_path):
    _setup(tmp_path, "CT...CT...")
    assert CGCFinder(_config(tmp_path)).run() is True
    out = pd.read_csv(tmp_path / "cgc_standard_out.tsv", sep="\t", dtype=str)
    assert list(out["CGC#"].unique()) == ["CGC1", "CGC2"]
