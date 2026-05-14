"""Tests for Prodigal GFF streaming annotation vs BCBio GFF.parse path."""
from pathlib import Path

import pytest

from dbcan.IO.gff import ProdigalProcessor, annotate_prodigal_gff_streaming
from dbcan.configs.base_config import GFFConfig
import dbcan.constants.gff_constants as G


MINIMAL_PRODIGAL_GFF = """##gff-version 3
##sequence-region c1 1 5000
c1\tProdigal_v2.6.3\tCDS\t2\t100\t.\t+\t0\tID=p1;partial=00;
c1\tProdigal_v2.6.3\tCDS\t200\t350\t.\t-\t0\tID=p2;partial=00;
"""


def _parse_path_output(path: Path) -> list[str]:
    return path.read_text(encoding="utf-8").strip().splitlines()


def test_annotate_prodigal_gff_streaming_cgc_lookup(tmp_path):
    inp = tmp_path / "in.gff"
    out = tmp_path / "out.gff"
    inp.write_text(MINIMAL_PRODIGAL_GFF, encoding="utf-8")
    cgc = {
        "p1": {G.GFF_CGC_ANNOTATION_COL: "TC|foo"},
        "p2": {G.GFF_CGC_ANNOTATION_COL: G.GFF_NULL_ANNOTATION},
    }
    n = annotate_prodigal_gff_streaming(inp, out, cgc)
    assert n == 2
    lines = _parse_path_output(out)
    assert any("p1" in ln and "TC|foo" in ln for ln in lines)
    assert sum(1 for ln in lines if "p2" in ln) == 1
    assert all("\t.\tCDS\t" in ln for ln in lines)


def test_prodigal_streaming_output_matches_expected_lines(tmp_path):
    """Golden output for minimal GFF (same as write_gff / BCBio semantics for these rows)."""
    inp = tmp_path / "in.gff"
    inp.write_text(MINIMAL_PRODIGAL_GFF, encoding="utf-8")
    cgc_data = {
        "p1": {G.GFF_CGC_ANNOTATION_COL: "GH1|test"},
        "p2": {G.GFF_CGC_ANNOTATION_COL: "TC|x"},
    }
    out = tmp_path / "stream.tsv"
    annotate_prodigal_gff_streaming(inp, out, cgc_data)
    lines = set(_parse_path_output(out))
    assert lines == {
        "c1\t.\tCDS\t2\t100\t.\t+\t.\tprotein_id=p1;CGC_annotation=GH1|test",
        "c1\t.\tCDS\t200\t350\t.\t-\t.\tprotein_id=p2;CGC_annotation=TC|x",
    }


def test_use_streaming_path_auto_threshold(tmp_path):
    small = tmp_path / "small.gff"
    small.write_text(MINIMAL_PRODIGAL_GFF, encoding="utf-8")
    cfg = GFFConfig(
        output_dir=str(tmp_path),
        input_gff=str(small),
        gff_type=G.GFF_FORMAT_PRODIGAL,
        prodigal_gff_streaming="auto",
        prodigal_streaming_threshold_mb=50,
    )
    proc = ProdigalProcessor(cfg)
    assert proc._use_streaming_path(small) is False

    cfg_zero = GFFConfig(
        output_dir=str(tmp_path),
        input_gff=str(small),
        gff_type=G.GFF_FORMAT_PRODIGAL,
        prodigal_gff_streaming="auto",
        prodigal_streaming_threshold_mb=0,
    )
    proc_zero = ProdigalProcessor(cfg_zero)
    assert proc_zero._use_streaming_path(small) is True


def test_use_streaming_path_on_off(tmp_path):
    p = tmp_path / "f.gff"
    p.write_text("x", encoding="utf-8")
    for mode, expected in (("on", True), ("off", False), ("ON", True)):
        cfg = GFFConfig(
            output_dir=str(tmp_path),
            input_gff=str(p),
            gff_type=G.GFF_FORMAT_PRODIGAL,
            prodigal_gff_streaming=mode,
        )
        assert ProdigalProcessor(cfg)._use_streaming_path(p) is expected
