"""Coverage filters of the DIAMOND signature searches.

Full-length protein databases (TCDB, TF, SulfAtlas) filter on query coverage. The
peptidase database holds MEROPS peptidase units (domains, median ~210 aa), so it filters
on subject coverage: with --query-cover 35 a peptidase unit inside a long multi-domain
protein (e.g. 160 aa in a 1,000 aa protein, 16 %) was always rejected.
"""
import subprocess
from unittest.mock import patch

import pytest

from dbcan.annotation.diamond import PeptidaseDiamondProcessor, TCDBDiamondProcessor
from dbcan.configs.diamond_config import DiamondPeptidaseConfig, DiamondTCConfig


def _cmd(processor_cls, config_cls, tmp_path, **overrides):
    (tmp_path / "uniInput.faa").write_text(">p1\nMKV\n")
    cfg = config_cls(db_dir=str(tmp_path), output_dir=str(tmp_path), **overrides)
    (tmp_path / cfg.db_file).write_text("")
    seen = {}

    def fake_run(cmd, **kwargs):
        seen["cmd"] = cmd
        return subprocess.CompletedProcess(cmd, 0, "", "")

    with patch("dbcan.annotation.diamond.subprocess.run", side_effect=fake_run), \
         patch.object(processor_cls, "format_results", lambda self: None):
        processor_cls(cfg).run()
    return seen["cmd"]


def test_peptidase_filters_on_subject_coverage(tmp_path):
    cmd = _cmd(PeptidaseDiamondProcessor, DiamondPeptidaseConfig, tmp_path)
    assert "--subject-cover" in cmd and "--query-cover" not in cmd
    assert cmd[cmd.index("--subject-cover") + 1] == "35"
    assert "scovhsp" in cmd and "qcovhsp" not in cmd


def test_full_length_databases_keep_query_coverage(tmp_path):
    cmd = _cmd(TCDBDiamondProcessor, DiamondTCConfig, tmp_path)
    assert "--query-cover" in cmd and "--subject-cover" not in cmd
    assert "qcovhsp" in cmd


def test_invalid_coverage_mode_is_rejected(tmp_path):
    with pytest.raises(ValueError, match="coverage_mode"):
        _cmd(PeptidaseDiamondProcessor, DiamondPeptidaseConfig, tmp_path, coverage_mode="both")


def test_sulfatase_filters_on_query_and_subject_coverage(tmp_path):
    """SulfAtlas entries are full-length proteins, many of them sulfatases fused to other domains
    (GH2, GH32, ...). A GH32 protein that only aligns to the GH32 part of such an entry (subject
    coverage ~30 %) must not take its S1_7 label, so the hit must also cover most of the entry."""
    from dbcan.annotation.diamond import SulfatlasDiamondProcessor
    from dbcan.configs.diamond_config import DiamondSulfataseConfig
    cmd = _cmd(SulfatlasDiamondProcessor, DiamondSulfataseConfig, tmp_path)
    assert cmd[cmd.index("--query-cover") + 1] == "35"
    assert cmd[cmd.index("--subject-cover") + 1] == "70"
    off = _cmd(SulfatlasDiamondProcessor, DiamondSulfataseConfig, tmp_path, subject_coverage_threshold=None)
    assert "--subject-cover" not in off


def test_sulfatase_subject_coverage_option_reaches_the_config():
    from dbcan.configs.diamond_config import DiamondSulfataseConfig
    from dbcan.parameter import create_config
    cfg = create_config(DiamondSulfataseConfig, namespace="sulfatase", db_dir="d", output_dir="o",
                        subject_coverage_threshold_sulfatase=80.0)
    assert cfg.subject_coverage_threshold == 80.0
