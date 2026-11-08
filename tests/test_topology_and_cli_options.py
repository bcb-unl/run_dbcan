"""Topology (SignalP/DeepTMHMM) runs per requested tool and fails loudly; deprecated CLI options; --threads."""
import sys
import textwrap
from unittest.mock import patch

import pandas as pd
import pytest
from click.testing import CliRunner

from dbcan.configs.base_config import GeneralConfig
from dbcan.configs.signalp_tmhmm_config import SignalPTMHMMConfig
from dbcan.core import DbcanStepError, run_dbCAN_topology_annotation
from dbcan.main import cli
from dbcan.parameter import create_config

PREDICT_OK = textwrap.dedent('''
    import argparse, os
    p = argparse.ArgumentParser(); p.add_argument("--fasta"); p.add_argument("--output-dir")
    a = p.parse_args(); os.makedirs(a.output_dir, exist_ok=True)
    ids = [l[1:].split()[0] for l in open(a.fasta) if l.startswith(">")]
    with open(os.path.join(a.output_dir, "predicted_topologies.3line"), "w") as fh:
        for i in ids:
            fh.write(f">{i} | TM\\nMKV\\nMMM\\n")
''')


@pytest.fixture
def outdir(tmp_path):
    pd.DataFrame({"Gene ID": ["p1", "p2"], "Recommend Results": ["GH5_2", "-"]}).to_csv(
        tmp_path / "overview.tsv", sep="\t", index=False)
    (tmp_path / "uniInput.faa").write_text(">p1\nMKVLA\n>p2\nMKKKK\n")
    return tmp_path


def _deeptmhmm_dir(tmp_path, script):
    d = tmp_path / "deeptmhmm"
    d.mkdir()
    (d / "predict.py").write_text(script)
    return d


def test_deeptmhmm_runs_without_signalp(outdir, tmp_path):
    d = _deeptmhmm_dir(tmp_path, PREDICT_OK)
    cfg = SignalPTMHMMConfig(output_dir=str(outdir), run_deeptmhmm=True, deeptmhmm_dir=str(d),
                             deeptmhmm_python=sys.executable)
    run_dbCAN_topology_annotation(cfg)
    ov = pd.read_csv(outdir / "overview.tsv", sep="\t", dtype=str).set_index("Gene ID")
    assert ov.loc["p1", "DeepTMHMM"] == "TM"
    assert ov.loc["p2", "DeepTMHMM"] == "-"
    assert "SignalP" not in ov.columns


def test_failed_deeptmhmm_raises(outdir, tmp_path):
    d = _deeptmhmm_dir(tmp_path, "import sys; sys.exit(3)\n")
    cfg = SignalPTMHMMConfig(output_dir=str(outdir), run_deeptmhmm=True, deeptmhmm_dir=str(d),
                             deeptmhmm_python=sys.executable)
    with pytest.raises(DbcanStepError, match="DeepTMHMM"):
        run_dbCAN_topology_annotation(cfg)


def test_missing_signalp_raises_but_keeps_deeptmhmm_columns(outdir, tmp_path):
    d = _deeptmhmm_dir(tmp_path, PREDICT_OK)
    cfg = SignalPTMHMMConfig(output_dir=str(outdir), run_signalp=True, run_deeptmhmm=True, deeptmhmm_dir=str(d),
                             deeptmhmm_python=sys.executable)
    import subprocess
    real_run = subprocess.run

    def fake(cmd, *a, **kw):
        if cmd[0] == "signalp6":
            raise FileNotFoundError("signalp6")
        return real_run(cmd, *a, **kw)

    with patch("dbcan.annotation.signalp_tmhmm.subprocess.run", side_effect=fake):
        with pytest.raises(DbcanStepError, match="SignalP6"):
            run_dbCAN_topology_annotation(cfg)
    ov = pd.read_csv(outdir / "overview.tsv", sep="\t", dtype=str).set_index("Gene ID")
    assert ov.loc["p1", "DeepTMHMM"] == "TM"


def test_missing_deeptmhmm_dir_raises(outdir):
    cfg = SignalPTMHMMConfig(output_dir=str(outdir), run_deeptmhmm=True)
    with pytest.raises(DbcanStepError, match="deeptmhmm_dir"):
        run_dbCAN_topology_annotation(cfg)


def test_cazyme_annotation_runs_topology_for_deeptmhmm_only(tmp_path):
    with patch("dbcan.core.run_dbCAN_input_process"), \
         patch("dbcan.core.run_dbCAN_CAZyme_annotation"), \
         patch("dbcan.core.run_dbCAN_topology_annotation") as topo:
        result = CliRunner().invoke(cli, [
            "CAZyme_annotation", "--input_raw_data", "x.faa", "--mode", "protein", "--output_dir", str(tmp_path),
            "--db_dir", str(tmp_path), "--run_deeptmhmm", "--deeptmhmm_dir", str(tmp_path)])
    assert result.exit_code == 0, result.output
    cfg = topo.call_args.args[0]
    assert cfg.run_deeptmhmm and not cfg.run_signalp


def test_cazyme_annotation_exits_nonzero_when_topology_fails(tmp_path):
    with patch("dbcan.core.run_dbCAN_input_process"), \
         patch("dbcan.core.run_dbCAN_CAZyme_annotation"), \
         patch("dbcan.core.run_dbCAN_topology_annotation", side_effect=DbcanStepError("SignalP6 failed")):
        result = CliRunner().invoke(cli, [
            "CAZyme_annotation", "--input_raw_data", "x.faa", "--mode", "protein", "--output_dir", str(tmp_path),
            "--db_dir", str(tmp_path), "--run_signalp"])
    assert result.exit_code != 0


# ---- deprecated / ignored options ---------------------------------------------------------------

def test_ignored_options_are_accepted_hidden_and_warned(tmp_path):
    import dbcan.parameter as P
    with patch("dbcan.core.run_dbCAN_CGCFinder") as finder, patch.object(P.logger, "warning") as warn:
        result = CliRunner().invoke(cli, ["cgc_finder", "--output_dir", str(tmp_path), "--feature_type", "gene"])
    assert result.exit_code == 0, result.output
    assert finder.called
    assert [c.args[1] for c in warn.call_args_list] == ["--feature_type"]
    assert "--feature_type" not in CliRunner().invoke(cli, ["cgc_finder", "--help"]).output

    with patch("dbcan.core.run_dbCAN_input_process"), patch("dbcan.core.run_dbCAN_CAZyme_annotation"), \
         patch("dbcan.core.run_dbCAN_topology_annotation"), patch.object(P.logger, "warning") as warn:
        result = CliRunner().invoke(cli, [
            "CAZyme_annotation", "--input_raw_data", "x.faa", "--mode", "protein", "--output_dir", str(tmp_path),
            "--db_dir", str(tmp_path), "--force_topology"])
    assert result.exit_code == 0, result.output
    assert [c.args[1] for c in warn.call_args_list] == ["--force_topology"]

    helptext = CliRunner().invoke(cli, ["Pfam_null_cgc", "--help"]).output
    assert "--run_pfam" not in helptext


def test_default_option_values_do_not_warn(tmp_path):
    import dbcan.parameter as P
    with patch("dbcan.core.run_dbCAN_CGCFinder"), patch.object(P.logger, "warning") as warn:
        result = CliRunner().invoke(cli, ["cgc_finder", "--output_dir", str(tmp_path)])
    assert result.exit_code == 0, result.output
    assert not warn.called


# ---- --threads reaches pyrodigal ----------------------------------------------------------------

def test_threads_reach_gene_prediction(tmp_path):
    from dbcan.IO.fasta import ProkProcessor
    fna = tmp_path / "x.fna"
    fna.write_text(">c\nATGAAATAA\n")
    cfg = create_config(GeneralConfig, input_raw_data=str(fna), output_dir=str(tmp_path), mode="prok", threads=3,
                        db_dir=str(tmp_path))
    assert cfg.threads == 3
    assert ProkProcessor(cfg).threads == 3
    with patch("dbcan.core.run_dbCAN_input_process") as inp, patch("dbcan.core.run_dbCAN_CAZyme_annotation"):
        result = CliRunner().invoke(cli, [
            "CAZyme_annotation", "--input_raw_data", str(fna), "--mode", "prok", "--output_dir", str(tmp_path),
            "--db_dir", str(tmp_path), "--threads", "2"])
    assert result.exit_code == 0, result.output
    assert inp.call_args.args[0].threads == 2
