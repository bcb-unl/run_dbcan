"""substrate_prediction reads --output_dir only; --mode/--input_raw_data are not required (#81)."""
from unittest.mock import patch

import pytest
from click.testing import CliRunner

from dbcan.main import cli


@pytest.fixture
def calls():
    recorded = {}
    with patch("dbcan.core.run_dbCAN_CGCFinder_substrate", side_effect=lambda cfg: recorded.setdefault("sub", cfg)), \
         patch("dbcan.core.run_dbcan_syn_plot", side_effect=lambda cfg: recorded.setdefault("syn", cfg)):
        yield recorded


def test_substrate_prediction_without_input_options(tmp_path, calls):
    result = CliRunner().invoke(
        cli, ["substrate_prediction", "--output_dir", str(tmp_path), "--db_dir", str(tmp_path)]
    )
    assert result.exit_code == 0, result.output
    assert calls["sub"].output_dir == str(tmp_path)
    assert calls["sub"].db_dir == str(tmp_path)


def test_substrate_prediction_still_accepts_legacy_input_options(tmp_path, calls):
    result = CliRunner().invoke(
        cli,
        ["substrate_prediction", "--mode", "protein", "--input_raw_data", "x.faa",
         "--output_dir", str(tmp_path), "--db_dir", str(tmp_path)],
    )
    assert result.exit_code == 0, result.output
    assert "sub" in calls
    # setup_logging() replaces the root handlers, so the warning is read from the captured stderr
    assert "not used by substrate_prediction" in (result.output + (result.stderr if result.stderr_bytes else ""))


def test_substrate_prediction_help_hides_input_options():
    result = CliRunner().invoke(cli, ["substrate_prediction", "--help"])
    assert result.exit_code == 0
    assert "--input_raw_data" not in result.output
    assert "--output_dir" in result.output


def test_easy_substrate_still_requires_input_raw_data(tmp_path):
    result = CliRunner().invoke(
        cli, ["easy_substrate", "--mode", "protein", "--output_dir", str(tmp_path), "--db_dir", str(tmp_path)]
    )
    assert result.exit_code != 0
    assert "input_raw_data" in result.output
