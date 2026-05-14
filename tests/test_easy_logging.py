"""Unit tests for easy_* composite commands (logging propagation, etc.)."""
from pathlib import Path
from unittest.mock import patch

import pytest
from click.testing import CliRunner

from dbcan.main import cli

TEST_ROOT = Path(__file__).parent
DATA_ROOT = TEST_ROOT / "_data"
TEST_PROTEIN = DATA_ROOT / "EscheriaColiK12MG1655.faa"


@pytest.fixture
def runner():
    return CliRunner()


def test_easy_cgc_propagates_log_level_to_subcommands(runner, tmp_path):
    """Sub-invokes must receive log_level/log_file/verbose so setup_logging is not reset to WARNING."""
    assert TEST_PROTEIN.exists(), f"Missing fixture {TEST_PROTEIN}"
    captured = []

    def capture_invoke_subset(ctx, cmd, all_kwargs):
        captured.append(
            {
                "cmd": cmd.name,
                "log_level": all_kwargs.get("log_level"),
                "verbose": all_kwargs.get("verbose"),
                "log_file": all_kwargs.get("log_file"),
            }
        )

    out_dir = tmp_path / "out"
    out_dir.mkdir()
    db_dir = tmp_path / "db"
    db_dir.mkdir()

    with patch("dbcan.main._invoke_subset", side_effect=capture_invoke_subset):
        result = runner.invoke(
            cli,
            [
                "easy_CGC",
                "--log-level",
                "DEBUG",
                "--input_raw_data",
                str(TEST_PROTEIN),
                "--output_dir",
                str(out_dir),
                "--mode",
                "prok",
                "--db_dir",
                str(db_dir),
            ],
        )

    assert result.exit_code == 0, result.output
    assert len(captured) == 3
    for row in captured:
        assert row["log_level"] == "DEBUG"
        assert row["verbose"] is False
        assert row["log_file"] is None

    assert [c["cmd"] for c in captured] == [
        "CAZyme_annotation",
        "gff_process",
        "cgc_finder",
    ]


def test_easy_substrate_propagates_log_level_to_subcommands(runner, tmp_path):
    assert TEST_PROTEIN.exists(), f"Missing fixture {TEST_PROTEIN}"
    captured = []

    def capture_invoke_subset(ctx, cmd, all_kwargs):
        captured.append(
            {
                "cmd": cmd.name,
                "log_level": all_kwargs.get("log_level"),
                "verbose": all_kwargs.get("verbose"),
            }
        )

    out_dir = tmp_path / "out"
    out_dir.mkdir()
    db_dir = tmp_path / "db"
    db_dir.mkdir()

    with patch("dbcan.main._invoke_subset", side_effect=capture_invoke_subset):
        result = runner.invoke(
            cli,
            [
                "easy_substrate",
                "--log-level",
                "DEBUG",
                "--input_raw_data",
                str(TEST_PROTEIN),
                "--output_dir",
                str(out_dir),
                "--mode",
                "prok",
                "--db_dir",
                str(db_dir),
            ],
        )

    assert result.exit_code == 0, result.output
    assert len(captured) == 4
    for row in captured:
        assert row["log_level"] == "DEBUG"
        assert row["verbose"] is False

    assert [c["cmd"] for c in captured] == [
        "CAZyme_annotation",
        "gff_process",
        "cgc_finder",
        "substrate_prediction",
    ]
