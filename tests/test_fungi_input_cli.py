"""CLI checks for eukaryotic (fungal) input: --mode/--gff_type choices, required --gff_type in
protein mode, the gene/protein match check (before step 1 in easy_*), protein FASTA cleaning and
the eukaryote warning after pyrodigal."""
import logging
from pathlib import Path
from unittest.mock import patch

import pytest
from click.testing import CliRunner

from dbcan.configs.base_config import GeneralConfig
from dbcan.IO.fasta import ProteinProcessor, clean_protein_sequence, warn_if_eukaryotic_proteome
from dbcan.IO.gff_resolver import GffResolutionError
from dbcan.main import cli

FX = Path(__file__).parent / "_data" / "gff_dialects"


@pytest.fixture
def runner():
    return CliRunner()


def _easy(command, tmp_path, *extra):
    return [command, "--input_raw_data", str(FX / "funannotate2.faa"),
            "--output_dir", str(tmp_path / "out"), "--db_dir", str(tmp_path), *extra]


def test_mode_is_a_choice(runner, tmp_path):
    result = runner.invoke(cli, _easy("easy_CGC", tmp_path, "--mode", "fungi"))
    assert result.exit_code == 2
    assert "'fungi' is not one of" in result.output and "protein" in result.output


def test_gff_type_is_a_choice(runner, tmp_path):
    result = runner.invoke(cli, _easy("easy_CGC", tmp_path, "--mode", "protein",
                                      "--input_gff", str(FX / "funannotate2.gff3"), "--gff_type", "maker"))
    assert result.exit_code == 2
    assert "'maker' is not one of" in result.output and "funannotate" in result.output


@pytest.mark.parametrize("command", ["easy_CGC", "easy_substrate"])
def test_protein_mode_requires_gff_type(runner, tmp_path, command):
    with patch("dbcan.main._invoke_subset") as sub:
        result = runner.invoke(cli, _easy(command, tmp_path, "--mode", "protein",
                                          "--input_gff", str(FX / "funannotate2.gff3")))
    assert result.exit_code == 2
    assert "--gff_type is required" in result.output
    for t in ("prodigal", "NCBI_euk", "funannotate", "BRAKER", "GFF3", "GTF"):
        assert t in result.output
    sub.assert_not_called()


@pytest.mark.parametrize("command", ["easy_CGC", "easy_substrate"])
def test_wrong_gff_type_fails_before_step_1(runner, tmp_path, command):
    with patch("dbcan.main._invoke_subset") as sub:
        result = runner.invoke(cli, _easy(command, tmp_path, "--mode", "protein",
                                          "--input_gff", str(FX / "funannotate2.gff3"), "--gff_type", "NCBI_euk"))
    assert result.exit_code == 1
    assert "GFF check failed before step 1" in result.output
    assert "Suggested: --gff_type" in result.output
    sub.assert_not_called()
    assert (tmp_path / "out" / "gff_resolution.tsv").exists()


def test_right_gff_type_passes_check_case_insensitive(runner, tmp_path):
    with patch("dbcan.main._invoke_subset") as sub:
        result = runner.invoke(cli, _easy("easy_CGC", tmp_path, "--mode", "Protein",
                                          "--input_gff", str(FX / "funannotate2.gff3"), "--gff_type", "FUNANNOTATE"))
    assert result.exit_code == 0, result.output
    assert sub.call_count == 3
    assert sub.call_args_list[1].args[2]["gff_type"] == "funannotate"
    tsv = (tmp_path / "out" / "gff_resolution.tsv").read_text()
    assert tsv.count("\tmatched\t") == 5


def test_gff_min_match_zero_lets_a_bad_gff_through(runner, tmp_path):
    with patch("dbcan.main._invoke_subset") as sub:
        result = runner.invoke(cli, _easy("easy_CGC", tmp_path, "--mode", "protein",
                                          "--input_gff", str(FX / "funannotate2.gff3"), "--gff_type", "NCBI_euk",
                                          "--gff_min_match", "0"))
    assert result.exit_code == 0, result.output
    assert sub.call_count == 3


def test_gff_process_requires_gff_type_for_a_user_gff(runner, tmp_path):
    result = runner.invoke(cli, ["gff_process", "--output_dir", str(tmp_path), "--db_dir", str(tmp_path),
                                 "--input_gff", str(FX / "braker.gff3")])
    assert result.exit_code == 2
    assert "--gff_type is required" in result.output


def test_gff_process_fails_before_searches(runner, tmp_path):
    (tmp_path / "uniInput.faa").write_text((FX / "funannotate2.faa").read_text())
    with patch("dbcan.core.run_dbCAN_tcdb_diamond") as tc:
        result = runner.invoke(cli, ["gff_process", "--output_dir", str(tmp_path), "--db_dir", str(tmp_path),
                                     "--input_gff", str(FX / "funannotate2.gff3"), "--gff_type", "BRAKER",
                                     "--gff_min_match", "0.9"])
    # BRAKER key order still finds funannotate IDs (generic fallback): the check passes, the
    # (mocked) TC search runs and the next search fails on the empty database directory.
    tc.assert_called_once()
    assert result.exit_code != 0
    assert not isinstance(result.exception, GffResolutionError)
    assert "GFF check failed" not in result.output
    assert (tmp_path / "gff_resolution.tsv").read_text().count("\tmatched\t") == 5
    with patch("dbcan.core.run_dbCAN_tcdb_diamond") as tc:
        result = runner.invoke(cli, ["gff_process", "--output_dir", str(tmp_path), "--db_dir", str(tmp_path),
                                     "--input_gff", str(FX / "funannotate2.gff3"), "--gff_type", "JGI"])
    # a failed match check is a clean click error: exit 1, message, no traceback
    assert result.exit_code == 1
    assert isinstance(result.exception, SystemExit)
    assert "GFF check failed before the signature searches" in result.output
    assert "Suggested: --gff_type" in result.output
    assert "Traceback" not in result.output
    tc.assert_not_called()


# ------------------------------------------------------------------ protein FASTA cleaning
@pytest.mark.parametrize("raw,clean,trailing,internal", [
    ("MKV*", "MKV", True, 0),
    ("MKV.", "MKV", True, 0),
    ("MK.LV.", "MKXLV", True, 1),
    ("MK*L*V", "MKXLXV", False, 2),
    ("MKV", "MKV", False, 0),
])
def test_clean_protein_sequence(raw, clean, trailing, internal):
    assert clean_protein_sequence(raw) == (clean, trailing, internal)


def test_protein_mode_cleans_stops_for_diamond(tmp_path, caplog):
    faa = tmp_path / "in.faa"
    faa.write_text(">a gffread\nMKV*\n>b\nMK.LV.\n>c\nMKLV\n")
    cfg = GeneralConfig(input_raw_data=str(faa), output_dir=str(tmp_path / "out"), mode="protein")
    with caplog.at_level(logging.INFO, logger="dbcan.IO.fasta"):
        assert ProteinProcessor(cfg).process_input()
    seqs = (tmp_path / "out" / "uniInput.faa").read_text().split()
    assert seqs == [">a", "MKV", ">b", "MKXLV", ">c", "MKLV"]
    assert "removed a trailing stop ('*' or '.') from 2 sequences" in caplog.text
    assert "with 'X' in 1 sequences" in caplog.text


# ------------------------------------------------------------------ eukaryote warning
def test_warns_on_fungal_like_pyrodigal_output(caplog):
    lengths = [60] * 120 + [140] * 60 + [400] * 120     # median 140 aa, 40% < 100 aa
    with caplog.at_level(logging.WARNING, logger="dbcan.IO.fasta"):
        assert warn_if_eukaryotic_proteome(lengths)
    assert "--mode protein" in caplog.text and "funannotate" in caplog.text


def test_no_warning_for_bacterial_lengths_or_tiny_inputs():
    assert not warn_if_eukaryotic_proteome([280] * 270 + [80] * 30)   # E. coli-like
    assert not warn_if_eukaryotic_proteome([60] * 50)                 # too few genes to judge
