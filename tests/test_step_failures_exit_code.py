"""Failed pipeline steps must make the command exit with a non-zero status.

Before, these failures were only logged and the commands exited with 0, so workflow
managers treated empty or stale outputs as a successful run.
"""
import pytest
from click.testing import CliRunner

from dbcan.configs.base_config import OverviewGeneratorConfig
from dbcan.IO.OverviewGenerator import OverviewGenerator
from dbcan.main import cli


@pytest.fixture
def runner():
    return CliRunner()


def test_cazyme_annotation_fails_on_fasta_without_sequences(runner, tmp_path):
    fasta = tmp_path / "no_records.faa"
    fasta.write_text("this is not a FASTA record\n")
    db = tmp_path / "db"
    db.mkdir()
    result = runner.invoke(cli, [
        "CAZyme_annotation", "--mode", "protein", "--input_raw_data", str(fasta),
        "--output_dir", str(tmp_path / "out"), "--db_dir", str(db), "--methods", "hmm",
    ])
    assert result.exit_code != 0
    assert "Input processing failed" in str(result.exception)


def test_cazyme_annotation_fails_when_gene_prediction_fails(runner, tmp_path):
    # pyrodigal needs >= 20 kb to train in prok mode
    fna = tmp_path / "short.fna"
    fna.write_text(">c1\n" + "ATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGA" * 10 + "\n")
    db = tmp_path / "db"
    db.mkdir()
    result = runner.invoke(cli, [
        "CAZyme_annotation", "--mode", "prok", "--input_raw_data", str(fna),
        "--output_dir", str(tmp_path / "out"), "--db_dir", str(db), "--methods", "hmm",
    ])
    assert result.exit_code != 0
    assert "Input processing failed" in str(result.exception)


def test_gff_process_requires_existing_gff_before_running_searches(runner, tmp_path):
    db = tmp_path / "db"
    db.mkdir()
    result = runner.invoke(cli, ["gff_process", "--output_dir", str(tmp_path), "--db_dir", str(db)])
    assert result.exit_code == 2  # click.UsageError
    assert "GFF file not found" in result.output


@pytest.mark.parametrize("command", ["easy_CGC", "easy_substrate"])
def test_easy_commands_require_gff_for_protein_input(runner, tmp_path, command):
    faa = tmp_path / "p.faa"
    faa.write_text(">p1\nMKV\n")
    result = runner.invoke(cli, [
        command, "--mode", "protein", "--input_raw_data", str(faa),
        "--output_dir", str(tmp_path / "out"), "--db_dir", str(tmp_path),
    ])
    assert result.exit_code == 2
    assert "--input_gff" in result.output


def test_cgc_finder_fails_without_cgc_gff(runner, tmp_path):
    result = runner.invoke(cli, ["cgc_finder", "--output_dir", str(tmp_path)])
    assert result.exit_code != 0
    assert "CGC identification failed" in str(result.exception)


def test_overview_generation_error_is_raised(tmp_path, monkeypatch):
    (tmp_path / "dbCAN_hmm_results.tsv").write_text(
        "HMM Name\tHMM Length\tTarget Name\tTarget Length\ti-Evalue\tHMM From\tHMM To\t"
        "Target From\tTarget To\tCoverage\tHMM File Name\n"
        "GH5.hmm\t300\tp1\t400\t1e-50\t1\t290\t10\t300\t0.96\tdbCAN\n"
    )
    gen = OverviewGenerator(OverviewGeneratorConfig(output_dir=str(tmp_path)))

    def boom(*a, **k):
        raise ValueError("boom")

    monkeypatch.setattr(gen, "aggregate_data", boom)
    with pytest.raises(ValueError, match="boom"):
        gen.run()
