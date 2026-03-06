import os
import pytest
from click.testing import CliRunner
from pathlib import Path

from dbcan.main import cli

# Test root directory and data directory
TEST_ROOT = Path(__file__).parent
DATA_ROOT = TEST_ROOT / "_data"

# Test data file paths
TEST_PROTEIN = DATA_ROOT / "EscheriaColiK12MG1655.faa"
TEST_NUCLEOTIDE = DATA_ROOT / "EscheriaColiK12MG1655.fna"
TEST_GFF = DATA_ROOT / "EscheriaColiK12MG1655.gff"


@pytest.fixture
def runner():
    return CliRunner()


@pytest.fixture
def actual_db_dir(tmp_path, runner):
    db_dir = tmp_path / "db"
    db_dir.mkdir()
    result = runner.invoke(cli, [
        'database',
        '--db_dir', str(db_dir),
        '--aws_s3'
    ])
    if result.exit_code != 0:
        print(f"Database command failed with exit code {result.exit_code}")
        print(f"Output: {result.output}")
        print(f"Exception: {result.exception}")
        pytest.skip("Failed to create database, skipping test")
    return str(db_dir)


class TestEasySubstrate:
    # NOTE: All patches removed to run the real pipeline.
    # If the environment lacks required external tools, mark or skip accordingly.
    @pytest.mark.integration
    def test_easy_substrate_protein(self, runner, actual_db_dir, tmp_path):
        """
        Integration test for easy_substrate command with protein input (no mocks).
        """
        output_dir = tmp_path / "output_protein"
        output_dir.mkdir()
        output_dir_str = str(output_dir)

        assert TEST_PROTEIN.exists(), f"Test protein file not found at {TEST_PROTEIN}"
        assert TEST_GFF.exists(), f"Test GFF file not found at {TEST_GFF}"

        print("Running test with:")
        print(f"  TEST_PROTEIN: {TEST_PROTEIN}")
        print(f"  TEST_GFF: {TEST_GFF}")
        print(f"  db_dir: {actual_db_dir}")
        print(f"  output_dir: {output_dir_str}")

        result = runner.invoke(cli, [
            'easy_substrate',
            '--mode', 'protein',
            '--input_raw_data', str(TEST_PROTEIN),
            '--input_gff', str(TEST_GFF),
            '--gff_type', 'NCBI_prok',
            '--output_dir', output_dir_str,
            '--db_dir', actual_db_dir,
            #'--threads', '2'
        ])

        if result.exit_code != 0:
            print("STDOUT/STDERR:")
            print(result.output)
            if result.exception:
                print("Exception:", result.exception)

        assert result.exit_code == 0, f"Command failed: {result.output}"

        # Key outputs
        for fn in ("overview.tsv", "cgc_standard_out.tsv", "substrate_prediction.tsv"):
            assert (output_dir / fn).exists(), f"{fn} not found"

    @pytest.mark.integration
    def test_easy_substrate_nucleotide(self, runner, actual_db_dir, tmp_path):
        """
        Integration test for easy_substrate command with nucleotide input (no mocks).
        """
        output_dir = tmp_path / "output_nucleotide"
        output_dir.mkdir()
        output_dir_str = str(output_dir)

        assert TEST_NUCLEOTIDE.exists(), f"Test nucleotide file not found at {TEST_NUCLEOTIDE}"

        print("Running test with:")
        print(f"  TEST_NUCLEOTIDE: {TEST_NUCLEOTIDE}")
        print(f"  db_dir: {actual_db_dir}")
        print(f"  output_dir: {output_dir_str}")

        # Provide path where GFF will be produced (pipeline should generate)
        intended_gff = output_dir / "uniInput.gff"

        result = runner.invoke(cli, [
            'easy_substrate',
            '--mode', 'prok',
            '--input_raw_data', str(TEST_NUCLEOTIDE),
            '--gff_type', 'prodigal',
            '--output_dir', output_dir_str,
            '--db_dir', actual_db_dir,
          #  '--threads', '2'
        ])

        if result.exit_code != 0:
            print("STDOUT/STDERR:")
            print(result.output)
            if result.exception:
                print("Exception:", result.exception)

        assert result.exit_code == 0, f"Command failed: {result.output}"

        for fn in ("overview.tsv", "cgc_standard_out.tsv", "substrate_prediction.tsv"):
            assert (output_dir / fn).exists(), f"{fn} not found"


class TestStructureMethod:
    """Tests for structure search method and overview integration."""

    def test_methods_accept_structure(self, runner):
        """CLI accepts --methods with structure (no 'Invalid method' error)."""
        result = runner.invoke(cli, [
            'CAZyme_annotation',
            '--methods', 'diamond,structure',
            '--input_raw_data', str(TEST_PROTEIN),
            '--output_dir', '/tmp/out',
            '--db_dir', '/tmp/db',
            '--mode', 'protein',
        ])
        assert 'Invalid method' not in (result.output or ''), "structure should be a valid method"

    def test_overview_constants_include_structure(self):
        """Overview constants define Structure column and field."""
        import dbcan.constants.OverviewGenerator_constants as O
        assert 'Structure' in O.OVERVIEW_COLUMNS
        assert O.STRUCTURE_FIELD == 'Structure'
        assert O.STRUCTURE_RESULT_FILE == 'structure_search_results.tsv'

    @pytest.mark.integration
    def test_cazyme_annotation_without_structure_creates_empty_structure_file(self, runner, actual_db_dir, tmp_path):
        """With --methods diamond,hmm,dbCANsub (no structure), empty structure result file is created and overview has Structure column."""
        output_dir = tmp_path / "output_no_structure"
        output_dir.mkdir()
        assert TEST_PROTEIN.exists()
        result = runner.invoke(cli, [
            'CAZyme_annotation',
            '--mode', 'protein',
            '--input_raw_data', str(TEST_PROTEIN),
            '--output_dir', str(output_dir),
            '--db_dir', actual_db_dir,
            '--methods', 'diamond,hmm,dbCANsub',
        ])
        assert result.exit_code == 0, f"Command failed: {result.output}"
        structure_file = output_dir / "structure_search_results.tsv"
        assert structure_file.exists(), "structure_search_results.tsv should exist (empty when structure not in methods)"
        overview = output_dir / "overview.tsv"
        assert overview.exists()
        with open(overview) as f:
            header = f.readline()
        assert "Structure" in header, "overview.tsv should contain Structure column"


