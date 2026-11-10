"""overview.tsv: CAZy family parsing of DIAMOND subject IDs, and merging only the selected --methods."""
import logging
from pathlib import Path

import pandas as pd
import pytest

import dbcan.constants.diamond_constants as D
from dbcan.configs.base_config import OverviewGeneratorConfig
from dbcan.IO.OverviewGenerator import OverviewGenerator
from dbcan.IO.table import read_tsv


@pytest.mark.parametrize(
    "sseqid, expected",
    [
        ("CAD7272370.1|GH33", "GH33"),
        # GenBank accessions that start with a class prefix are not families
        ("AAK49401.1|GH19", "GH19"),
        ("CEA12345.1|CE4", "CE4"),
        ("GTB00001.1|GT2", "GT2"),
        # new CAZyDB headers repeat the subfamily and carry an EC number
        ("1803517A|GH13_6|GH13_6|3.2.1.1", "GH13_6"),
        ("XYZ.1|GH5_2|CBM1|GH5_2", "GH5_2+CBM1"),
        # JGI headers: families, numeric protein id, source file
        ("CE16|390012|Krezon1_GeneCatalog_proteins_20220127.aa.fasta", "CE16"),
        ("GH15|CBM20|716895|Lacbi81306_1_GeneCatalog_proteins_20140617.aa.fasta", "GH15+CBM20"),
        ("CBM13|CBM13|CBM13|PL1|1000032|Neolan1_GeneCatalog_proteins_20200610.aa.fasta", "CBM13+PL1"),
        # nothing that looks like a family: unchanged
        ("Q12345|unknown", "Q12345|unknown"),
        ("", ""),
    ],
)
def test_extract_cazy_id(sseqid, expected):
    assert OverviewGenerator.extract_cazy_id(sseqid) == expected


HMM_COLS = ["HMM Name", "HMM Length", "Target Name", "Target Length", "i-Evalue",
            "HMM From", "HMM To", "Target From", "Target To", "Coverage", "HMM File Name"]
SUB_COLS = ["Subfam Name", "Subfam Composition", "Subfam EC", "Substrate", "HMM Length", "Target Name",
            "Target Length", "i-Evalue", "HMM From", "HMM To", "Target From", "Target To", "Coverage",
            "HMM File Name"]


def _write_all_results(out: Path):
    """Result files of all three methods, as left behind by an earlier full run."""
    pd.DataFrame([["p1", "AAK49401.1|GH5_2", 50, 300, 1, 0, 1, 300, 1, 300, 1e-120, 500]],
                 columns=D.CAZY_COLUMN_NAMES).to_csv(out / "diamond.out", sep="\t", index=False)
    pd.DataFrame([["GH5_2", 300, "p1", 400, 1e-50, 1, 290, 10, 300, 0.96, "dbCAN"]],
                 columns=HMM_COLS).to_csv(out / "dbCAN_hmm_results.tsv", sep="\t", index=False)
    pd.DataFrame([["GH5_e12", "GH5:10", "3.2.1.4:3", "cellulose", 300, "p1", 400, 1e-50, 1, 290, 10, 300,
                   0.96, "dbCAN-sub"]],
                 columns=SUB_COLS).to_csv(out / "dbCANsub_hmm_results.tsv", sep="\t", index=False)


def test_overview_all_methods_by_default(tmp_path):
    _write_all_results(tmp_path)
    OverviewGenerator(OverviewGeneratorConfig(output_dir=str(tmp_path))).run()
    row = read_tsv(tmp_path / "overview.tsv").set_index("Gene ID").loc["p1"]
    assert row["#ofTools"] == 3
    assert row["DIAMOND"] == "GH5_2"


def test_overview_ignores_stale_files_of_unselected_methods(tmp_path, caplog):
    _write_all_results(tmp_path)
    caplog.set_level(logging.INFO, logger="dbcan.IO.OverviewGenerator")
    OverviewGenerator(OverviewGeneratorConfig(output_dir=str(tmp_path), methods=("hmm",))).run()
    row = read_tsv(tmp_path / "overview.tsv").set_index("Gene ID").loc["p1"]
    assert row["#ofTools"] == 1
    assert row["DIAMOND"] == "-"
    assert row["dbCAN_sub"] == "-"
    assert row["dbCAN_hmm"] == "GH5_2(10-300)"
    assert "Ignoring diamond.out" in caplog.text
    assert "Ignoring dbCANsub_hmm_results.tsv" in caplog.text


def test_cazyme_annotation_passes_methods_to_overview(tmp_path):
    """run_dbCAN_CAZyme_annotation hands its methods to the overview step."""
    from unittest.mock import patch
    from types import SimpleNamespace
    import dbcan.core as core

    _write_all_results(tmp_path)
    cfg = OverviewGeneratorConfig(output_dir=str(tmp_path))
    hmm_cfg = SimpleNamespace(output_dir=str(tmp_path), output_file="dbCAN_hmm_results.tsv")
    hmm_table = (tmp_path / "dbCAN_hmm_results.tsv").read_text()

    def fake_hmmer(_cfg):   # the selected method rewrites its result file
        (tmp_path / "dbCAN_hmm_results.tsv").write_text(hmm_table)

    with patch.object(core, "run_dbCAN_hmmer", side_effect=fake_hmmer):
        core.run_dbCAN_CAZyme_annotation(
            SimpleNamespace(output_dir=str(tmp_path), output_file="diamond.out"),
            hmm_cfg,
            SimpleNamespace(output_dir=str(tmp_path)),
            cfg,
            ("hmm",),
        )
    assert cfg.methods == ("hmm",)
    assert read_tsv(tmp_path / "overview.tsv").set_index("Gene ID").loc["p1", "#ofTools"] == 1


def test_unknown_method_is_an_error(tmp_path):
    with pytest.raises(ValueError):
        OverviewGenerator(OverviewGeneratorConfig(output_dir=str(tmp_path), methods=("blast",))).run()


def test_failed_selected_method_does_not_merge_its_stale_file(tmp_path, caplog):
    """DIAMOND is selected but fails now: its diamond.out from an earlier run must not be merged."""
    from unittest.mock import patch
    from types import SimpleNamespace
    import dbcan.core as core

    _write_all_results(tmp_path)
    hmm_table = (tmp_path / "dbCAN_hmm_results.tsv").read_text()
    cfg = OverviewGeneratorConfig(output_dir=str(tmp_path))
    with patch.object(core, "run_dbCAN_cazy_diamond", side_effect=RuntimeError("diamond crashed")), \
         patch.object(core, "run_dbCAN_hmmer",
                      side_effect=lambda c: (tmp_path / "dbCAN_hmm_results.tsv").write_text(hmm_table)):
        with pytest.raises(RuntimeError, match="DIAMOND CAZy: diamond crashed"):
            core.run_dbCAN_CAZyme_annotation(
                SimpleNamespace(output_dir=str(tmp_path), output_file="diamond.out", column_names=None),
                SimpleNamespace(output_dir=str(tmp_path), output_file="dbCAN_hmm_results.tsv"),
                SimpleNamespace(output_dir=str(tmp_path)),
                cfg,
                ("diamond", "hmm"),
            )
    row = read_tsv(tmp_path / "overview.tsv").set_index("Gene ID").loc["p1"]
    assert row["DIAMOND"] == "-"
    assert row["#ofTools"] == 1
    # the earlier run's diamond.out was replaced by an empty table
    assert len(read_tsv(tmp_path / "diamond.out")) == 0


def test_ignoring_unselected_method_is_a_warning(tmp_path, caplog):
    _write_all_results(tmp_path)
    caplog.set_level(logging.WARNING, logger="dbcan.IO.OverviewGenerator")
    OverviewGenerator(OverviewGeneratorConfig(output_dir=str(tmp_path), methods=("hmm",))).run()
    assert any(r.levelno == logging.WARNING and "Ignoring diamond.out" in r.getMessage() for r in caplog.records)


def test_empty_placeholder_of_unselected_method_is_not_warned(tmp_path, caplog):
    """core writes a header-only diamond.out when DIAMOND is not selected; that is no stale result."""
    _write_all_results(tmp_path)
    pd.DataFrame(columns=D.CAZY_COLUMN_NAMES).to_csv(tmp_path / "diamond.out", sep="\t", index=False)
    caplog.set_level(logging.WARNING, logger="dbcan.IO.OverviewGenerator")
    OverviewGenerator(OverviewGeneratorConfig(output_dir=str(tmp_path), methods=("hmm", "dbCANsub"))).run()
    assert "Ignoring diamond.out" not in caplog.text
    assert read_tsv(tmp_path / "overview.tsv").set_index("Gene ID").loc["p1", "DIAMOND"] == "-"


def test_unreadable_result_of_selected_method_is_an_error(tmp_path):
    """A selected method's table that pandas cannot parse must fail the overview, not be skipped."""
    _write_all_results(tmp_path)
    (tmp_path / "diamond.out").write_text("Gene ID\tCAZy ID\np1\tGH5\np2\tGH5\textra\tfields\n")
    with pytest.raises(RuntimeError, match="Error loading diamond results"):
        OverviewGenerator(OverviewGeneratorConfig(output_dir=str(tmp_path))).run()


def test_selected_result_missing_columns_is_an_error(tmp_path):
    _write_all_results(tmp_path)
    pd.DataFrame([["p1", "GH5_2"]], columns=["Target Name", "HMM Name"]).to_csv(
        tmp_path / "dbCAN_hmm_results.tsv", sep="\t", index=False)
    with pytest.raises(RuntimeError, match="Missing columns"):
        OverviewGenerator(OverviewGeneratorConfig(output_dir=str(tmp_path), methods=("hmm", "diamond"))).run()


@pytest.mark.parametrize("content", ["", "Target Name\tHMM Name\n"], ids=["zero-byte", "header-only-other-columns"])
def test_empty_result_of_selected_method_is_not_an_error(tmp_path, content):
    """A selected method that found nothing leaves an empty table: no hits are lost, so no error."""
    _write_all_results(tmp_path)
    (tmp_path / "dbCAN_hmm_results.tsv").write_text(content)
    OverviewGenerator(OverviewGeneratorConfig(output_dir=str(tmp_path))).run()
    row = read_tsv(tmp_path / "overview.tsv").set_index("Gene ID").loc["p1"]
    assert row["dbCAN_hmm"] == "-"
    assert row["#ofTools"] == 2


def test_header_only_result_of_selected_method_is_not_an_error(tmp_path):
    _write_all_results(tmp_path)
    pd.DataFrame(columns=D.CAZY_COLUMN_NAMES).to_csv(tmp_path / "diamond.out", sep="\t", index=False)
    OverviewGenerator(OverviewGeneratorConfig(output_dir=str(tmp_path), methods=("diamond", "hmm"))).run()
    row = read_tsv(tmp_path / "overview.tsv").set_index("Gene ID").loc["p1"]
    assert row["DIAMOND"] == "-"
    assert row["#ofTools"] == 1


def test_unreadable_result_of_unselected_method_is_ignored(tmp_path):
    _write_all_results(tmp_path)
    (tmp_path / "diamond.out").write_text("Gene ID\tCAZy ID\np1\tGH5\np2\tGH5\textra\tfields\n")
    OverviewGenerator(OverviewGeneratorConfig(output_dir=str(tmp_path), methods=("hmm",))).run()
    assert read_tsv(tmp_path / "overview.tsv").set_index("Gene ID").loc["p1", "#ofTools"] == 1


def test_cazyme_annotation_fails_on_unreadable_selected_result(tmp_path):
    """core records the overview step as failed, so the command exits non-zero."""
    from unittest.mock import patch
    from types import SimpleNamespace
    import dbcan.core as core

    _write_all_results(tmp_path)
    cfg = OverviewGeneratorConfig(output_dir=str(tmp_path))

    def corrupt_hmmer(_cfg):
        (tmp_path / "dbCAN_hmm_results.tsv").write_text("Target Name\tHMM Name\np1\tGH5\np2\tGH5\tx\ty\n")

    with patch.object(core, "run_dbCAN_hmmer", side_effect=corrupt_hmmer):
        with pytest.raises(RuntimeError, match="CAZyme overview: Error loading dbcan_hmm results"):
            core.run_dbCAN_CAZyme_annotation(
                SimpleNamespace(output_dir=str(tmp_path), output_file="diamond.out", column_names=None),
                SimpleNamespace(output_dir=str(tmp_path), output_file="dbCAN_hmm_results.tsv"),
                SimpleNamespace(output_dir=str(tmp_path)),
                cfg,
                ("hmm",),
            )
