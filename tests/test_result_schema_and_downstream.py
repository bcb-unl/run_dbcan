"""Public-schema and downstream-regression tests for domain resolution.

The public result tables keep their historical columns (also when empty), and
the domains accepted by pyHMMER resolution survive overview, CGC and Pfam
post-processing.
"""

from types import SimpleNamespace

import pandas as pd
import pytest

import dbcan.constants.base_constants as B
import dbcan.constants.process_dbcan_sub_constants as PS
import dbcan.constants.process_utils_constants as PU
from dbcan.configs.base_config import OverviewGeneratorConfig
from dbcan.IO.OverviewGenerator import OverviewGenerator
from dbcan.process.process_utils import (
    process_cgc_null_pfam_annotation,
    process_cgc_sig_results,
    process_results,
)


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def hmm_row(hmm_name, target, ali_from, ali_to, i_evalue, hmm_len=400, target_len=600):
    return [
        hmm_name,
        hmm_len,
        target,
        target_len,
        i_evalue,
        1,
        hmm_len,
        ali_from,
        ali_to,
        0.9,
        "dbCAN",
    ]


def _architecture_rows():
    """GH18 (1-400) and a CBM18 (300-334) inside its alignment span."""
    return [
        hmm_row("GH18.hmm", "prot1", 1, 400, 1e-120),
        hmm_row("CBM18.hmm", "prot1", 300, 334, 1e-30, hmm_len=40),
    ]


def _write_tsv(path, columns, rows):
    pd.DataFrame(rows, columns=columns).to_csv(path, sep="\t", index=False)


# ---------------------------------------------------------------------------
# Output schema
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "result_file",
    [
        B.DBCAN_HMM_RESULT_FILE,
        B.DBCAN_SUB_HMM_RAW_FILE,
        B.TF_HMM_RESULT_FILE,
        B.STP_HMM_RESULT_FILE,
        B.PFAM_HMM_RESULT_FILE,
    ],
)
@pytest.mark.parametrize("empty", [True, False])
def test_public_columns_are_identical_with_and_without_hits(tmp_path, result_file, empty):
    out = tmp_path / result_file
    rows = [] if empty else _architecture_rows()
    process_results(rows, out)

    written = pd.read_csv(out, sep="\t")
    assert list(written.columns) == PU.HMMER_COLUMN_NAMES
    # the CBM18 inside the GH18 alignment span is removed (alignment-span overlap)
    assert len(written) == (0 if empty else 1)
    # only the result table is written (no side files)
    assert [p.name for p in tmp_path.iterdir()] == [result_file]


# ---------------------------------------------------------------------------
# Phase 11 Test J: OverviewGenerator regression
# ---------------------------------------------------------------------------
def _build_overview_dir(tmp_path):
    _write_tsv(
        tmp_path / B.DBCAN_HMM_RESULT_FILE,
        PU.HMMER_COLUMN_NAMES,
        [
            hmm_row("GH18.hmm", "prot1", 1, 400, 1e-120),
            hmm_row("CBM18.hmm", "prot1", 300, 334, 1e-30, hmm_len=40),
        ],
    )
    _write_tsv(
        tmp_path / B.DBCAN_SUB_HMM_RESULT_FILE,
        PS.DBCAN_SUB_COLUMN_NAMES,
        [["GH18_e1", "-", "-", "-", 400, "prot1", 600, 1e-118, 1, 400, 1, 400, 0.9, "dbCAN_sub"]],
    )
    return OverviewGenerator(OverviewGeneratorConfig(output_dir=str(tmp_path)))


def test_overview_groups_hmm_domains_under_a_spanning_subfamily(tmp_path):
    """One dbCAN-sub subfamily over both dbCAN HMM domains: one recommended name (span overlap)."""
    _build_overview_dir(tmp_path).run()
    overview = pd.read_csv(tmp_path / B.OVERVIEW_FILE, sep="\t")
    row = overview.loc[overview["Gene ID"] == "prot1"].iloc[0]
    assert str(row["Recommend Results"]).split("|") == ["GH18_e1"]
    # the per-tool column still reports both HMM domains
    assert "CBM18(300-334)" in str(row["dbCAN_hmm"])


def test_overview_no_longer_regroups_already_resolved_hmm_hits(tmp_path):
    """Two dbCAN-HMM hits are never merged with each other by the overview.

    They only reach ``dbCAN_hmm_results.tsv`` after domain resolution accepted
    them, so re-running an overlap rule over them would undo that decision.
    """
    generator = _build_overview_dir(tmp_path)
    # Spans overlap heavily (CBM18 lies inside GH18's span); the old HMM-vs-HMM
    # grouping loop collapsed these two to a single representative.
    annotations = [("GH18", 1, 400, "hmm"), ("CBM18", 300, 334, "hmm")]
    grouped = generator.graph_based_grouping(annotations)
    assert {g[0] for g in grouped} == {"GH18", "CBM18"}


# ---------------------------------------------------------------------------
# CGC downstream regression
# ---------------------------------------------------------------------------
def test_k_cgc_post_processing_keeps_multidomain_and_cross_class_annotations(tmp_path):
    # STP comes from pyhmmer and is already resolved: whatever STP_hmm_results.tsv
    # holds (here two overlapping domains) is passed through.
    _write_tsv(
        tmp_path / B.STP_HMM_RESULT_FILE,
        PU.HMMER_COLUMN_NAMES,
        [
            hmm_row("STP_big.hmm", "prot1", 1, 400, 1e-120),
            hmm_row("STP_small.hmm", "prot1", 300, 334, 1e-30, hmm_len=40),
        ],
    )
    # TC comes from DIAMOND and is not pre-resolved: two competing TC hits over
    # the same region, one of which must win.
    _write_tsv(
        tmp_path / B.TCDB_DIAMOND_OUTPUT,
        PU.CGC_SIG_RESULT_COLUMN,
        [
            ["TC_a", 400, "prot1", 600, 1e-50, 1, 400, 1, 400, 0.9, "tcdb"],
            ["TC_b", 400, "prot1", 600, 1e-20, 1, 400, 50, 380, 0.9, "tcdb"],
        ],
    )

    cfg = SimpleNamespace(output_dir=str(tmp_path))
    process_cgc_sig_results(cfg, None, SimpleNamespace(fungi=False), cfg, cfg, cfg)

    out = pd.read_csv(tmp_path / B.GFF_CGC_SIG_FILE, sep="\t")
    observed = set(zip(out["Type"], out["Target From"], out["Target To"]))
    # both STP domains survive (no second overlap filter of a resolved source)
    assert ("STP", 1, 400) in observed
    assert ("STP", 300, 334) in observed
    # TC and STP cover the same region but are different signature classes
    assert ("TC", 1, 400) in observed
    # the two competing TC hits are still resolved within their own source
    assert ("TC", 50, 380) not in observed
    assert len(out) == 3


def test_cgc_post_processing_handles_no_inputs(tmp_path):
    cfg = SimpleNamespace(output_dir=str(tmp_path))
    process_cgc_sig_results(cfg, None, SimpleNamespace(fungi=False), cfg, cfg, cfg)
    out = pd.read_csv(tmp_path / B.GFF_CGC_SIG_FILE, sep="\t")
    # 'Source' (TC / TF_prok / TF_fungi / STP / ...) keeps TF DIAMOND and TF HMM apart for the CGC step
    assert list(out.columns) == PU.CGC_SIG_RESULT_COLUMN + ["Type", PU.SIGNATURE_SOURCE_COLUMN,
                                                          PU.SIGNATURE_BITSCORE_COLUMN]
    assert out.empty


# ---------------------------------------------------------------------------
# Pfam downstream regression
# ---------------------------------------------------------------------------
def test_l_pfam_post_processing_does_not_refilter_resolved_domains(tmp_path):
    pfam_path = tmp_path / B.PFAM_HMM_RESULT_FILE
    _write_tsv(
        pfam_path,
        PU.HMMER_COLUMN_NAMES,
        [
            hmm_row("PF00704.hmm", "prot1", 1, 400, 1e-120),
            hmm_row("PF00187.hmm", "prot1", 300, 334, 1e-30, hmm_len=40),
        ],
    )

    process_cgc_null_pfam_annotation(SimpleNamespace(output_dir=str(tmp_path)))

    out = pd.read_csv(pfam_path, sep="\t")
    assert list(out.columns) == PU.HMMER_COLUMN_NAMES
    assert set(out["HMM Name"]) == {"PF00704.hmm", "PF00187.hmm"}


def test_pfam_post_processing_creates_an_empty_table_when_input_is_missing(tmp_path):
    process_cgc_null_pfam_annotation(SimpleNamespace(output_dir=str(tmp_path)))
    out = pd.read_csv(tmp_path / B.PFAM_HMM_RESULT_FILE, sep="\t")
    assert list(out.columns) == PU.HMMER_COLUMN_NAMES
    assert out.empty
