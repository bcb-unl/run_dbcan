"""Public-schema and downstream-regression tests for domain resolution.

Phase 9 (public schemas), Phase 11 Tests I/J/K/L: the architecture accepted by
PyHMMER resolution must survive overview, CGC and Pfam post-processing, and the
public result tables must keep their historical columns in every mode.
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
    alignment_details_path,
    process_cgc_null_pfam_annotation,
    process_cgc_sig_results,
    process_results,
)


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def hmm_row(hmm_name, target, ali_from, ali_to, i_evalue, matched=None, hmm_len=400, target_len=600):
    row = [
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
    if matched is not None:
        row.append(matched)
    return row


#: The GH18/CBM18 architecture: the CBM sits inside the catalytic domain's
#: insert-state region, so the spans overlap but the match states do not.
GH18_MATCHED = "1:250,340:400"
CBM18_MATCHED = "300:334"


def _architecture_rows(with_matched=True):
    if with_matched:
        return [
            hmm_row("GH18.hmm", "prot1", 1, 400, 1e-120, GH18_MATCHED),
            hmm_row("CBM18.hmm", "prot1", 300, 334, 1e-30, CBM18_MATCHED, hmm_len=40),
        ]
    return [
        hmm_row("GH18.hmm", "prot1", 1, 400, 1e-120),
        hmm_row("CBM18.hmm", "prot1", 300, 334, 1e-30, hmm_len=40),
    ]


def _write_tsv(path, columns, rows):
    pd.DataFrame(rows, columns=columns).to_csv(path, sep="\t", index=False)


# ---------------------------------------------------------------------------
# Phase 11 Test I: output schema consistency
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
@pytest.mark.parametrize("matched_overlap", [True, False])
@pytest.mark.parametrize("empty", [True, False])
def test_i_public_columns_are_identical_across_modes_and_emptiness(
    tmp_path, result_file, matched_overlap, empty
):
    out = tmp_path / result_file
    rows = [] if empty else _architecture_rows(with_matched=matched_overlap)
    process_results(rows, out, matched_overlap=matched_overlap)

    written = pd.read_csv(out, sep="\t")
    assert list(written.columns) == PU.HMMER_COLUMN_NAMES
    assert PU.MATCHED_INTERVALS_COLUMN not in written.columns
    assert len(written) == (0 if empty else 2 if matched_overlap else 1)


def test_alignment_details_sidecar_only_exists_in_matched_mode(tmp_path):
    out = tmp_path / B.DBCAN_HMM_RESULT_FILE
    details = alignment_details_path(out)
    assert details.name == "dbCAN_hmm_alignment_details.tsv"

    process_results(_architecture_rows(with_matched=False), out, matched_overlap=False)
    assert not details.exists()

    process_results(_architecture_rows(), out, matched_overlap=True)
    assert details.exists()
    df = pd.read_csv(details, sep="\t")
    assert list(df.columns) == PU.ALIGNMENT_DETAILS_COLUMNS
    assert dict(zip(df["HMM Name"], df[PU.MATCHED_INTERVALS_COLUMN])) == {
        "GH18.hmm": GH18_MATCHED,
        "CBM18.hmm": CBM18_MATCHED,
    }


def test_empty_result_still_writes_an_empty_sidecar(tmp_path):
    out = tmp_path / B.DBCAN_HMM_RESULT_FILE
    process_results([], out, matched_overlap=True)
    details = pd.read_csv(alignment_details_path(out), sep="\t")
    assert list(details.columns) == PU.ALIGNMENT_DETAILS_COLUMNS
    assert details.empty


def test_sidecar_naming_for_the_dbcan_sub_raw_table(tmp_path):
    out = tmp_path / B.DBCAN_SUB_HMM_RAW_FILE
    assert alignment_details_path(out).name == "dbCANsub_hmm_alignment_details.tsv"


# ---------------------------------------------------------------------------
# Phase 11 Test J: OverviewGenerator regression
# ---------------------------------------------------------------------------
def _build_overview_dir(tmp_path, with_sidecars=True):
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
    if with_sidecars:
        _write_tsv(
            tmp_path / "dbCAN_hmm_alignment_details.tsv",
            PU.ALIGNMENT_DETAILS_COLUMNS,
            [
                ["GH18.hmm", "prot1", 600, 1, 400, GH18_MATCHED, 1e-120, 0.9, "dbCAN"],
                ["CBM18.hmm", "prot1", 600, 300, 334, CBM18_MATCHED, 1e-30, 0.9, "dbCAN"],
            ],
        )
        _write_tsv(
            tmp_path / "dbCANsub_hmm_alignment_details.tsv",
            PU.ALIGNMENT_DETAILS_COLUMNS,
            [["GH18_e1.hmm", "prot1", 600, 1, 400, GH18_MATCHED, 1e-118, 0.9, "dbCAN_sub"]],
        )
    return OverviewGenerator(OverviewGeneratorConfig(output_dir=str(tmp_path)))


def test_j_overview_preserves_catalytic_plus_accessory_architecture(tmp_path):
    _build_overview_dir(tmp_path).run()
    overview = pd.read_csv(tmp_path / B.OVERVIEW_FILE, sep="\t")

    row = overview.loc[overview["Gene ID"] == "prot1"].iloc[0]
    recommended = str(row["Recommend Results"]).split("|")
    assert "CBM18" in recommended, "accessory domain lost from Recommend Results"
    assert recommended == ["GH18_e1", "CBM18"]
    # the per-tool columns must still report both HMM domains
    assert "CBM18(300-334)" in str(row["dbCAN_hmm"])


def test_j_without_sidecars_the_overview_falls_back_to_span_overlap(tmp_path):
    """Legacy runs write no sidecar; the historical span behaviour is retained."""
    _build_overview_dir(tmp_path, with_sidecars=False).run()
    overview = pd.read_csv(tmp_path / B.OVERVIEW_FILE, sep="\t")
    row = overview.loc[overview["Gene ID"] == "prot1"].iloc[0]
    assert str(row["Recommend Results"]).split("|") == ["GH18_e1"]


def test_overview_no_longer_regroups_already_resolved_hmm_hits(tmp_path):
    """Two dbCAN-HMM hits are never merged with each other by the overview.

    They only reach ``dbCAN_hmm_results.tsv`` after domain resolution accepted
    them, so re-running an overlap rule over them would undo that decision.
    """
    generator = _build_overview_dir(tmp_path, with_sidecars=False)
    # Spans overlap heavily (CBM18 lies inside GH18's span); the old HMM-vs-HMM
    # grouping loop collapsed these two to a single representative.
    annotations = [("GH18", 1, 400, "hmm"), ("CBM18", 300, 334, "hmm")]
    grouped = generator.graph_based_grouping(annotations)
    assert {g[0] for g in grouped} == {"GH18", "CBM18"}


# ---------------------------------------------------------------------------
# Phase 11 Test K: CGC downstream regression
# ---------------------------------------------------------------------------
def test_k_cgc_post_processing_keeps_multidomain_and_cross_class_annotations(tmp_path):
    # STP comes from pyhmmer and is already resolved: a catalytic domain plus an
    # accessory domain inside its insertion.
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
    # both STP domains survive (no legacy span re-filter of a resolved source)
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
    assert list(out.columns) == PU.CGC_SIG_RESULT_COLUMN + ["Type"]
    assert out.empty


# ---------------------------------------------------------------------------
# Phase 11 Test L: Pfam downstream regression
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


# ---------------------------------------------------------------------------
# Phase 3: legacy mode must be reachable for every pyhmmer library
# ---------------------------------------------------------------------------
def test_legacy_overlap_switch_reaches_every_pyhmmer_library():
    from dbcan.configs.pyhmmer_config import (
        DBCANSUBConfig,
        PyHMMERDBCANConfig,
        PyHMMERPfamConfig,
        PyHMMERSTPConfig,
        PyHMMERTFConfig,
    )
    from dbcan.parameter import create_config

    libraries = {
        "dbcan": PyHMMERDBCANConfig,
        "dbsub": DBCANSUBConfig,
        "tf": PyHMMERTFConfig,
        "stp": PyHMMERSTPConfig,
        "pfam": PyHMMERPfamConfig,
    }
    for namespace, config_class in libraries.items():
        default = create_config(config_class, namespace=namespace, db_dir="db", output_dir="out")
        assert default.matched_overlap is True, namespace
        legacy = create_config(
            config_class, namespace=namespace, db_dir="db", output_dir="out", matched_overlap=False
        )
        assert legacy.matched_overlap is False, namespace
