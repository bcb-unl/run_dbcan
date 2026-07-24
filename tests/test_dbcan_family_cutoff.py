import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from dbcan.annotation.pyhmmer_search import PyHMMERDBCANProcessor
from dbcan.configs.family_cutoff import (
    DbcanFamilyCutoffPolicy,
    FamilyCutoffRule,
    load_dbcan_family_cutoff_policy,
    load_dbcan_overlap_coexistence_policy,
)
from dbcan.configs.pyhmmer_config import PyHMMERDBCANConfig
from dbcan.parameter import create_config


def _policy_payload(**updates):
    payload = {
        "schema_version": 1,
        "profile_id": "test-policy",
        "deployment_status": "experimental",
        "dbcan_hmm_sha256": "a" * 64,
        "hmm_profile_count": 2,
        "search_domain_evalue_max": 1e-3,
        "families": {
            "CBM1.hmm": {
                "min_domain_bitscore": 15.871,
                "min_hmm_coverage": 0.4,
                "decision_status": "candidate",
            }
        },
    }
    payload.update(updates)
    return payload


def _write_policy(tmp_path, payload):
    policy_path = tmp_path / "policy.json"
    policy_path.write_text(json.dumps(payload), encoding="utf-8")
    hmm_path = tmp_path / "dbCAN.hmm"
    hmm_path.write_text("test HMM placeholder\n", encoding="utf-8")
    return policy_path, hmm_path


def _mock_library(monkeypatch):
    monkeypatch.setattr(
        "dbcan.configs.family_cutoff._hmm_names",
        lambda path: {"CBM1.hmm", "GH5.hmm"},
    )
    monkeypatch.setattr("dbcan.configs.family_cutoff.sha256_file", lambda path: "a" * 64)


def test_experimental_policy_requires_explicit_authorization(tmp_path, monkeypatch):
    _mock_library(monkeypatch)
    policy_path, hmm_path = _write_policy(tmp_path, _policy_payload())

    with pytest.raises(ValueError, match="allow-experimental"):
        load_dbcan_family_cutoff_policy(policy_path, hmm_path)

    policy = load_dbcan_family_cutoff_policy(
        policy_path,
        hmm_path,
        allow_experimental=True,
    )
    assert policy.profile_id == "test-policy"
    assert policy.search_domain_evalue_max == pytest.approx(1e-3)
    assert policy.rules["CBM1.hmm"].min_domain_bitscore == pytest.approx(15.871)


@pytest.mark.parametrize(
    ("updates", "message"),
    [
        ({"dbcan_hmm_sha256": "b" * 64}, "SHA-256 mismatch"),
        ({"hmm_profile_count": 3}, "profile count mismatch"),
        (
            {
                "families": {
                    "UNKNOWN.hmm": {
                        "min_domain_bitscore": 10,
                        "min_hmm_coverage": 0.3,
                        "decision_status": "candidate",
                    }
                }
            },
            "unknown exact HMM NAME",
        ),
    ],
)
def test_policy_mismatch_fails_closed(tmp_path, monkeypatch, updates, message):
    _mock_library(monkeypatch)
    policy_path, hmm_path = _write_policy(tmp_path, _policy_payload(**updates))
    with pytest.raises(ValueError, match=message):
        load_dbcan_family_cutoff_policy(policy_path, hmm_path, allow_experimental=True)


def test_overlap_policy_normalizes_subfamilies_and_rejects_duplicates(tmp_path, monkeypatch):
    _mock_library(monkeypatch)
    payload = {
        "schema_version": 1,
        "policy_id": "overlap-test",
        "deployment_status": "experimental",
        "dbcan_hmm_sha256": "a" * 64,
        "hmm_profile_count": 2,
        "pairs": [
            {
                "family_a": "CBM1.hmm",
                "family_b": "GH5_1.hmm",
                "decision_status": "candidate",
            }
        ],
    }
    policy_path, hmm_path = _write_policy(tmp_path, payload)
    policy = load_dbcan_overlap_coexistence_policy(
        policy_path,
        hmm_path,
        allow_experimental=True,
    )
    assert policy.pairs == frozenset({frozenset({"CBM1", "GH5"})})

    payload["pairs"].append(
        {"family_a": "GH5", "family_b": "CBM1", "decision_status": "candidate"}
    )
    policy_path.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(ValueError, match="duplicate unordered"):
        load_dbcan_overlap_coexistence_policy(
            policy_path,
            hmm_path,
            allow_experimental=True,
        )


def _processor(tmp_path):
    config = PyHMMERDBCANConfig(db_dir=str(tmp_path), output_dir=str(tmp_path))
    return PyHMMERDBCANProcessor(config)


def test_no_policy_preserves_legacy_strict_boundaries(tmp_path):
    processor = _processor(tmp_path)
    assert processor.search_domain_evalue_threshold == pytest.approx(1e-15)
    assert processor.hmmsearch_extra_options == {}
    assert processor._passes_domain("GH5.hmm", 1e-20, 0.36, 0.0)
    assert not processor._passes_domain("GH5.hmm", 1e-15, 0.36, 999.0)
    assert not processor._passes_domain("GH5.hmm", 1e-20, 0.35, 999.0)
    domains_with_included = SimpleNamespace(included=["reported"])
    hit = SimpleNamespace(domains=domains_with_included)
    assert processor._domains_for_hit(hit) == ["reported"]


def test_family_policy_is_inclusive_and_other_families_remain_legacy(tmp_path):
    processor = _processor(tmp_path)
    processor.family_cutoff_policy = DbcanFamilyCutoffPolicy(
        profile_id="test",
        deployment_status="experimental",
        search_domain_evalue_max=1e-3,
        dbcan_hmm_sha256="a" * 64,
        hmm_profile_count=2,
        rules={"CBM1.hmm": FamilyCutoffRule(15.871, 0.4, "candidate")},
        source_path=Path("policy.json"),
        source_sha256="b" * 64,
    )

    assert processor.search_domain_evalue_threshold == pytest.approx(1e-3)
    assert processor.hmmsearch_extra_options == {}
    domains = ["reported", "tail"]
    hit = SimpleNamespace(domains=domains)
    assert processor._domains_for_hit(hit) is domains
    assert processor._passes_domain("CBM1.hmm", 1.0, 0.4, 15.871)
    assert not processor._passes_domain("CBM1.hmm", 1e-100, 0.4, 15.870)
    assert not processor._passes_domain("CBM1.hmm", 1e-100, 0.399, 999.0)
    assert processor._passes_domain("GH5.hmm", 1e-20, 0.36, 0.0)


def test_namespaced_cli_values_reach_dbcan_config():
    config = create_config(
        PyHMMERDBCANConfig,
        namespace="dbcan",
        db_dir="db",
        output_dir="out",
        family_cutoff_file_dbcan="cutoffs.json",
        overlap_coexistence_file_dbcan="overlap.json",
        allow_experimental_dbcan_cutoffs_dbcan=True,
    )
    assert config.family_cutoff_file == "cutoffs.json"
    assert config.overlap_coexistence_file == "overlap.json"
    assert config.allow_experimental_dbcan_cutoffs is True
