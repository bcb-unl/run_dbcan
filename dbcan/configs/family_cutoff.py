"""Fail-closed loader for opt-in dbCAN family-specific domain cutoffs."""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping

import pyhmmer


ALLOWED_DECISION_STATUSES = {"approved", "candidate", "provisional", "exploratory"}


@dataclass(frozen=True)
class FamilyCutoffRule:
    min_domain_bitscore: float
    min_hmm_coverage: float
    decision_status: str


@dataclass(frozen=True)
class DbcanFamilyCutoffPolicy:
    profile_id: str
    deployment_status: str
    search_domain_evalue_max: float
    dbcan_hmm_sha256: str
    hmm_profile_count: int
    rules: Mapping[str, FamilyCutoffRule]
    source_path: Path
    source_sha256: str


@dataclass(frozen=True)
class DbcanOverlapCoexistencePolicy:
    policy_id: str
    deployment_status: str
    dbcan_hmm_sha256: str
    hmm_profile_count: int
    pairs: frozenset[frozenset[str]]
    source_path: Path
    source_sha256: str


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _hmm_names(path: Path) -> set[str]:
    names: set[str] = set()
    with pyhmmer.plan7.HMMFile(str(path)) as handle:
        for hmm in handle:
            if hmm.name is None:
                raise ValueError(f"HMM without NAME in {path}")
            name = hmm.name.decode("utf-8") if isinstance(hmm.name, bytes) else str(hmm.name)
            if name in names:
                raise ValueError(f"duplicate HMM NAME in library: {name}")
            names.add(name)
    return names


def _finite_number(value, field: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{field} must be a finite number")
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{field} must be a finite number") from exc
    if not math.isfinite(number):
        raise ValueError(f"{field} must be a finite number")
    return number


def load_dbcan_family_cutoff_policy(
    policy_path: str | Path,
    hmm_path: str | Path,
    *,
    allow_experimental: bool = False,
) -> DbcanFamilyCutoffPolicy:
    """Load and validate a policy against the exact HMM library in use."""

    source_path = Path(policy_path).expanduser().resolve()
    library_path = Path(hmm_path).expanduser().resolve()
    if not source_path.is_file():
        raise FileNotFoundError(f"dbCAN family cutoff policy not found: {source_path}")
    if not library_path.is_file():
        raise FileNotFoundError(f"dbCAN HMM library not found: {library_path}")

    try:
        payload = json.loads(source_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"invalid dbCAN family cutoff JSON: {source_path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise ValueError("dbCAN family cutoff policy root must be an object")
    if payload.get("schema_version") != 1:
        raise ValueError("dbCAN family cutoff policy schema_version must equal 1")

    profile_id = payload.get("profile_id")
    if not isinstance(profile_id, str) or not profile_id.strip():
        raise ValueError("profile_id must be a non-empty string")
    deployment_status = payload.get("deployment_status")
    if deployment_status not in {"experimental", "approved"}:
        raise ValueError("deployment_status must be experimental or approved")
    if deployment_status == "experimental" and not allow_experimental:
        raise ValueError(
            "experimental dbCAN cutoff policy requires --allow-experimental-dbcan-cutoffs"
        )

    expected_digest = payload.get("dbcan_hmm_sha256")
    if not isinstance(expected_digest, str) or len(expected_digest) != 64:
        raise ValueError("dbcan_hmm_sha256 must be a 64-character SHA-256 hex digest")
    actual_digest = sha256_file(library_path)
    if actual_digest.lower() != expected_digest.lower():
        raise ValueError(
            f"dbCAN HMM SHA-256 mismatch: policy={expected_digest} actual={actual_digest}"
        )

    names = _hmm_names(library_path)
    expected_count = payload.get("hmm_profile_count")
    if isinstance(expected_count, bool) or not isinstance(expected_count, int):
        raise ValueError("hmm_profile_count must be an integer")
    if expected_count != len(names):
        raise ValueError(
            f"dbCAN HMM profile count mismatch: policy={expected_count} actual={len(names)}"
        )

    search_evalue = _finite_number(
        payload.get("search_domain_evalue_max"), "search_domain_evalue_max"
    )
    if not 0.0 < search_evalue <= 10.0:
        raise ValueError("search_domain_evalue_max must be in (0, 10]")

    raw_rules = payload.get("families")
    if not isinstance(raw_rules, dict) or not raw_rules:
        raise ValueError("families must be a non-empty object")
    unknown = sorted(set(raw_rules) - names)
    if unknown:
        raise ValueError(f"unknown exact HMM NAME(s) in cutoff policy: {', '.join(unknown)}")

    rules: dict[str, FamilyCutoffRule] = {}
    for name, raw_rule in raw_rules.items():
        if not isinstance(raw_rule, dict):
            raise ValueError(f"family rule must be an object: {name}")
        bitscore = _finite_number(raw_rule.get("min_domain_bitscore"), f"{name}.min_domain_bitscore")
        coverage = _finite_number(raw_rule.get("min_hmm_coverage"), f"{name}.min_hmm_coverage")
        status = raw_rule.get("decision_status")
        if bitscore <= 0.0:
            raise ValueError(f"{name}.min_domain_bitscore must be > 0")
        if not 0.0 <= coverage <= 1.0:
            raise ValueError(f"{name}.min_hmm_coverage must be in [0, 1]")
        if status not in ALLOWED_DECISION_STATUSES:
            raise ValueError(
                f"{name}.decision_status must be one of {sorted(ALLOWED_DECISION_STATUSES)}"
            )
        rules[name] = FamilyCutoffRule(bitscore, coverage, status)

    return DbcanFamilyCutoffPolicy(
        profile_id=profile_id,
        deployment_status=deployment_status,
        search_domain_evalue_max=search_evalue,
        dbcan_hmm_sha256=actual_digest,
        hmm_profile_count=len(names),
        rules=rules,
        source_path=source_path,
        source_sha256=sha256_file(source_path),
    )


def _base_family(hmm_name: str) -> str:
    name = hmm_name[:-4] if hmm_name.endswith(".hmm") else hmm_name
    return name.split("_", 1)[0]


def load_dbcan_overlap_coexistence_policy(
    policy_path: str | Path,
    hmm_path: str | Path,
    *,
    allow_experimental: bool = False,
) -> DbcanOverlapCoexistencePolicy:
    """Load an opt-in unordered base-family coexistence allowlist."""

    source_path = Path(policy_path).expanduser().resolve()
    library_path = Path(hmm_path).expanduser().resolve()
    if not source_path.is_file():
        raise FileNotFoundError(f"dbCAN overlap coexistence policy not found: {source_path}")
    if not library_path.is_file():
        raise FileNotFoundError(f"dbCAN HMM library not found: {library_path}")
    try:
        payload = json.loads(source_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"invalid dbCAN overlap coexistence JSON: {source_path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise ValueError("dbCAN overlap coexistence policy root must be an object")
    if payload.get("schema_version") != 1:
        raise ValueError("dbCAN overlap coexistence schema_version must equal 1")

    policy_id = payload.get("policy_id")
    if not isinstance(policy_id, str) or not policy_id.strip():
        raise ValueError("policy_id must be a non-empty string")
    deployment_status = payload.get("deployment_status")
    if deployment_status not in {"experimental", "approved"}:
        raise ValueError("deployment_status must be experimental or approved")
    if deployment_status == "experimental" and not allow_experimental:
        raise ValueError(
            "experimental dbCAN overlap policy requires --allow-experimental-dbcan-cutoffs"
        )

    expected_digest = payload.get("dbcan_hmm_sha256")
    if not isinstance(expected_digest, str) or len(expected_digest) != 64:
        raise ValueError("dbcan_hmm_sha256 must be a 64-character SHA-256 hex digest")
    actual_digest = sha256_file(library_path)
    if actual_digest.lower() != expected_digest.lower():
        raise ValueError(
            f"dbCAN HMM SHA-256 mismatch: policy={expected_digest} actual={actual_digest}"
        )

    names = _hmm_names(library_path)
    expected_count = payload.get("hmm_profile_count")
    if isinstance(expected_count, bool) or not isinstance(expected_count, int):
        raise ValueError("hmm_profile_count must be an integer")
    if expected_count != len(names):
        raise ValueError(
            f"dbCAN HMM profile count mismatch: policy={expected_count} actual={len(names)}"
        )
    known_bases = {_base_family(name) for name in names}

    raw_pairs = payload.get("pairs")
    if not isinstance(raw_pairs, list) or not raw_pairs:
        raise ValueError("pairs must be a non-empty array")
    pairs: set[frozenset[str]] = set()
    for index, raw_pair in enumerate(raw_pairs):
        if not isinstance(raw_pair, dict):
            raise ValueError(f"pairs[{index}] must be an object")
        family_a = raw_pair.get("family_a")
        family_b = raw_pair.get("family_b")
        status = raw_pair.get("decision_status")
        if not isinstance(family_a, str) or not isinstance(family_b, str):
            raise ValueError(f"pairs[{index}] family_a/family_b must be strings")
        family_a, family_b = _base_family(family_a), _base_family(family_b)
        if family_a == family_b:
            raise ValueError(f"pairs[{index}] cannot contain the same base family twice")
        unknown = sorted({family_a, family_b} - known_bases)
        if unknown:
            raise ValueError(f"pairs[{index}] contains unknown base family: {', '.join(unknown)}")
        if status not in ALLOWED_DECISION_STATUSES:
            raise ValueError(
                f"pairs[{index}].decision_status must be one of "
                f"{sorted(ALLOWED_DECISION_STATUSES)}"
            )
        pair = frozenset((family_a, family_b))
        if pair in pairs:
            raise ValueError(f"duplicate unordered overlap coexistence pair: {family_a}, {family_b}")
        pairs.add(pair)

    return DbcanOverlapCoexistencePolicy(
        policy_id=policy_id,
        deployment_status=deployment_status,
        dbcan_hmm_sha256=actual_digest,
        hmm_profile_count=len(names),
        pairs=frozenset(pairs),
        source_path=source_path,
        source_sha256=sha256_file(source_path),
    )
