"""Regression tests for the conflict-graph HMM domain resolution.

Overlap is measured on the alignment span with the 5.2.9 arithmetic; what
changed from 5.2.9 is that conflicts are resolved on a per-protein conflict
graph instead of the ``keep[-1]`` chain, so the result no longer depends on row
order.
"""

import itertools

import pandas as pd

from dbcan.process.domain_resolution import (
    DomainHit,
    geometric_overlap,
    resolve_domain_hits,
)
from dbcan.process.process_utils import filter_overlaps


def hit(name, ali_from, ali_to, i_evalue=1.0, index=0, source=""):
    return DomainHit(
        hmm_name=name,
        target_name="prot1",
        source=source,
        ali_from=ali_from,
        ali_to=ali_to,
        i_evalue=i_evalue,
        index=index,
    )


def resolved_names(hits, **kwargs):
    return [h.hmm_name for h in resolve_domain_hits(hits, **kwargs)]


# ---------------------------------------------------------------------------
# Geometry: alignment span, 5.2.9 arithmetic
# ---------------------------------------------------------------------------
def test_geometry_keeps_the_5_2_9_off_by_one():
    # overlap = min(to) - max(from); length = to - from (no +1 on either)
    assert geometric_overlap(hit("A", 1, 201), hit("B", 151, 251)) == (50, 200, 100)


def test_exactly_half_is_not_a_conflict_under_the_5_2_9_arithmetic():
    """The rule is strict ``ratio > 0.5``; exactly 50% (as 5.2.9 counts it) is kept."""
    hits = [hit("A", 1, 201, i_evalue=1e-50, index=0), hit("B", 151, 251, i_evalue=1e-10, index=1)]
    assert resolved_names(hits) == ["A", "B"]
    # B one residue shorter: 50 / 99 > 0.5, a conflict
    hits_over = [hit("A", 1, 201, i_evalue=1e-50, index=0), hit("B", 151, 250, i_evalue=1e-10, index=1)]
    assert resolved_names(hits_over) == ["A"]


def test_accessory_domain_inside_a_long_span_is_removed():
    """The span rule (as in 5.2.9): a CBM inside a long catalytic alignment conflicts with it."""
    hits = [hit("GH18", 1, 400, i_evalue=1e-120, index=0), hit("CBM18", 300, 334, i_evalue=1e-30, index=1)]
    assert resolved_names(hits) == ["GH18"]


def test_true_competitors_resolve_by_evidence_not_row_order():
    hits = [hit("A", 100, 300, i_evalue=1e-40, index=0), hit("B", 150, 280, i_evalue=1e-20, index=1)]
    assert resolved_names(hits) == ["A"]
    hits_b_better = [hit("A", 100, 300, i_evalue=1e-20, index=0), hit("B", 150, 280, i_evalue=1e-40, index=1)]
    assert resolved_names(hits_b_better) == ["B"]


# ---------------------------------------------------------------------------
# Conflict graph vs the 5.2.9 keep[-1] chain
# ---------------------------------------------------------------------------
def test_kinase_split_by_an_insert_keeps_both_halves():
    """Xylona heveae KZF26007.1 (STP.hmm): Pkinase is split by an insert into two halves.

    Rows in Target From order: Pkinase 34-161, Pkinase_Tyr 34-364, Pkinase 215-363.
    The 5.2.9 chain replaced the first half by Pkinase_Tyr (better i-Evalue), then
    Pkinase_Tyr by the second half, and so lost the first half, which conflicts
    only with the discarded Pkinase_Tyr.
    """
    hits = [
        hit("Pkinase", 34, 161, i_evalue=1.21e-15, index=0),
        hit("Pkinase_Tyr", 34, 364, i_evalue=1.01e-16, index=1),
        hit("Pkinase", 215, 363, i_evalue=3.10e-22, index=2),
    ]
    kept = resolve_domain_hits(hits)
    assert [(h.hmm_name, h.ali_from) for h in kept] == [("Pkinase", 34), ("Pkinase", 215)]


def test_resolution_is_independent_of_input_order():
    base = [
        hit("A", 100, 300, i_evalue=1e-50),
        hit("B", 600, 800, i_evalue=1e-40),
        hit("C", 150, 280, i_evalue=1e-30),
    ]
    outcomes = set()
    for permutation in itertools.permutations(range(3)):
        hits = [
            DomainHit(
                hmm_name=base[src].hmm_name,
                target_name=base[src].target_name,
                ali_from=base[src].ali_from,
                ali_to=base[src].ali_to,
                i_evalue=base[src].i_evalue,
                index=position,
            )
            for position, src in enumerate(permutation)
        ]
        outcomes.add(frozenset(h.hmm_name for h in resolve_domain_hits(hits)))
    assert outcomes == {frozenset({"A", "B"})}


def test_non_adjacent_conflict_is_not_missed():
    """The ``keep[-1]`` chain compared C only with B (the last kept hit) and kept all three."""
    hits = [
        hit("A", 100, 300, i_evalue=1e-50, index=0),
        hit("B", 305, 400, i_evalue=1e-40, index=1),
        hit("C", 120, 290, i_evalue=1e-30, index=2),
    ]
    assert resolved_names(hits) == ["A", "B"]


def test_non_conflicting_third_hit_survives_its_component():
    """A-B and B-C conflict, A and C do not: both ends are kept."""
    hits = [
        hit("A", 1, 100, i_evalue=1e-60, index=0),
        hit("B", 45, 144, i_evalue=1e-10, index=1),
        hit("C", 90, 189, i_evalue=1e-50, index=2),
    ]
    assert resolved_names(hits) == ["A", "C"]


# ---------------------------------------------------------------------------
# Coexistence policy
# ---------------------------------------------------------------------------
def test_coexistence_policy_keeps_an_overlapping_pair():
    hits = [
        hit("GH189.hmm", 1251, 1529, i_evalue=8.7e-133, index=0),
        hit("GT84.hmm", 1301, 1515, i_evalue=2.8e-104, index=1),
    ]
    assert resolved_names(hits) == ["GH189.hmm"]
    assert resolved_names(hits, coexistence_pairs={frozenset({"GT84", "GH189"})}) == [
        "GH189.hmm",
        "GT84.hmm",
    ]


def test_coexistence_is_pair_specific_not_a_blanket_exemption():
    hits = [
        hit("GH192.hmm", 1251, 1529, i_evalue=8.7e-133, index=0),
        hit("GT84.hmm", 1301, 1515, i_evalue=2.8e-104, index=1),
    ]
    assert resolved_names(hits, coexistence_pairs={frozenset({"GT84", "GH189"})}) == ["GH192.hmm"]


# ---------------------------------------------------------------------------
# Source awareness
# ---------------------------------------------------------------------------
def test_hits_from_different_sources_do_not_compete_by_default():
    hits = [
        hit("STP_kinase", 1, 400, i_evalue=1e-50, index=0, source="STP"),
        hit("TC_transporter", 10, 390, i_evalue=1e-10, index=1, source="TC"),
    ]
    assert resolved_names(hits) == ["STP_kinase", "TC_transporter"]
    assert resolved_names(hits, cross_source_conflicts=True) == ["STP_kinase"]


# ---------------------------------------------------------------------------
# DataFrame adapter
# ---------------------------------------------------------------------------
def _kinase_rows():
    return [
        {"HMM Name": "Pkinase", "Target Name": "KZF26007.1", "Target From": 34, "Target To": 161,
         "i-Evalue": 1.21e-15},
        {"HMM Name": "Pkinase_Tyr", "Target Name": "KZF26007.1", "Target From": 34, "Target To": 364,
         "i-Evalue": 1.01e-16},
        {"HMM Name": "Pkinase", "Target Name": "KZF26007.1", "Target From": 215, "Target To": 363,
         "i-Evalue": 3.10e-22},
    ]


def test_filter_overlaps_matches_the_resolver_and_is_order_invariant():
    rows = _kinase_rows()
    forward = filter_overlaps(pd.DataFrame(rows))
    reverse = filter_overlaps(pd.DataFrame(list(reversed(rows))))
    assert sorted(forward["Target From"]) == [34, 215]
    assert sorted(reverse["Target From"]) == [34, 215]
    assert set(forward["HMM Name"]) == {"Pkinase"}


def test_filter_overlaps_does_not_mutate_the_caller_frame():
    df = pd.DataFrame(
        [{"HMM Name": "A.hmm", "Target Name": "p", "Target From": "1", "Target To": "100", "i-Evalue": "1e-5"}]
    )
    filter_overlaps(df)
    assert df["Target From"].iloc[0] == "1"


def test_filter_overlaps_returns_input_when_required_columns_are_missing():
    df = pd.DataFrame([{"HMM Name": "A.hmm"}])
    assert filter_overlaps(df) is df
