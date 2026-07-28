"""Regression tests for match-state-aware HMM domain resolution.

Covers the distinction the legacy filter collapsed:

    continuous alignment span != profile match-state evidence != domain conflict
"""

import itertools

import pandas as pd
import pytest

from dbcan.annotation.pyhmmer_search import PyHMMERProcessor
from dbcan.process.domain_resolution import (
    DomainHit,
    format_intervals,
    intervals_length,
    intervals_overlap,
    parse_intervals,
    resolve_domain_hits,
)
from dbcan.process.process_utils import filter_overlaps


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
class FakeAlignment:
    """Minimal stand-in for a pyhmmer DomainAlignment.

    HMMER renders an aligned model line where insert columns are ``.`` and an
    aligned target line where delete columns are ``-``; insert-aligned target
    residues are lowercase. Verified against pyhmmer 0.11.0 (see
    ``test_matched_intervals_against_real_pyhmmer_alignment``).
    """

    def __init__(self, hmm_sequence, target_sequence, target_from):
        self.hmm_sequence = hmm_sequence
        self.target_sequence = target_sequence
        self.target_from = target_from


def hit(name, ali_from, ali_to, matched=None, i_evalue=1.0, index=0, source=""):
    return DomainHit(
        hmm_name=name,
        target_name="prot1",
        source=source,
        ali_from=ali_from,
        ali_to=ali_to,
        matched_intervals=parse_intervals(matched, ali_from, ali_to),
        i_evalue=i_evalue,
        index=index,
    )


def resolved_names(hits, **kwargs):
    kwargs.setdefault("use_matched", True)
    return [h.hmm_name for h in resolve_domain_hits(hits, **kwargs)]


# ---------------------------------------------------------------------------
# Phase 7: interval arithmetic
# ---------------------------------------------------------------------------
def test_inclusive_interval_length_and_overlap():
    assert intervals_length([(1, 250), (340, 400)]) == 250 + 61
    assert intervals_overlap([(1, 10)], [(10, 20)]) == 1
    assert intervals_overlap([(1, 10)], [(11, 20)]) == 0
    assert intervals_overlap([(1, 250), (340, 400)], [(300, 334)]) == 0
    assert intervals_overlap([(1, 250), (340, 400)], [(200, 350)]) == 51 + 11


def test_intervals_overlap_is_symmetric_and_handles_many_intervals():
    a = [(1, 20), (41, 60), (81, 100)]
    b = [(15, 45), (55, 85)]
    assert intervals_overlap(a, b) == intervals_overlap(b, a)
    assert intervals_overlap(a, b) == 6 + 5 + 6 + 5


def test_parse_intervals_falls_back_to_alignment_span():
    assert parse_intervals(None, 10, 20) == ((10, 20),)
    assert parse_intervals("", 10, 20) == ((10, 20),)
    assert parse_intervals(float("nan"), 10, 20) == ((10, 20),)
    assert parse_intervals("garbage", 10, 20) == ((10, 20),)
    assert parse_intervals("1:250,340:400", 1, 400) == ((1, 250), (340, 400))
    assert format_intervals(((1, 250), (340, 400))) == "1:250,340:400"


# ---------------------------------------------------------------------------
# Phase 11 Test A: accessory domain inside an HMM insertion
# ---------------------------------------------------------------------------
def _test_a_hits():
    return [
        hit("GH18", 1, 400, "1:250,340:400", i_evalue=1e-120, index=0),
        hit("CBM18", 300, 334, "300:334", i_evalue=1e-30, index=1),
    ]


def test_a_matched_mode_keeps_accessory_domain_in_insertion():
    assert resolved_names(_test_a_hits()) == ["GH18", "CBM18"]


def test_a_legacy_mode_still_deletes_accessory_domain():
    assert resolved_names(_test_a_hits(), use_matched=False) == ["GH18"]


# ---------------------------------------------------------------------------
# Phase 11 Test B: true competing profiles
# ---------------------------------------------------------------------------
def test_b_true_competitors_resolve_to_one():
    hits = [
        hit("A", 100, 300, "100:300", i_evalue=1e-40, index=0),
        hit("B", 150, 280, "150:280", i_evalue=1e-20, index=1),
    ]
    assert resolved_names(hits) == ["A"]
    # evidence, not row order, decides the winner
    hits_b_better = [
        hit("A", 100, 300, "100:300", i_evalue=1e-20, index=0),
        hit("B", 150, 280, "150:280", i_evalue=1e-40, index=1),
    ]
    assert resolved_names(hits_b_better) == ["B"]


# ---------------------------------------------------------------------------
# Phase 11 Test C: overlapping ali spans, zero matched overlap
# ---------------------------------------------------------------------------
def test_c_overlapping_spans_with_zero_matched_overlap_keep_both():
    hits = [
        hit("A", 1, 500, "1:200,400:500", i_evalue=1e-90, index=0),
        hit("B", 250, 350, "250:350", i_evalue=1e-10, index=1),
    ]
    assert resolved_names(hits) == ["A", "B"]
    assert resolved_names(hits, use_matched=False) == ["A"]


# ---------------------------------------------------------------------------
# Phase 11 Tests D/E + Phase 12: PyHMMER alignment -> intervals
# ---------------------------------------------------------------------------
def test_e_consecutive_insert_state_residues_are_excluded():
    aln = FakeAlignment("ABC.....DEF", "ABCxyzuvDEF", 1)
    # 5 model insert columns but only 5 inserted residues (xyzuv) at 4-8
    assert PyHMMERProcessor._matched_interval_pairs(aln) == [(1, 3), (9, 11)]
    assert PyHMMERProcessor._matched_intervals(aln) == "1:3,9:11"


def test_d_delete_state_columns_consume_no_target_coordinate():
    # Columns 4-6 are delete states: the model has consensus positions there,
    # the target has none, so target coordinates stay contiguous.
    aln = FakeAlignment("ABCDEFGHI", "ABC---DEF", 1)
    assert PyHMMERProcessor._matched_interval_pairs(aln) == [(1, 6)]


def test_d_deletes_and_inserts_combined():
    aln = FakeAlignment("ABCDEF..GHI", "ABC--Fxy GHI".replace(" ", ""), 10)
    # target residues: A(10) B(11) C(12) F(13) x(14) y(15) G(16) H(17) I(18)
    # '-' at model D,E consume nothing; x,y sit in insert columns
    assert PyHMMERProcessor._matched_interval_pairs(aln) == [(10, 13), (16, 18)]


def test_matched_intervals_of_a_fully_matched_alignment_is_the_span():
    aln = FakeAlignment("ABCDEF", "ABCDEF", 21)
    assert PyHMMERProcessor._matched_interval_pairs(aln) == [(21, 26)]


def test_empty_alignment_yields_no_intervals():
    aln = FakeAlignment("", "", 1)
    assert PyHMMERProcessor._matched_interval_pairs(aln) == []
    assert PyHMMERProcessor._matched_intervals(aln) == ""


def test_matched_intervals_against_real_pyhmmer_alignment():
    """Phase 12: exercise the real pyhmmer alignment representation.

    Builds a tiny profile from an in-memory MSA and searches a target that
    carries a 10-residue insertion, so HMMER must emit insert-state columns.
    No external database is involved.
    """
    pyhmmer = pytest.importorskip("pyhmmer")
    from pyhmmer.easel import Alphabet, DigitalSequenceBlock, TextMSA, TextSequence

    alphabet = Alphabet.amino()
    core = "MKAILVVLLYTFATANADTLCIGYHANNSTDTVDTVLEKNVTVTHSVNLLEDKHNGKLCKL"
    msa = TextMSA(
        name=b"toy",
        sequences=[
            TextSequence(name=f"s{i}".encode(), sequence=core[: 3 + i] + "V" + core[4 + i:])
            for i in range(6)
        ],
    ).digitize(alphabet)
    hmm, _, _ = pyhmmer.plan7.Builder(alphabet).build_msa(
        msa, pyhmmer.plan7.Background(alphabet)
    )

    half = len(core) // 2
    insertion = "GSGSGSGSGS"
    target = core[:half] + insertion + core[half:]
    block = DigitalSequenceBlock(
        alphabet, [TextSequence(name=b"tgt", sequence=target).digitize(alphabet)]
    )

    alignments = [
        domain.alignment
        for hits in pyhmmer.hmmsearch([hmm], block, cpus=1, E=1000.0, domE=1000.0)
        for hit_ in hits
        for domain in hit_.domains
    ]
    assert len(alignments) == 1, "expected one domain spanning the insertion"
    aln = alignments[0]

    assert "." in aln.hmm_sequence, "expected insert-state columns in the model line"
    intervals = PyHMMERProcessor._matched_interval_pairs(aln)

    # The continuous span covers the insertion; the match-state intervals do not.
    inserted = set(range(half + 1, half + len(insertion) + 1))
    assert aln.target_from <= min(inserted) and aln.target_to >= max(inserted)
    covered = {pos for start, end in intervals for pos in range(start, end + 1)}
    assert not (covered & inserted)
    assert intervals == [(1, half), (half + len(insertion) + 1, len(target))]


# ---------------------------------------------------------------------------
# Phase 11 Test F: exact threshold boundary
# ---------------------------------------------------------------------------
def test_f_exactly_half_matched_overlap_is_not_a_conflict():
    """The historical rule is strict ``ratio > 0.5``; exactly 50% is kept."""
    # B is 100 residues long and shares exactly 50 of them with A.
    hits = [
        hit("A", 1, 200, "1:200", i_evalue=1e-50, index=0),
        hit("B", 151, 250, "151:250", i_evalue=1e-10, index=1),
    ]
    assert intervals_overlap([(1, 200)], [(151, 250)]) == 50
    assert resolved_names(hits) == ["A", "B"]

    # One residue more and it becomes a conflict.
    hits_over = [
        hit("A", 1, 201, "1:201", i_evalue=1e-50, index=0),
        hit("B", 151, 250, "151:250", i_evalue=1e-10, index=1),
    ]
    assert resolved_names(hits_over) == ["A"]


# ---------------------------------------------------------------------------
# Phase 11 Test G: input-order invariance
# ---------------------------------------------------------------------------
def test_g_resolution_is_independent_of_input_order():
    base = [
        hit("A", 100, 300, "100:300", i_evalue=1e-50),
        hit("B", 600, 800, "600:800", i_evalue=1e-40),
        hit("C", 150, 280, "150:280", i_evalue=1e-30),
    ]
    outcomes = set()
    for permutation in itertools.permutations(range(3)):
        hits = [
            DomainHit(
                hmm_name=base[src].hmm_name,
                target_name=base[src].target_name,
                ali_from=base[src].ali_from,
                ali_to=base[src].ali_to,
                matched_intervals=base[src].matched_intervals,
                i_evalue=base[src].i_evalue,
                index=position,
            )
            for position, src in enumerate(permutation)
        ]
        outcomes.add(frozenset(h.hmm_name for h in resolve_domain_hits(hits, use_matched=True)))
    assert outcomes == {frozenset({"A", "B"})}


def test_g_non_adjacent_conflict_is_not_missed():
    """The old ``keep[-1]`` chain only ever compared against the last kept hit.

    A conflicts with C, B does not. With rows ordered A, B, C the legacy filter
    compared C against B only and kept all three.
    """
    hits = [
        hit("A", 100, 300, "100:300", i_evalue=1e-50, index=0),
        hit("B", 305, 400, "305:400", i_evalue=1e-40, index=1),
        hit("C", 120, 290, "120:290", i_evalue=1e-30, index=2),
    ]
    assert resolved_names(hits) == ["A", "B"]


def test_non_conflicting_third_hit_survives_its_component():
    """A-B conflict and B-C conflict, but A and C do not touch each other.

    All three land in one connected component; resolving it best-evidence
    first keeps both non-conflicting ends instead of only one winner.
    """
    hits = [
        hit("A", 1, 100, "1:100", i_evalue=1e-60, index=0),
        hit("B", 45, 144, "45:144", i_evalue=1e-10, index=1),
        hit("C", 90, 189, "90:189", i_evalue=1e-50, index=2),
    ]
    assert resolved_names(hits) == ["A", "C"]


# ---------------------------------------------------------------------------
# Phase 11 Test H: coexistence policy
# ---------------------------------------------------------------------------
def test_h_coexistence_policy_keeps_genuinely_overlapping_pair():
    hits = [
        hit("GH189.hmm", 1251, 1529, "1251:1529", i_evalue=8.7e-133, index=0),
        hit("GT84.hmm", 1301, 1515, "1301:1515", i_evalue=2.8e-104, index=1),
    ]
    assert resolved_names(hits) == ["GH189.hmm"]
    assert resolved_names(hits, coexistence_pairs={frozenset({"GT84", "GH189"})}) == [
        "GH189.hmm",
        "GT84.hmm",
    ]


def test_coexistence_is_pair_specific_not_a_blanket_exemption():
    hits = [
        hit("GH192.hmm", 1251, 1529, "1251:1529", i_evalue=8.7e-133, index=0),
        hit("GT84.hmm", 1301, 1515, "1301:1515", i_evalue=2.8e-104, index=1),
    ]
    assert resolved_names(hits, coexistence_pairs={frozenset({"GT84", "GH189"})}) == ["GH192.hmm"]


# ---------------------------------------------------------------------------
# Source awareness
# ---------------------------------------------------------------------------
def test_hits_from_different_sources_do_not_compete_by_default():
    hits = [
        hit("STP_kinase", 1, 400, "1:400", i_evalue=1e-50, index=0, source="STP"),
        hit("TC_transporter", 10, 390, "10:390", i_evalue=1e-10, index=1, source="TC"),
    ]
    assert resolved_names(hits) == ["STP_kinase", "TC_transporter"]
    assert resolved_names(hits, cross_source_conflicts=True) == ["STP_kinase"]


# ---------------------------------------------------------------------------
# DataFrame adapter parity
# ---------------------------------------------------------------------------
def _frame(rows):
    return pd.DataFrame(rows)


def test_filter_overlaps_matches_the_resolver_and_is_order_invariant():
    rows = [
        {
            "HMM Name": "GH18.hmm",
            "Target Name": "prot1",
            "Target From": 1,
            "Target To": 400,
            "i-Evalue": 1e-120,
            "Matched Intervals": "1:250,340:400",
        },
        {
            "HMM Name": "CBM18.hmm",
            "Target Name": "prot1",
            "Target From": 300,
            "Target To": 334,
            "i-Evalue": 1e-30,
            "Matched Intervals": "300:334",
        },
    ]
    forward = filter_overlaps(_frame(rows), matched_overlap=True)
    reverse = filter_overlaps(_frame(list(reversed(rows))), matched_overlap=True)
    assert list(forward["HMM Name"]) == ["GH18.hmm", "CBM18.hmm"]
    assert list(reverse["HMM Name"]) == ["GH18.hmm", "CBM18.hmm"]


def test_filter_overlaps_does_not_mutate_the_caller_frame():
    df = _frame(
        [
            {
                "HMM Name": "A.hmm",
                "Target Name": "p",
                "Target From": "1",
                "Target To": "100",
                "i-Evalue": "1e-5",
            }
        ]
    )
    filter_overlaps(df)
    assert df["Target From"].iloc[0] == "1"


def test_filter_overlaps_returns_input_when_required_columns_are_missing():
    df = _frame([{"HMM Name": "A.hmm"}])
    assert filter_overlaps(df) is df
