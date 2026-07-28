"""Centralised HMM domain-resolution logic.

Three concepts are deliberately kept apart here, because conflating them is
exactly what the historical overlap filter got wrong:

``alignment geometry``
    The continuous ``Target From``/``Target To`` span HMMER reports for a
    domain alignment.

``match-state evidence``
    The subset of that span whose residues are aligned to a profile *match*
    state.  A target residue aligned to an HMM *insert* state sits between two
    consensus match positions: it is inside the alignment span, but the profile
    carries no consensus-column evidence over it.  A long catalytic profile can
    therefore geometrically span an inserted accessory domain (a CBM, say)
    without having any match-state evidence there.

``biological domain conflict``
    Two models genuinely competing for the same conserved region.  Geometric
    overlap is only the *evidence* for a conflict; an explicit coexistence
    policy can still allow two substantially overlapping models to be reported
    together, and hits from unrelated annotation sources never compete at all.

Note that an insert state is an HMM alignment concept, not a biological one: a
lineage-specific loop or an inserted structural region may still belong to the
domain biologically.  Match-state occupancy is used here because it is the
evidence the *profile* actually provides, not because inserted residues are
"not part of the domain".
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Iterable, List, Optional, Sequence, Set, Tuple

Interval = Tuple[int, int]
Intervals = Tuple[Interval, ...]

#: Default fraction of the shorter domain that must be covered before two hits
#: are treated as competitors.  The historical rule is a strict ``>``: an
#: overlap of exactly this fraction is *not* a conflict.
DEFAULT_OVERLAP_RATIO_THRESHOLD = 0.5

_INTERVAL_SEP = ","
_BOUND_SEP = ":"


# ---------------------------------------------------------------------------
# Interval arithmetic (1-based, inclusive)
# ---------------------------------------------------------------------------
def format_intervals(intervals: Iterable[Sequence[int]]) -> str:
    """Serialise 1-based inclusive intervals as ``"start:end,start:end"``."""
    return _INTERVAL_SEP.join(f"{int(a)}{_BOUND_SEP}{int(b)}" for a, b in intervals)


def parse_intervals(value, ali_from, ali_to) -> Intervals:
    """Parse a ``"start:end,..."`` string into 1-based inclusive intervals.

    Falls back to the single alignment span when the value is absent or empty,
    so a hit without match-state metadata behaves exactly like the legacy
    ali-span rule instead of silently disappearing.
    """
    text = "" if value is None else str(value).strip()
    if not text or text.lower() == "nan":
        return ((int(ali_from), int(ali_to)),)
    intervals: List[Interval] = []
    for part in text.split(_INTERVAL_SEP):
        if _BOUND_SEP not in part:
            continue
        a, b = part.split(_BOUND_SEP, 1)
        try:
            intervals.append((int(float(a)), int(float(b))))
        except (TypeError, ValueError):
            continue
    if not intervals:
        return ((int(ali_from), int(ali_to)),)
    return tuple(intervals)


def intervals_length(intervals: Iterable[Sequence[int]]) -> int:
    """Total residues covered by 1-based inclusive intervals."""
    return sum(int(b) - int(a) + 1 for a, b in intervals)


def intervals_overlap(a_intervals: Sequence[Sequence[int]], b_intervals: Sequence[Sequence[int]]) -> int:
    """Residues shared by two sets of 1-based inclusive intervals.

    Both inputs are treated as disjoint interval sets (which is how
    ``_matched_intervals`` builds them).  A linear merge is used instead of
    per-residue sets so long domains stay cheap.
    """
    if not a_intervals or not b_intervals:
        return 0
    a_sorted = sorted((int(a), int(b)) for a, b in a_intervals)
    b_sorted = sorted((int(a), int(b)) for a, b in b_intervals)
    total = 0
    i = j = 0
    while i < len(a_sorted) and j < len(b_sorted):
        a0, a1 = a_sorted[i]
        b0, b1 = b_sorted[j]
        total += max(0, min(a1, b1) - max(a0, b0) + 1)
        if a1 <= b1:
            i += 1
        else:
            j += 1
    return total


def base_hmm_family(value) -> str:
    """``"GH5_39.hmm"`` -> ``"GH5"``; coexistence policies are keyed on this."""
    name = str(value)
    if name.endswith(".hmm"):
        name = name[:-4]
    return name.split("_", 1)[0]


def _as_float(value, default: float) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return default
    if math.isnan(result):
        return default
    return result


# ---------------------------------------------------------------------------
# Domain hit
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class DomainHit:
    """One HMM domain alignment, with geometry and match-state evidence apart.

    ``index`` is the stable identity of the row this hit came from.  It is used
    to restore input order in the output and as the last resort tie-break, so
    resolution never depends on the order hits are handed to the resolver.
    """

    hmm_name: str
    target_name: str
    source: str = ""

    hmm_from: int = 0
    hmm_to: int = 0

    ali_from: int = 0
    ali_to: int = 0

    matched_intervals: Intervals = field(default=())

    i_evalue: float = math.inf
    bitscore: float = 0.0
    coverage: float = 0.0

    index: int = 0

    @property
    def effective_intervals(self) -> Intervals:
        """Match-state intervals, or the alignment span if none are known."""
        return self.matched_intervals or ((self.ali_from, self.ali_to),)

    @property
    def matched_length(self) -> int:
        return intervals_length(self.effective_intervals)

    @property
    def base_family(self) -> str:
        return base_hmm_family(self.hmm_name)

    @property
    def evidence_rank(self) -> tuple:
        """Deterministic "better first" ordering.

        i-Evalue is the historical evidence measure and stays primary.  The
        remaining terms exist only to make ties total and therefore
        order-independent; they reproduce the legacy tie-break, which kept the
        earlier-starting hit when i-Evalues were equal.
        """
        return (
            self.i_evalue,
            self.ali_from,
            self.ali_to,
            str(self.hmm_name),
            self.index,
        )


# ---------------------------------------------------------------------------
# Conflict determination
# ---------------------------------------------------------------------------
def geometric_overlap(a: DomainHit, b: DomainHit, use_matched: bool) -> Tuple[int, int, int]:
    """Return ``(overlap, length_a, length_b)`` under the requested semantics."""
    if use_matched:
        ia = a.effective_intervals
        ib = b.effective_intervals
        return (
            intervals_overlap(ia, ib),
            max(1, intervals_length(ia)),
            max(1, intervals_length(ib)),
        )
    # Legacy alignment-span arithmetic, preserved verbatim (including its
    # historical off-by-one: neither the overlap nor the lengths add the +1 an
    # inclusive span would need).  Legacy mode exists for reproducibility, so
    # this must not be "fixed" here.
    overlap = min(a.ali_to, b.ali_to) - max(a.ali_from, b.ali_from)
    return (overlap, max(1, a.ali_to - a.ali_from), max(1, b.ali_to - b.ali_from))


def overlaps_substantially(
    a: DomainHit,
    b: DomainHit,
    overlap_ratio_threshold: float = DEFAULT_OVERLAP_RATIO_THRESHOLD,
    use_matched: bool = False,
) -> bool:
    """Pure geometry: does either hit lose more than the threshold fraction?

    The comparison is strict (``>``), so an overlap of exactly the threshold
    fraction is not substantial.  This is the historical dbCAN behaviour.
    """
    overlap, len_a, len_b = geometric_overlap(a, b, use_matched)
    if overlap <= 0:
        return False
    return (overlap / len_a) > overlap_ratio_threshold or (overlap / len_b) > overlap_ratio_threshold


def may_coexist(a: DomainHit, b: DomainHit, coexistence_pairs: Optional[Set[frozenset]]) -> bool:
    """Policy: is this unordered base-family pair explicitly allowed to overlap?"""
    if not coexistence_pairs:
        return False
    return frozenset((a.base_family, b.base_family)) in coexistence_pairs


def competes(
    a: DomainHit,
    b: DomainHit,
    overlap_ratio_threshold: float = DEFAULT_OVERLAP_RATIO_THRESHOLD,
    coexistence_pairs: Optional[Set[frozenset]] = None,
    use_matched: bool = False,
    cross_source_conflicts: bool = False,
) -> bool:
    """True competition = substantial overlap AND no policy allowing both.

    Hits produced by different annotation sources do not compete by default.
    CGC signature classes (CAZyme, TC, TF, STP, sulfatase, peptidase) are not
    mutually exclusive, so a generic family-competition rule must not be
    applied across them.
    """
    if not cross_source_conflicts and a.source != b.source:
        return False
    if not overlaps_substantially(a, b, overlap_ratio_threshold, use_matched):
        return False
    return not may_coexist(a, b, coexistence_pairs)


# ---------------------------------------------------------------------------
# Resolution
# ---------------------------------------------------------------------------
def build_conflict_graph(
    hits: Sequence[DomainHit],
    overlap_ratio_threshold: float = DEFAULT_OVERLAP_RATIO_THRESHOLD,
    coexistence_pairs: Optional[Set[frozenset]] = None,
    use_matched: bool = False,
    cross_source_conflicts: bool = False,
) -> List[Set[int]]:
    """Adjacency list over hit positions; an edge means "true competitors"."""
    n = len(hits)
    adjacency: List[Set[int]] = [set() for _ in range(n)]
    if n < 2:
        return adjacency

    # Match-state intervals are always contained in the alignment span, so
    # disjoint ali spans can never overlap in match-state terms either.  A
    # sweep over ali_from therefore only tests genuinely nearby pairs instead
    # of all n^2 of them.
    order = sorted(range(n), key=lambda k: (hits[k].ali_from, hits[k].ali_to, k))
    for pos, i in enumerate(order):
        for j in order[pos + 1:]:
            if hits[j].ali_from > hits[i].ali_to:
                break
            if competes(
                hits[i],
                hits[j],
                overlap_ratio_threshold=overlap_ratio_threshold,
                coexistence_pairs=coexistence_pairs,
                use_matched=use_matched,
                cross_source_conflicts=cross_source_conflicts,
            ):
                adjacency[i].add(j)
                adjacency[j].add(i)
    return adjacency


def _connected_components(adjacency: Sequence[Set[int]]) -> List[List[int]]:
    seen = [False] * len(adjacency)
    components: List[List[int]] = []
    for start in range(len(adjacency)):
        if seen[start]:
            continue
        stack = [start]
        seen[start] = True
        component = []
        while stack:
            node = stack.pop()
            component.append(node)
            for neighbour in adjacency[node]:
                if not seen[neighbour]:
                    seen[neighbour] = True
                    stack.append(neighbour)
        components.append(component)
    return components


def resolve_domain_hits(
    hits: Sequence[DomainHit],
    overlap_ratio_threshold: float = DEFAULT_OVERLAP_RATIO_THRESHOLD,
    coexistence_pairs: Optional[Set[frozenset]] = None,
    use_matched: bool = False,
    cross_source_conflicts: bool = False,
) -> List[DomainHit]:
    """Resolve one target protein's domain architecture.

    Each hit is a node; an edge marks two hits as true competitors.  Every
    connected component is then resolved by walking its members best-evidence
    first and keeping a hit whenever it does not conflict with something
    already kept.  That keeps a non-conflicting third hit that the old
    ``keep[-1]`` chain could drop or spare purely by row order, and the result
    depends only on the evidence ranking, never on input order.

    Returned hits are ordered by their original ``index``.
    """
    hits = list(hits)
    if len(hits) < 2:
        return hits

    adjacency = build_conflict_graph(
        hits,
        overlap_ratio_threshold=overlap_ratio_threshold,
        coexistence_pairs=coexistence_pairs,
        use_matched=use_matched,
        cross_source_conflicts=cross_source_conflicts,
    )

    kept: List[int] = []
    for component in _connected_components(adjacency):
        if len(component) == 1:
            kept.append(component[0])
            continue
        accepted: List[int] = []
        for node in sorted(component, key=lambda k: hits[k].evidence_rank):
            if any(other in adjacency[node] for other in accepted):
                continue
            accepted.append(node)
        kept.extend(accepted)

    kept.sort(key=lambda k: hits[k].index)
    return [hits[k] for k in kept]


def hit_from_mapping(row, index: int = 0, source: str = "", columns=None) -> DomainHit:
    """Build a :class:`DomainHit` from a row-like mapping (e.g. a pandas Series).

    ``columns`` maps logical field names to the row's column names; missing
    columns fall back to neutral values so partially-populated tables (the CGC
    signature table has no ``HMM Name``, for instance) still resolve.
    """
    columns = columns or {}

    def _get(field_name: str, default):
        key = columns.get(field_name, field_name)
        try:
            if key not in row:
                return default
        except TypeError:
            return default
        value = row[key]
        return default if value is None else value

    ali_from = int(_as_float(_get("ali_from", 0), 0.0))
    ali_to = int(_as_float(_get("ali_to", 0), 0.0))
    return DomainHit(
        hmm_name=str(_get("hmm_name", "")),
        target_name=str(_get("target_name", "")),
        source=source,
        hmm_from=int(_as_float(_get("hmm_from", 0), 0.0)),
        hmm_to=int(_as_float(_get("hmm_to", 0), 0.0)),
        ali_from=ali_from,
        ali_to=ali_to,
        matched_intervals=parse_intervals(_get("matched_intervals", None), ali_from, ali_to),
        i_evalue=_as_float(_get("i_evalue", math.inf), math.inf),
        bitscore=_as_float(_get("bitscore", 0.0), 0.0),
        coverage=_as_float(_get("coverage", 0.0), 0.0),
        index=index,
    )
