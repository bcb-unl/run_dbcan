"""Which CGC signature labels of a protein describe the same region as its CAZyme domain.

A domain-level signature (STP or TF HMM domain, MEROPS peptidase unit) that lies on the
CAZyme domain of the same protein is usually one domain classified by two databases, e.g.
CE1 (dbCAN) and S09 (MEROPS) on an alpha/beta-hydrolase, or CE9 and M38 on an
amidohydrolase. Such a protein must count once, as a CAZyme, when CGCs are called.

Only domain-level sources are compared by coordinates. TC, TF (DIAMOND) and sulfatase hits
come from full-length protein databases: their query interval is an alignment span, not a
domain, so they are never treated as the same region as a CAZyme domain. E-values are not
used (DIAMOND and HMMER scores are not comparable).
"""
from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import Dict, List, Set, Tuple

import dbcan.constants.base_constants as B
import dbcan.constants.process_utils_constants as PU
from dbcan.IO.table import read_tsv

logger = logging.getLogger(__name__)

# overview.tsv domain strings: "GH5_4(20-300)+CBM2(320-400)"
_DOMAIN_RE = re.compile(r"\((\d+)-(\d+)\)")


def _overlap_fraction(a: Tuple[int, int], b: Tuple[int, int]) -> float:
    """Overlap of two inclusive intervals relative to the shorter one."""
    ov = min(a[1], b[1]) - max(a[0], b[0]) + 1
    if ov <= 0:
        return 0.0
    return ov / min(a[1] - a[0] + 1, b[1] - b[0] + 1)


def cazyme_domain_intervals(overview_path: Path) -> Dict[str, List[Tuple[int, int]]]:
    """CAZyme domain coordinates (dbCAN HMM and dbCAN-sub) of proteins used as CGC CAZymes."""
    df = read_tsv(overview_path)
    if "#ofTools" not in df.columns:
        return {}
    df = df[df["#ofTools"].astype(int) >= B.GFF_MIN_TOOL_COUNT]
    out: Dict[str, List[Tuple[int, int]]] = {}
    for _, row in df.iterrows():
        spans = []
        for col in ("dbCAN_hmm", "dbCAN_sub"):
            spans += [(int(a), int(b)) for a, b in _DOMAIN_RE.findall(str(row.get(col, "")))]
        if spans:
            out[str(row["Gene ID"])] = spans
    return out


def same_region_categories(output_dir, threshold: float = PU.SAME_REGION_OVERLAP_THRESHOLD) -> Dict[str, Set[str]]:
    """protein -> signature categories whose every hit lies on a CAZyme domain of that protein.

    A category is returned only if all of its hits on the protein are domain-level and
    overlap a CAZyme domain by more than ``threshold`` of the shorter interval; a single
    independent hit (another STP domain elsewhere, or any full-length DIAMOND label) keeps
    the category as a real, separate signature.
    """
    out_dir = Path(output_dir)
    overview, sig = out_dir / B.OVERVIEW_FILE, out_dir / B.GFF_CGC_SIG_FILE
    if not overview.exists() or not sig.exists():
        logger.warning("CGC same-region check skipped: %s or %s not found", overview.name, sig.name)
        return {}
    cazy = cazyme_domain_intervals(overview)
    df = read_tsv(sig)
    if "Source" not in df.columns:
        logger.warning("CGC same-region check skipped: %s has no Source column (re-run gff_process)", sig.name)
        return {}
    df = df[df[PU.TARGET_NAME_COLUMN].isin(cazy.keys())]
    result: Dict[str, Set[str]] = {}
    for (prot, cat), rows in df.groupby([PU.TARGET_NAME_COLUMN, "Type"], sort=False):
        independent = False
        for _, r in rows.iterrows():
            if r["Source"] not in PU.DOMAIN_LEVEL_SIGNATURE_SOURCES:
                independent = True
                break
            span = (int(r[PU.TARGET_FROM_COLUMN]), int(r[PU.TARGET_TO_COLUMN]))
            if not any(_overlap_fraction(span, c) > threshold for c in cazy[prot]):
                independent = True
                break
        if not independent:
            result.setdefault(str(prot), set()).add(str(cat))
    if result:
        logger.info("CGC: %d proteins have signature labels on their CAZyme domain (counted as CAZyme only)", len(result))
    return result
