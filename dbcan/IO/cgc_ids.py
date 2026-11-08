"""CGC identifiers and CGC_annotation roles, shared by substrate prediction, abundance and plots.

* ``contig|CGC12`` is a CGC ID; the contig may itself contain ``|`` (``gi|123|ref|NC_1.1|``), so
  IDs are parsed with patterns anchored on the ``CGC<n>`` token, never split by position.
* ``<contig>|<CGC#>|<protein>|<type>|<family>`` is the ID of a CGC protein in ``CGC.faa`` and
  the query ID of ``PUL_blast.out`` (:func:`cgc_query_id` builds it from a cgc_standard_out row,
  :func:`parse_query_id` reads it back). Type and
  family never contain ``|`` (see :func:`role_family`), so the protein may.
* A ``CGC_annotation`` / ``Gene Annotation`` value lists one or more roles,
  ``CAZyme|GH13_22|GT5_e25+TC|8.A.9.1.1``; :func:`parse_roles` returns them ordered by the
  CGCFinder Gene Type priority CAZyme > TC > TF > STP > Sulfatase > Peptidase.
"""
from __future__ import annotations

import re
from typing import Dict, List, NamedTuple, Optional, Tuple

import dbcan.constants.base_constants as B

OTHER = B.CGC_OTHER_GENE_TYPE

# 'contig|CGC12' or 'contig|CGC12|protein|type|family'
_CGC_ID_RE = re.compile(r"^(?P<contig>.+?)\|(?P<cgc>CGC\d+)(?:\|(?P<rest>.*))?$")
_QUERY_ID_RE = re.compile(
    r"^(?P<contig>.+?)\|(?P<cgc>CGC\d+)\|(?P<protein>.+)\|(?P<type>[^|]+)\|(?P<family>[^|]*)$")
# separators inside one role's detail: '|' (CAZyme families) and ';' (merged roles of one type)
_DETAIL_SEP = re.compile(r"[|;]")


def split_cgc_id(text: str) -> Optional[Tuple[str, str, str]]:
    """'NC_1|CGC3|NP_1.1|TC|2.A' -> ('NC_1', 'CGC3', 'NP_1.1|TC|2.A'); None if no CGC token."""
    m = _CGC_ID_RE.match(text or "")
    if not m:
        return None
    return m.group("contig"), m.group("cgc"), m.group("rest") or ""


def split_cgc_key(cgcid: str) -> Optional[Tuple[str, str]]:
    """'contig|CGC12' (a CGC ID, nothing after the CGC token) -> ('contig', 'CGC12').

    The CGC token is the last field, so this also works when the contig itself contains a
    '|CGC<n>' token ('scaffold|CGC7|v2|CGC1' -> ('scaffold|CGC7|v2', 'CGC1')).
    """
    contig, sep, cgc = (cgcid or "").rpartition("|")
    if not sep or not contig or not re.fullmatch(r"CGC\d+", cgc):
        return None
    return contig, cgc


class QueryID(NamedTuple):
    """Parts of a CGC protein's DIAMOND query ID, ``<contig>|<CGC#>|<protein>|<type>|<family>``."""
    contig: str
    cgc: str
    protein: str
    type: str
    family: str

    @property
    def cgc_id(self) -> str:
        return f"{self.contig}|{self.cgc}"

    def format(self) -> str:
        return format_query_id(self.contig, self.cgc, self.protein, self.type, self.family)


def format_query_id(contig: str, cgc: str, protein: str, gene_type: str, family: str) -> str:
    return f"{contig}|{cgc}|{protein}|{gene_type}|{family}"


def cgc_query_id(contig: str, cgc: str, protein: str, gene_type: str, annotation: str) -> QueryID:
    """Query ID of a cgc_standard_out.tsv row (as written to CGC.faa by substrate_prediction).

    The family is that of the Gene Type role (the CGCFinder priority role) of the annotation,
    TC IDs cut to the subfamily level; 'null' for genes without a signature.
    """
    if (gene_type or "").lower() == B.NULL.lower():
        family = B.NULL
    else:
        family = role_family(annotation, gene_type,
                             tc_levels=B.CGC_QUERY_TC_LEVELS if gene_type == B.TC else None)
    return QueryID(contig, cgc, protein, gene_type, family)


def parse_query_id(qseqid: str, known: Optional[Dict[str, QueryID]] = None) -> Optional[QueryID]:
    """Parse a CGC protein query ID (CGC.faa header, PUL_blast.out column 1).

    ``known`` (query ID -> QueryID, built from cgc_standard_out.tsv) is used first; the pattern
    only for IDs that are not in it. Contig and protein may contain '|'.
    """
    if known and qseqid in known:
        return known[qseqid]
    m = _QUERY_ID_RE.match(qseqid or "")
    if not m:
        return None
    return QueryID(m["contig"], m["cgc"], m["protein"], m["type"], m["family"])


def normalize_gene_type(name: str) -> str:
    """Map a signature type spelling (SULFATLAS, PEPTIDASE, cazyme ...) to its canonical name."""
    return B.CGC_GENE_TYPE_ALIASES.get((name or "").strip().lower(), OTHER)


def parse_roles(annotation: str) -> List[Tuple[str, str]]:
    """'CAZyme|GH13_31+TC|8.A.9.1.1' -> [('CAZyme', 'GH13_31'), ('TC', '8.A.9.1.1')].

    Roles of the same type are merged (details joined by ';'); roles are ordered by the
    Gene Type priority; 'null' and unknown types (Other|tRNA) are dropped.
    """
    merged: Dict[str, List[str]] = {}
    for part in (annotation or "").split("+"):
        part = part.strip()
        if not part or part.lower() in ("null", "na", "-"):
            continue
        name, _, detail = part.partition("|")
        gtype = normalize_gene_type(name)
        if gtype == OTHER:
            continue
        details = merged.setdefault(gtype, [])
        detail = detail.strip()
        if detail and detail not in details:
            details.append(detail)
    order = B.CGC_GENE_TYPE_PRIORITY
    return [(t, ";".join(d)) for t, d in sorted(merged.items(), key=lambda kv: order.get(kv[0], 99))]


def role_family(annotation: str, gene_type: str = "", tc_levels: Optional[int] = None) -> str:
    """Family/ID of the role ``gene_type`` (default: the highest-priority role) of an annotation.

    'CAZyme|GH32_e132+Sulfatase|S1_7' -> 'GH32_e132'; 'STP|HTH_1+STP|LysR_substrate' ->
    'HTH_1+LysR_substrate'; 'CAZyme|GH13_22|GT5_e25' -> 'GH13_22+GT5_e25'; with tc_levels=3,
    'TC|1.A.1.1.3' -> '1.A.1'. Several families of one role are joined by '+', so the result
    never contains '|' or ';'. Annotations without a signature role ('Other|tRNA') give their
    text after the first '|'; 'null' gives 'null'.
    """
    annotation = (annotation or "").strip()
    roles = parse_roles(annotation)
    if not roles:
        if not annotation or annotation.lower() == B.NULL.lower():
            return B.NULL
        detail = annotation.partition("|")[2] or annotation
        return "+".join(p for p in _DETAIL_SEP.split(detail) if p) or B.NULL
    wanted = normalize_gene_type(gene_type) if gene_type else roles[0][0]
    gtype, detail = next(((t, d) for t, d in roles if t == wanted), roles[0])
    parts: List[str] = []
    for p in _DETAIL_SEP.split(detail):
        p = p.strip()
        if gtype == B.TC and tc_levels:
            p = ".".join(p.split(".")[:tc_levels])
        if p and p not in parts:
            parts.append(p)
    return "+".join(parts) or gtype
