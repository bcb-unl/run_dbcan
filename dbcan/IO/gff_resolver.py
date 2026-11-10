"""Resolve the genes of an annotation GFF3/GTF to the proteins of the protein FASTA.

``cgc.gff`` needs one row per gene with the ID of that gene's protein in ``uniInput.faa``.
Eukaryotic annotation pipelines (funannotate/funannotate2, BRAKER/AUGUSTUS, Helixer, Ensembl,
NCBI GTF, ...) describe genes as ``gene -> transcript (mRNA) -> CDS`` and name the protein after
the transcript, the CDS or the gene, never through one fixed attribute. This module:

* reads GFF3 (``key=value``, URL-decoded, several ``Parent`` values) and GTF
  (``key "value";``) line by line, without loading exon/UTR rows;
* groups CDS rows into transcripts (any feature that is the ``Parent`` of a CDS, or the GTF
  ``transcript_id``) and transcripts into genes;
* tries, per transcript, the identifiers that pipelines use for the protein (CDS
  ``protein_id``/``proteinId``/``Name``/``ID`` without ``cds-``/``CDS:`` prefixes and ``.cds``/
  ``.CDS.n`` suffixes, transcript ``ID``/``transcript_id``, gene ``ID``/``Name``/``locus_tag``),
  each passed through the same ID cleaning that ``--mode protein`` applies to FASTA headers,
  and keeps the first one found in the FASTA (version-less accessions and JGI ``proteinId``
  values only when they are unique);
* keeps one protein per gene: the isoform annotated as a CAZyme / CGC signature, else the
  longest CDS;
* keeps genes without a CDS (tRNA, rRNA, ...) as ``protein_id=NA`` rows, as ``NCBI_euk`` does,
  so that the null-gene distances used by ``cgc_finder`` stay the same.

It also reports, for every ``--gff_type`` (including the original prodigal / NCBI / JGI
readers, whose ``cgc.gff`` code is unchanged), how many protein-coding genes were matched, so
that a wrong ``--gff_type`` stops the run with a clear message instead of giving an empty CGC
table (see :func:`check_gff_resolution`).
"""
from __future__ import annotations

import gzip
import logging
import re
from collections import Counter, OrderedDict, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Set, Tuple
from urllib.parse import unquote

import dbcan.constants.gff_constants as G
from dbcan.core import DbcanStepError
from dbcan.IO.fasta import clean_sequence_id

logger = logging.getLogger(__name__)

_GTF_ATTR = re.compile(r'\s*([^\s"=;]+)\s+"([^"]*)"\s*;?')
_GTF_HINT = re.compile(r'^\s*[A-Za-z_][\w.]*\s+"')
_SUFFIXES = tuple(re.compile(p) for p in G.GFF_RESOLVER_ID_SUFFIX_PATTERNS)
_GNL = re.compile(r"^gnl\|[^|]+\|")
_VERSION = re.compile(r"^(.*)\.\d+$")
_CDS_TYPES = frozenset(G.GFF_RESOLVER_CDS_TYPES)
_GENE_TYPES = frozenset(G.GFF_RESOLVER_GENE_TYPES)
_TRANSCRIPT_TYPES = frozenset(G.GFF_RESOLVER_TRANSCRIPT_TYPES)
_IGNORED_TYPES = frozenset(G.GFF_RESOLVER_IGNORED_TYPES)
_PSEUDO_TYPES = frozenset(G.GFF_RESOLVER_PSEUDO_TYPES)
_STRONG_ATTRS = frozenset(G.GFF_RESOLVER_STRONG_KEY_ATTRS)

# Where each --gff_type looks for the protein ID first. Every resolver type falls back to the
# generic GFF3 order, so a key used by another pipeline is still found; the order only matters
# when two keys of one transcript name two different proteins of the FASTA.
_GENERIC_KEYS: Tuple[Tuple[str, Tuple[str, ...]], ...] = (
    ("CDS", ("protein_id", "proteinId", "Name", "ID", "transcript_id", "orig_protein_id")),
    ("transcript", ("protein_id", "proteinId", "ID", "Name", "transcript_id", "orig_protein_id",
                    "orig_transcript_id", "locus_tag")),
    ("gene", ("proteinId", "ID", "Name", "locus_tag", "gene_id")),
    ("CDS", ("gene_id", "locus_tag")),
)
_KEY_PROFILES: Dict[str, Tuple[Tuple[str, Tuple[str, ...]], ...]] = {
    G.GFF_FORMAT_GFF3: _GENERIC_KEYS,
    # funannotate / funannotate2: proteins are named by the mRNA ID (FUN_000001-T1); CDS rows
    # carry ID=FUN_000001-T1.cds; after `funannotate annotate` the locus tag may replace FUN_.
    G.GFF_FORMAT_FUNANNOTATE: (
        ("transcript", ("ID",)),
        ("CDS", ("ID", "protein_id", "Name")),
        ("transcript", ("Name", "locus_tag")),
        ("gene", ("ID", "Name", "locus_tag")),
    ) + _GENERIC_KEYS,
    # BRAKER / AUGUSTUS GFF3: gene g1 -> mRNA or transcript g1.t1 -> CDS g1.t1.CDS1 / g1.t1.cds;
    # braker.aa and augustus.aa are named by the transcript ID (g1.t1).
    G.GFF_FORMAT_BRAKER: (
        ("transcript", ("ID", "transcript_id")),
        ("CDS", ("ID",)),
        ("gene", ("ID", "gene_id")),
    ) + _GENERIC_KEYS,
    # GTF: CDS rows carry transcript_id / gene_id (and protein_id in NCBI/Ensembl GTF).
    G.GFF_FORMAT_GTF: (
        ("CDS", ("protein_id",)),
        ("transcript", ("transcript_id",)),
        ("CDS", ("transcript_id",)),
        ("gene", ("gene_id",)),
        ("CDS", ("gene_id",)),
    ) + _GENERIC_KEYS,
}


class GffResolutionError(DbcanStepError):
    """Too few protein-coding genes of the GFF could be matched to a protein of the FASTA."""


# --------------------------------------------------------------------------------------------
# --gff_type helpers
# --------------------------------------------------------------------------------------------
def normalize_gff_type(gff_type: Optional[str]) -> Optional[str]:
    """Return the canonical spelling of a --gff_type value (case-insensitive), or None."""
    if gff_type is None:
        return None
    wanted = str(gff_type).strip().lower()
    for choice in G.GFF_FORMAT_CHOICES:
        if choice.lower() == wanted:
            return choice
    raise ValueError(
        f"Unsupported GFF type: {gff_type!r}. Choose one of: {', '.join(G.GFF_FORMAT_CHOICES)}"
    )


def is_resolver_type(gff_type: str) -> bool:
    return normalize_gff_type(gff_type) in G.GFF_RESOLVER_FORMATS


# --------------------------------------------------------------------------------------------
# Parsing
# --------------------------------------------------------------------------------------------
@dataclass
class Feature:
    seqid: str
    ftype: str
    start: int
    end: int
    strand: str
    attrs: Dict[str, str]
    order: int
    gtf: bool = False

    @property
    def fid(self) -> Optional[str]:
        return self.attrs.get("ID")

    @property
    def parents(self) -> List[str]:
        return [p for p in self.attrs.get("Parent", "").split(",") if p]


def parse_attributes(col9: str) -> Tuple[Dict[str, str], bool]:
    """Parse GFF3 or GTF column 9. Returns (attributes, is_gtf).

    AUGUSTUS GTF writes a bare ID ('g1', 'g1.t1') on its gene/transcript rows; it is returned
    as ``ID`` with is_gtf=True.
    """
    col9 = col9.strip()
    out: Dict[str, str] = {}
    if not col9 or col9 == ".":
        return out, False
    if _GTF_HINT.match(col9) and '="' not in col9:
        for m in _GTF_ATTR.finditer(col9):
            out.setdefault(m.group(1), m.group(2))
        return out, True
    if "=" not in col9:
        out["ID"] = col9.rstrip(";").strip()
        return out, True
    for part in col9.split(";"):
        part = part.strip()
        if not part or "=" not in part:
            continue
        key, _, val = part.partition("=")
        key = key.strip()
        if key and key not in out:
            out[key] = unquote(val.strip())
    return out, False


def _open_text(path: Path):
    """Open a GFF/GTF for reading; gzip (``.gz`` suffix or gzip magic bytes) is read transparently."""
    path = Path(path)
    is_gz = path.suffix == ".gz"
    if not is_gz:
        with path.open("rb") as fh:
            is_gz = fh.read(2) == b"\x1f\x8b"
    if is_gz:
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open(encoding="utf-8", errors="replace")


def iter_features(path: Path, keep_ignored: bool = False) -> Iterator[Feature]:
    """Yield feature rows of a GFF3/GTF file (comments skipped, stops at ##FASTA)."""
    with _open_text(path) as fh:
        for n, line in enumerate(fh):
            if line.startswith("##FASTA") or line.startswith(">"):
                break
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip("\n\r").split("\t")
            if len(p) < 9:
                continue
            if not keep_ignored and p[2] in _IGNORED_TYPES:
                continue
            try:
                s, e = int(p[3]), int(p[4])
            except ValueError:
                continue
            attrs, is_gtf = parse_attributes(p[8])
            yield Feature(p[0], p[2], min(s, e), max(s, e), p[6], attrs, n, is_gtf)


# --------------------------------------------------------------------------------------------
# FASTA side
# --------------------------------------------------------------------------------------------
def read_fasta_ids(path: Path, clean: bool = True) -> List[str]:
    """Protein IDs of a FASTA file (first header word), optionally cleaned like --mode protein."""
    ids: List[str] = []
    with _open_text(path) as fh:
        for line in fh:
            if line.startswith(">"):
                head = line[1:].strip()
                rid = head.split()[0] if head else ""
                ids.append(clean_sequence_id(rid) if clean else rid)
    return ids


class ProteinIndex:
    """Exact protein IDs plus unambiguous aliases (version-less accession, JGI proteinId, '|' parts)."""

    def __init__(self, ids: Iterable[str]):
        self.ids: List[str] = list(ids)
        self.exact: Set[str] = set(self.ids)
        self._alias: Optional[Dict[str, str]] = None
        self._ambiguous: Set[str] = set()

    def _build_aliases(self) -> None:
        # built on first use only: the prodigal/NCBI/JGI reports need exact IDs (metagenome scale)
        hits: Dict[str, Set[str]] = defaultdict(set)
        for pid in self.ids:
            for alias in self._aliases(pid):
                if alias != pid:
                    hits[alias].add(pid)
        self._alias = {a: next(iter(s)) for a, s in hits.items() if len(s) == 1 and a not in self.exact}
        self._ambiguous = {a for a, s in hits.items() if len(s) > 1}

    @property
    def alias(self) -> Dict[str, str]:
        if self._alias is None:
            self._build_aliases()
        return self._alias

    @property
    def ambiguous(self) -> Set[str]:
        if self._alias is None:
            self._build_aliases()
        return self._ambiguous

    @staticmethod
    def _aliases(pid: str) -> List[str]:
        out: List[str] = []
        m = _VERSION.match(pid)                       # KZF25850.1 -> KZF25850
        if m:
            out.append(m.group(1))
        for sep in ("|", "-"):                        # jgi|Xylhe1|242238|... (cleaned: jgi-Xylhe1-242238-...)
            if pid.lower().startswith("jgi" + sep):
                toks = pid.split(sep)
                if len(toks) >= 3:
                    out.append(toks[2])
                    out.append(f"{toks[1]}{sep}{toks[2]}")
        if "|" in pid:
            out.extend(t for t in pid.split("|") if t)
        return out

    def lookup(self, key: str, allow_alias: bool = True) -> Optional[str]:
        if not key or any(c.isspace() for c in key):
            return None
        key = clean_sequence_id(key)
        if key in self.exact:
            return key
        return self.alias.get(key) if allow_alias else None


def _id_variants(value: str) -> List[str]:
    """value, then the value without gnl|db| / known prefixes / known suffixes (unique, in order)."""
    out = [value]
    y = _GNL.sub("", value)
    out.append(y)
    for pre in G.GFF_RESOLVER_ID_PREFIXES:
        if y.startswith(pre):
            y = y[len(pre):]
            out.append(y)
            break
    for rx in _SUFFIXES:
        z = rx.sub("", y)
        if z != y:
            out.append(z)
    seen: Set[str] = set()
    res: List[str] = []
    for v in out:
        if v and v not in seen:
            seen.add(v)
            res.append(v)
    return res


# --------------------------------------------------------------------------------------------
# Result types
# --------------------------------------------------------------------------------------------
@dataclass
class GeneResolution:
    gene_id: str
    contig: str
    start: int
    end: int
    strand: str
    protein_id: str              # FASTA ID; 'NA' for non-coding genes; 'unknown' if unmatched
    status: str                  # matched | unmatched | noncoding
    matched_by: str = ""         # e.g. 'transcript.ID', 'CDS.ID~' (~: after prefix/suffix removal)
    feature_type: str = "gene"   # non-coding rows: child feature type (tRNA, rRNA, ...)
    n_transcripts: int = 0
    isoform_proteins: Tuple[str, ...] = ()
    candidate_keys: Tuple[str, ...] = ()
    order: int = 0


@dataclass
class ResolutionResult:
    gff_type: str
    rows: List[GeneResolution]
    fasta_total: int
    ambiguous_aliases: int = 0
    cross_gene_collisions: int = 0
    gtf_rows: int = 0
    gff3_rows: int = 0
    # Streamed results (prodigal type, metagenome scale) keep no rows, only counts:
    # (coding, matched, noncoding, proteins placed), up to 5 unmatched rows, and the report path.
    streamed_counts: Optional[Tuple[int, int, int, int]] = None
    streamed_examples: List[GeneResolution] = field(default_factory=list)
    streamed_tsv: Optional[Path] = None

    def _count(self, status: str) -> int:
        return sum(1 for r in self.rows if r.status == status)

    @property
    def n_coding(self) -> int:
        if self.streamed_counts is not None:
            return self.streamed_counts[0]
        return len(self.rows) - self._count(G.GFF_STATUS_NONCODING)

    @property
    def n_matched(self) -> int:
        if self.streamed_counts is not None:
            return self.streamed_counts[1]
        return self._count(G.GFF_STATUS_MATCHED)

    @property
    def n_noncoding(self) -> int:
        if self.streamed_counts is not None:
            return self.streamed_counts[2]
        return self._count(G.GFF_STATUS_NONCODING)

    @property
    def match_rate(self) -> float:
        n = self.n_coding
        return self.n_matched / n if n else 0.0

    @property
    def fasta_placed(self) -> int:
        if self.streamed_counts is not None:
            return self.streamed_counts[3]
        placed: Set[str] = set()
        for r in self.rows:
            if r.isoform_proteins:
                placed.update(r.isoform_proteins)
            elif r.status == G.GFF_STATUS_MATCHED:
                placed.add(r.protein_id)
        return len(placed)

    def matched_by_counts(self) -> Counter:
        if self.streamed_counts is not None:
            return Counter({"ID": self.streamed_counts[1]}) if self.streamed_counts[1] else Counter()
        return Counter(r.matched_by for r in self.rows if r.status == G.GFF_STATUS_MATCHED)

    def unmatched_examples(self, n: int = 5) -> List[GeneResolution]:
        if self.streamed_counts is not None:
            return self.streamed_examples[:n]
        return [r for r in self.rows if r.status == G.GFF_STATUS_UNMATCHED][:n]

    def summary(self) -> str:
        keys = ", ".join(f"{k}={v}" for k, v in self.matched_by_counts().most_common(4)) or "none"
        multi = sum(1 for r in self.rows if r.n_transcripts > 1)
        return (
            f"GFF resolution (--gff_type {self.gff_type}): {self.n_matched}/{self.n_coding} protein-coding "
            f"genes matched to a protein ({100 * self.match_rate:.1f}%); {self.n_noncoding} non-coding genes "
            f"kept as null rows; {multi} genes with several transcripts; "
            f"{self.fasta_total - self.fasta_placed}/{self.fasta_total} FASTA proteins not placed in the GFF; "
            f"matched by: {keys}"
        )


# --------------------------------------------------------------------------------------------
# Resolver (funannotate, BRAKER, GFF3, GTF)
# --------------------------------------------------------------------------------------------
@dataclass
class _Transcript:
    tid: str
    feat: Optional[Feature]
    gene: Optional[Feature]
    gene_key: str
    cds: List[Feature] = field(default_factory=list)
    # False when ``gene`` is not a gene row (no gene parent: the transcript itself, or a CDS
    # parented directly by a gene-less ID); the gene span is then taken from the transcripts.
    gene_is_row: bool = True

    @property
    def cds_len(self) -> int:
        return sum(c.end - c.start + 1 for c in self.cds)

    @property
    def is_pseudo(self) -> bool:
        """Pseudogene: a 'pseudogene'/'pseudogenic_transcript' row, or CDS rows marked pseudo=true
        that carry no protein ID (NCBI). The NCBI readers skip such genes."""
        for f in (self.gene, self.feat):
            if f is not None and f.ftype in _PSEUDO_TYPES:
                return True
        return bool(self.cds) and all(
            c.attrs.get("pseudo", "").lower() == "true" and not _protein_id_attr(c.attrs) for c in self.cds
        )


def _protein_id_attr(attrs: Dict[str, str]) -> Optional[str]:
    return attrs.get("protein_id") or attrs.get("proteinId")


def _gene_key_fallback(attrs: Dict[str, str]) -> Optional[str]:
    """gffread-style transcripts without a gene parent name their gene in geneID/gene_id/gene."""
    for name in G.GFF_RESOLVER_GENE_KEY_ATTRS:
        value = attrs.get(name)
        if value:
            return value.split(",")[0]
    return None


def _fallback_loci(items: List[Tuple[str, Feature, Optional[str]]]) -> List[Tuple[str, List[Tuple[Feature, Optional[str]]]]]:
    """Gene loci of gffread-style transcripts (no gene rows) named by geneID/gene_id/gene.

    ``items`` are (gene value, transcript row, transcript key or None for a transcript without CDS).
    The value alone does not identify a gene ('gene=' holds a symbol, and symbols repeat, e.g. two
    copies of tnpA): transcripts are grouped only when they share contig, strand and value AND
    overlap (transitively). Returns (value, members) per locus, in file order of the first member."""
    by_key: Dict[Tuple[str, str, str], List[Tuple[Feature, Optional[str]]]] = defaultdict(list)
    for value, feat, tid in items:
        by_key[(feat.seqid, feat.strand, value)].append((feat, tid))
    loci: List[Tuple[str, List[Tuple[Feature, Optional[str]]]]] = []
    for (_, _, value), members in by_key.items():
        members.sort(key=lambda m: (m[0].start, m[0].order))
        cur: List[Tuple[Feature, Optional[str]]] = []
        cur_end = -1
        for m in members:
            if not cur or m[0].start > cur_end:
                cur = []
                loci.append((value, cur))
                cur_end = m[0].end
            else:
                cur_end = max(cur_end, m[0].end)
            cur.append(m)
    loci.sort(key=lambda loc: min(f.order for f, _ in loc[1]))
    return loci


def _candidate_keys(tr: _Transcript, profile) -> List[Tuple[str, str, bool]]:
    """(identifier, provenance, strong). Strong: an ID-like attribute (protein_id, ID,
    transcript_id, ...) of the CDS or of a transcript row, or the transcript key itself; weak:
    Name/locus_tag and gene-level identifiers, which isoforms and neighbours can share."""
    tr_is_gene = tr.feat is not None and tr.feat.ftype in _GENE_TYPES
    sources = {
        "CDS": tr.cds[0].attrs if tr.cds else {},
        "transcript": tr.feat.attrs if tr.feat is not None else {},
        "gene": tr.gene.attrs if tr.gene is not None else {},
    }
    keys: List[Tuple[str, str, bool]] = []
    seen: Set[str] = set()
    for src, names in profile:
        attrs = sources[src]
        for name in names:
            value = attrs.get(name)
            if not value:
                continue
            strong = name in _STRONG_ATTRS and (src == "CDS" or (src == "transcript" and not tr_is_gene))
            for i, var in enumerate(_id_variants(value.split(",")[0])):
                if var not in seen:
                    seen.add(var)
                    keys.append((var, f"{src}.{name}" + ("~" if i else ""), strong))
    for var in _id_variants(tr.tid):
        if var not in seen:
            seen.add(var)
            keys.append((var, "transcript.key", not tr_is_gene))
    return keys


def _isoform_score(tr: _Transcript, pid: Optional[str], annotations: Optional[Dict[str, str]]):
    """Isoform choice: has a protein > annotated as CAZyme > annotated as another CGC signature >
    longest CDS > first in the file. The number of signatures is not counted: a shorter isoform
    can pick up an extra (coverage-based) hit, which must not outweigh the longest CDS."""
    has_cazyme = has_sig = 0
    if pid and annotations:
        ann = annotations.get(pid) or ""
        if ann and ann not in (G.GFF_NULL_ANNOTATION, G.GFF_UNKNOWN_ANNOTATION):
            has_sig = 1
            has_cazyme = int(any(p.startswith(G.GFF_CAZYME_PREFIX) for p in ann.split("+")))
    first = tr.cds[0].order if tr.cds else 0
    return (1 if pid else 0, has_cazyme, has_sig, tr.cds_len, -first)


def _resolve_transcripts(path: Path, gff_type: str, idx: ProteinIndex,
                         annotations: Optional[Dict[str, str]]) -> ResolutionResult:
    profile = _KEY_PROFILES[gff_type]
    by_id: Dict[str, Feature] = {}
    children: Dict[str, List[Feature]] = defaultdict(list)
    cds_rows: List[Feature] = []
    top_level: List[Feature] = []
    gtf_tx: Dict[str, Feature] = {}
    gtf_gene: Dict[str, Feature] = {}
    gtf_tx_of_gene: Dict[str, List[str]] = defaultdict(list)
    n_gtf = n_gff3 = 0

    for f in iter_features(path):
        a = f.attrs
        if f.gtf:
            n_gtf += 1
            if f.ftype in _CDS_TYPES:
                cds_rows.append(f)
                continue
            if f.ftype == "gene":
                gtf_gene.setdefault(a.get("gene_id") or a.get("ID", ""), f)
            elif f.ftype in _TRANSCRIPT_TYPES:
                tid = a.get("transcript_id") or a.get("ID", "")
                if tid not in gtf_tx:
                    gtf_tx[tid] = f
                    gid = a.get("gene_id") or tid.rsplit(".", 1)[0]
                    gtf_tx_of_gene[gid].append(tid)
            continue
        n_gff3 += 1
        if f.ftype in _CDS_TYPES:
            cds_rows.append(f)
            continue
        fid = f.fid
        if fid:
            by_id.setdefault(fid, f)
        parents = f.parents
        for p in parents:
            children[p].append(f)
        if not parents:
            top_level.append(f)

    transcripts: "OrderedDict[str, _Transcript]" = OrderedDict()
    fallback_tx: Dict[str, str] = {}      # gffread-style transcript -> geneID/gene_id/gene value
    for c in cds_rows:
        a = c.attrs
        if c.gtf:
            tid = a.get("transcript_id") or a.get("gene_id") or f"__cds{c.order}"
            gid = a.get("gene_id") or tid
            tr = transcripts.get(tid)
            if tr is None:
                tr = transcripts[tid] = _Transcript(tid, gtf_tx.get(tid), gtf_gene.get(gid), gid)
            tr.cds.append(c)
            continue
        parents = c.parents or [c.fid or f"__cds{c.order}"]
        for tid in parents:
            tr = transcripts.get(tid)
            if tr is None:
                tfeat = by_id.get(tid)
                gparents = tfeat.parents if tfeat is not None else []
                gfeat = by_id.get(gparents[0]) if gparents else None
                if gfeat is not None:
                    tr = _Transcript(tid, tfeat, gfeat, gparents[0])
                elif gparents:
                    # Parent=<gene> without a gene row: group the isoforms by that ID
                    tr = _Transcript(tid, tfeat, tfeat, gparents[0], gene_is_row=False)
                elif tfeat is not None and tfeat.ftype not in _GENE_TYPES and _gene_key_fallback(tfeat.attrs):
                    # gffread-style transcript (no gene rows): geneID=/gene_id=/gene= names the gene;
                    # the gene key is set per locus below
                    tr = _Transcript(tid, tfeat, tfeat, tid, gene_is_row=False)
                    fallback_tx[tid] = _gene_key_fallback(tfeat.attrs)
                else:
                    tr = _Transcript(tid, tfeat, tfeat, tid, gene_is_row=tfeat is not None)
                transcripts[tid] = tr
            tr.cds.append(c)

    # gffread-style loci: coding transcripts and transcripts without CDS that share contig, strand,
    # geneID/gene_id/gene value and overlap form one gene. A locus with a coding transcript is a
    # coding gene (its non-coding isoforms add nothing); otherwise one non-coding null row spanning
    # all its transcripts. The value names the gene when it names one locus only; a repeated value
    # (gene symbol on several loci) names each locus after its first transcript.
    fb_nc_feats = [
        f for f in top_level
        if f.ftype in _TRANSCRIPT_TYPES and f.ftype not in _GENE_TYPES
        and (f.fid or "") not in transcripts and _gene_key_fallback(f.attrs)
    ]
    loci = _fallback_loci(
        [(v, transcripts[t].feat, t) for t, v in fallback_tx.items()]
        + [(_gene_key_fallback(f.attrs), f, None) for f in fb_nc_feats]
    )
    fb_handled: Set[int] = {f.order for f in fb_nc_feats}
    fb_noncoding: List[GeneResolution] = []
    n_loci = Counter(value for value, _ in loci)
    taken: Set[str] = {tr.gene_key for t, tr in transcripts.items() if t not in fallback_tx}
    taken.update(by_id)
    taken.difference_update(t for t in fallback_tx)
    taken.difference_update(f.fid for f in fb_nc_feats if f.fid)
    for value, members in loci:
        first = min(members, key=lambda m: m[0].order)[0]
        key = value if n_loci[value] == 1 and value not in taken else (first.fid or f"__line{first.order}")
        while key in taken:
            key += "_"
        taken.add(key)
        coding = [t for _, t in members if t is not None]
        if coding:
            for t in coding:
                transcripts[t].gene_key = key
            continue
        feats = [f for f, _ in members]
        ftype = first.ftype
        if any(f.ftype in _PSEUDO_TYPES for f in feats):
            ftype = G.GFF_PSEUDOGENE_FEATURE
        fb_noncoding.append(GeneResolution(
            key, first.seqid, min(f.start for f in feats), max(f.end for f in feats), first.strand,
            G.GFF_NA_PROTEIN_ID, G.GFF_STATUS_NONCODING, "", ftype, order=first.order))

    # transcript -> protein, in two passes: first only strong keys found exactly in the FASTA
    # (CDS protein_id, transcript ID, ...), then every key (Name, gene IDs, version-less/JGI
    # aliases). A weak key of an earlier gene can thus never take the protein that a later gene
    # names exactly. Pseudogenes only take a protein named by their CDS protein ID.
    tx_protein: Dict[str, Tuple[str, str]] = {}
    claimed: Dict[str, str] = {}          # protein -> gene key of the transcript that uses it
    cross_gene = 0
    cand = {tid: _candidate_keys(tr, profile) for tid, tr in transcripts.items()}
    pseudo = {tid for tid, tr in transcripts.items() if tr.is_pseudo}
    for strong_pass in (True, False):
        for tid, tr in transcripts.items():
            if tid in tx_protein:
                continue
            for key, prov, strong in cand[tid]:
                if strong_pass and not strong:
                    continue
                if tid in pseudo and not prov.startswith(("CDS.protein_id", "CDS.proteinId")):
                    continue
                pid = idx.lookup(key, allow_alias=not strong_pass)
                if pid is None:
                    continue
                owner = claimed.get(pid)
                if owner is not None:
                    # already used by another transcript (e.g. isoforms matched through their shared
                    # gene ID); try the next key of this transcript
                    if owner != tr.gene_key and not strong_pass:
                        cross_gene += 1
                    continue
                claimed[pid] = tr.gene_key
                tx_protein[tid] = (pid, prov)
                break

    genes: "OrderedDict[str, List[str]]" = OrderedDict()
    for tid, tr in transcripts.items():
        genes.setdefault(tr.gene_key, []).append(tid)

    rows: List[GeneResolution] = []
    coding_keys: Set[str] = set()
    for gkey, tids in genes.items():
        best = max(tids, key=lambda t: _isoform_score(transcripts[t], tx_protein.get(t, (None,))[0], annotations))
        tr = transcripts[best]
        span = tr.gene if tr.gene is not None and tr.gene_is_row else None
        if span is not None:
            contig, start, end, strand, order = span.seqid, span.start, span.end, span.strand, span.order
        else:
            # no gene row: span of the transcript rows (or of the CDS rows) of every isoform
            parts = [f for t in tids for f in ([transcripts[t].feat] if transcripts[t].feat is not None
                                              else transcripts[t].cds)]
            first = min(parts, key=lambda f: f.order)
            contig, strand, order = first.seqid, first.strand, first.order
            start = min(f.start for f in parts)
            end = max(f.end for f in parts)
        pid, prov = tx_protein.get(best, (G.GFF_UNKNOWN_ANNOTATION, ""))
        status = G.GFF_STATUS_MATCHED if best in tx_protein else G.GFF_STATUS_UNMATCHED
        if status == G.GFF_STATUS_UNMATCHED and all(t in pseudo for t in tids):
            # pseudogene without a protein: a non-coding null row (Other|pseudogene), like tRNAs
            rows.append(GeneResolution(gkey, contig, start, end, strand, G.GFF_NA_PROTEIN_ID,
                                       G.GFF_STATUS_NONCODING, "", G.GFF_PSEUDOGENE_FEATURE, order=order))
        else:
            cands = () if status == G.GFF_STATUS_MATCHED else tuple(k for k, _, _ in cand[best])[:6]
            rows.append(GeneResolution(
                gkey, contig, start, end, strand, pid, status, prov, "gene", len(tids),
                tuple(tx_protein[t][0] for t in tids if t in tx_protein), cands, order,
            ))
        coding_keys.add(gkey)
        if tr.gene is not None and tr.gene.fid:
            coding_keys.add(tr.gene.fid)

    # non-coding genes (tRNA, rRNA, ncRNA, pseudogenes without CDS): null rows
    rows.extend(fb_noncoding)
    for f in top_level:
        if f.ftype not in _GENE_TYPES and f.ftype not in _TRANSCRIPT_TYPES:
            continue
        if f.order in fb_handled or (f.fid or "") in fallback_tx:
            continue                      # gffread-style transcript: handled per locus above
        gid = f.fid or f"__line{f.order}"
        if gid in coding_keys or gid in genes:
            continue
        kids = children.get(f.fid or "", [])
        ftype = kids[0].ftype if kids else f.ftype
        if ftype in _PSEUDO_TYPES or f.ftype in _PSEUDO_TYPES:
            ftype = G.GFF_PSEUDOGENE_FEATURE
        rows.append(GeneResolution(gid, f.seqid, f.start, f.end, f.strand, G.GFF_NA_PROTEIN_ID,
                                   G.GFF_STATUS_NONCODING, "", ftype, order=f.order))
    if gtf_gene or gtf_tx:
        for gid, txs in list(gtf_tx_of_gene.items()) + [(g, []) for g in gtf_gene if g not in gtf_tx_of_gene]:
            if gid in genes or gid in coding_keys or any(t in transcripts for t in txs):
                continue
            gfeat = gtf_gene.get(gid)
            tfeats = [gtf_tx[t] for t in txs]
            src = gfeat or (tfeats[0] if tfeats else None)
            if src is None:
                continue
            start = gfeat.start if gfeat else min(t.start for t in tfeats)
            end = gfeat.end if gfeat else max(t.end for t in tfeats)
            ftype = next((f.attrs[k] for f in tfeats + [gfeat] if f is not None
                          for k in ("transcript_biotype", "gene_biotype", "gene_type") if f.attrs.get(k)),
                         tfeats[0].ftype if tfeats else gfeat.ftype)
            rows.append(GeneResolution(gid, src.seqid, start, end, src.strand, G.GFF_NA_PROTEIN_ID,
                                       G.GFF_STATUS_NONCODING, "", ftype, order=src.order))

    rows.sort(key=lambda r: r.order)
    return ResolutionResult(gff_type, rows, len(idx.ids), len(idx._ambiguous), cross_gene, n_gtf, n_gff3)


# --------------------------------------------------------------------------------------------
# Original readers (prodigal, NCBI_prok, NCBI_euk, JGI): same rule as dbcan.IO.gff, report only
# --------------------------------------------------------------------------------------------
def _resolution_line(r: GeneResolution) -> str:
    return "\t".join([
        r.gene_id, r.contig, str(r.start), str(r.end), r.strand, r.protein_id,
        r.matched_by or ("-" if r.status != G.GFF_STATUS_NONCODING else r.feature_type),
        r.status, str(r.n_transcripts), ",".join(r.isoform_proteins) or "-",
    ]) + "\n"


def _resolve_prodigal(path: Path, idx: ProteinIndex, out_tsv: Optional[Path]) -> ResolutionResult:
    """prodigal: every feature row is a gene, its ID is the protein. Streamed (metagenome-size GFFs):
    no rows are kept; the report is written while reading."""
    coding = matched = 0
    placed: Set[str] = set()
    examples: List[GeneResolution] = []
    fh = None
    if out_tsv is not None:
        out_tsv = Path(out_tsv)
        out_tsv.parent.mkdir(parents=True, exist_ok=True)
        fh = out_tsv.open("w", encoding="utf-8")
        fh.write("\t".join(G.GFF_RESOLUTION_COLUMNS) + "\n")
    try:
        for f in iter_features(path, keep_ignored=True):
            if f.ftype.lower() == "remark" or f.parents:
                continue
            pid = f.attrs.get(G.GFF_ID_ATTR_PRODIGAL) or G.GFF_UNKNOWN_ANNOTATION
            ok = pid in idx.exact
            coding += 1
            if ok:
                matched += 1
                placed.add(pid)
            row = GeneResolution(pid, f.seqid, f.start, f.end, f.strand, pid,
                                 G.GFF_STATUS_MATCHED if ok else G.GFF_STATUS_UNMATCHED,
                                 "ID" if ok else "", order=f.order,
                                 candidate_keys=() if ok else (pid if f.attrs.get("ID") else "<no ID>",))
            if not ok and len(examples) < 5:
                examples.append(row)
            if fh is not None:
                fh.write(_resolution_line(row))
    finally:
        if fh is not None:
            fh.close()
    return ResolutionResult(G.GFF_FORMAT_PRODIGAL, [], len(idx.ids),
                            streamed_counts=(coding, matched, 0, len(placed)),
                            streamed_examples=examples, streamed_tsv=out_tsv)


def _resolve_legacy(path: Path, gff_type: str, idx: ProteinIndex) -> ResolutionResult:
    """NCBI_prok, NCBI_euk, JGI: the rule of the processors in dbcan.IO.gff (report only).
    All rows are kept: those processors look at every direct child of a gene (exon included)."""
    feats = list(iter_features(path, keep_ignored=True))
    rows: List[GeneResolution] = []
    n_gtf = sum(1 for f in feats if f.gtf)
    children: Dict[str, List[Feature]] = defaultdict(list)
    for f in feats:
        for p in f.parents:
            children[p].append(f)
    def jgi_simple(x: str) -> str:   # JGIProcessor: 3rd '-'-separated field (jgi-Xylhe1-242238-...)
        parts = x.split("-")
        return parts[2] if len(parts) >= 3 else x

    jgi_map: Dict[str, str] = {}
    if gff_type == G.GFF_FORMAT_JGI:
        for pid in idx.ids:
            jgi_map.setdefault(jgi_simple(pid), pid)

    for f in feats:
        if f.ftype != G.GFF_GENE_FEATURE or f.parents:
            continue
        gid = f.fid or f"__line{f.order}"
        kids = children.get(gid, [])
        pid, prov = G.GFF_UNKNOWN_ANNOTATION, ""
        noncoding_type = None
        if gff_type == G.GFF_FORMAT_JGI:
            raw = f.attrs.get(G.GFF_JGI_PROTEIN_ID_ATTR, G.GFF_UNKNOWN_ANNOTATION).split(",")[0]
            simple = jgi_simple(raw)
            if simple in jgi_map:
                pid, prov = jgi_map[simple], "gene.proteinId"
            else:
                pid = raw
        elif gff_type == G.GFF_FORMAT_NCBI_PROK:
            bad = next((k for k in kids if G.GFF_CDS_FEATURE not in k.ftype), None)
            if bad is not None:
                noncoding_type = bad.ftype
            else:
                for k in kids:
                    if k.ftype == G.GFF_CDS_FEATURE and k.attrs.get(G.GFF_PROTEIN_ID_ATTR_NCBI):
                        pid, prov = k.attrs[G.GFF_PROTEIN_ID_ATTR_NCBI].split(",")[0], "CDS.protein_id"
                        break
        else:  # NCBI_euk
            bad = next((k for k in kids if k.ftype != G.GFF_MRNA_FEATURE), None)
            if bad is not None:
                noncoding_type = bad.ftype
            else:
                for k in kids:
                    cds = [c for c in children.get(k.fid or "", []) if c.ftype == G.GFF_CDS_FEATURE]
                    if cds:
                        value = cds[0].attrs.get(G.GFF_PROTEIN_ID_ATTR_NCBI)
                        if value:
                            pid, prov = value.split(",")[0], "CDS.protein_id"
                            break
        if noncoding_type is None and not kids and gff_type != G.GFF_FORMAT_JGI:
            # A childless gene row whose biotype is not protein-coding (e.g. E. coli ralA,
            # gene_biotype=other) is a non-coding locus, not a protein missing from the FASTA.
            # The processors still write it with protein_id=unknown; only the report changes.
            biotype = f.attrs.get("gene_biotype") or f.attrs.get("gene_type") or ""
            if biotype and biotype != G.GFF_PROTEIN_CODING_BIOTYPE:
                rows.append(GeneResolution(gid, f.seqid, f.start, f.end, f.strand, G.GFF_UNKNOWN_ANNOTATION,
                                           G.GFF_STATUS_NONCODING, "", biotype, order=f.order))
                continue
        if noncoding_type is not None:
            rows.append(GeneResolution(gid, f.seqid, f.start, f.end, f.strand, G.GFF_NA_PROTEIN_ID,
                                       G.GFF_STATUS_NONCODING, "", noncoding_type, order=f.order))
            continue
        ok = bool(prov) and (gff_type == G.GFF_FORMAT_JGI or pid in idx.exact)
        cands: Tuple[str, ...] = ()
        if not ok:
            # what this reader looked at: only one attribute, never the gene/mRNA IDs
            attr = "gene proteinId" if gff_type == G.GFF_FORMAT_JGI else "CDS protein_id"
            cands = (pid,) if pid != G.GFF_UNKNOWN_ANNOTATION else (f"<no {attr}>",)
        rows.append(GeneResolution(gid, f.seqid, f.start, f.end, f.strand, pid,
                                   G.GFF_STATUS_MATCHED if ok else G.GFF_STATUS_UNMATCHED,
                                   prov if ok else "", order=f.order, n_transcripts=len(kids),
                                   candidate_keys=cands))
    return ResolutionResult(gff_type, rows, len(idx.ids), gtf_rows=n_gtf, gff3_rows=len(feats) - n_gtf)


# --------------------------------------------------------------------------------------------
# Public API
# --------------------------------------------------------------------------------------------
def resolve_gff(gff_path, gff_type: str, protein_ids: Iterable[str],
                annotations: Optional[Dict[str, str]] = None, out_tsv=None) -> ResolutionResult:
    """Match every gene of ``gff_path`` to a protein of ``protein_ids`` (cleaned FASTA IDs).

    ``annotations`` (protein -> CGC_annotation) is used to pick the annotated isoform; it only
    matters for the resolver types (funannotate, BRAKER, GFF3, GTF). ``out_tsv``: write the
    per-gene report (gff_resolution.tsv); the prodigal type writes it while reading.
    """
    gff_type = normalize_gff_type(gff_type)
    idx = ProteinIndex(protein_ids)
    path = Path(gff_path)
    if gff_type == G.GFF_FORMAT_PRODIGAL:
        return _resolve_prodigal(path, idx, out_tsv)
    if gff_type in G.GFF_RESOLVER_FORMATS:
        result = _resolve_transcripts(path, gff_type, idx, annotations)
    else:
        result = _resolve_legacy(path, gff_type, idx)
    if out_tsv is not None:
        write_resolution_tsv(result, out_tsv)
    return result


def write_resolution_tsv(result: ResolutionResult, path) -> Path:
    path = Path(path)
    if result.streamed_counts is not None:
        if result.streamed_tsv is None or Path(result.streamed_tsv).resolve() != path.resolve():
            raise ValueError("a streamed (prodigal) result keeps no rows; pass out_tsv to resolve_gff")
        return path
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        fh.write("\t".join(G.GFF_RESOLUTION_COLUMNS) + "\n")
        for r in result.rows:
            fh.write(_resolution_line(r))
    return path


def write_cgc_gff(result: ResolutionResult, path, annotations: Optional[Dict[str, str]] = None) -> int:
    """Write cgc.gff (one 'gene' row per gene, contigs in file order, genes sorted by position)."""
    contig_order: Dict[str, int] = {}
    for r in result.rows:
        contig_order.setdefault(r.contig, len(contig_order))
    rows = sorted(result.rows, key=lambda r: (contig_order[r.contig], r.start, r.end))
    annotations = annotations or {}
    with Path(path).open("w", encoding="utf-8") as fh:
        for r in rows:
            if r.status == G.GFF_STATUS_NONCODING:
                ann = G.GFF_OTHER_PREFIX + r.feature_type
            else:
                ann = annotations.get(r.protein_id) or G.GFF_NULL_ANNOTATION
            strand = r.strand if r.strand in ("+", "-") else "."
            fh.write(
                f"{r.contig}\t.\t{G.GFF_GENE_FEATURE}\t{r.start}\t{r.end}\t.\t{strand}\t.\t"
                f"{G.GFF_PROTEIN_ID_ATTR_NCBI}={r.protein_id};{G.GFF_CGC_ANNOTATION_COL}={ann}\n"
            )
    return len(rows)


def rank_gff_types(gff_path, protein_ids: Sequence[str]) -> List[Tuple[str, float, int]]:
    """Match rate of every --gff_type on this GFF + FASTA: [(type, rate, matched genes)], best first."""
    ids = list(protein_ids)
    ranked = []
    for t in G.GFF_FORMAT_CHOICES:
        try:
            res = resolve_gff(gff_path, t, ids)
        except Exception as e:  # pragma: no cover - only used for the error message
            logger.debug(f"rank_gff_types: {t} failed: {e}")
            continue
        ranked.append((t, res.match_rate, res.n_matched))
    # among equal match rates prefer the original readers, then the specific resolver types
    pref = {t: i for i, t in enumerate(G.GFF_FORMAT_CHOICES)}
    ranked.sort(key=lambda x: (-x[2], -x[1], pref[x[0]]))
    return ranked


def format_resolution_failure(result: ResolutionResult, min_match: float, gff_path,
                              protein_ids: Sequence[str],
                              ranked: Optional[List[Tuple[str, float, int]]] = None) -> str:
    lines = [
        f"Only {result.n_matched}/{result.n_coding} protein-coding genes of {gff_path} "
        f"({100 * result.match_rate:.1f}%) were matched to a protein of the FASTA with "
        f"--gff_type {result.gff_type} (required: --gff_min_match {min_match:g}).",
    ]
    if result.n_coding == 0:
        lines.append("No protein-coding gene (gene/mRNA/transcript with CDS rows) was found in the GFF.")
    ex = result.unmatched_examples(5)
    if ex:
        lines.append("Unmatched genes (gene: identifiers looked up in the FASTA):")
        for r in ex:
            lines.append(f"  {r.gene_id} ({r.contig}:{r.start}-{r.end}): {', '.join(r.candidate_keys) or '-'}")
    if protein_ids:
        lines.append("First FASTA IDs: " + ", ".join(list(protein_ids)[:5]))
    if ranked:
        best = [(t, rate) for t, rate, n in ranked if n > 0 and t != result.gff_type]
        if best and best[0][1] > result.match_rate:
            t, rate = best[0]
            lines.append(f"Suggested: --gff_type {t} matches {100 * rate:.1f}% of the genes of this file.")
        else:
            lines.append(
                "No --gff_type matches this GFF and FASTA. The protein FASTA must be the one written by the "
                "same annotation run (FASTA IDs = transcript, CDS, protein or gene IDs of the GFF)."
            )
    lines.append(f"Supported --gff_type values: {', '.join(G.GFF_FORMAT_CHOICES)}. "
                 f"See {G.GFF_RESOLUTION_FILE} for every gene; use --gff_min_match 0 to continue anyway.")
    return "\n".join(lines)


def check_gff_resolution(gff_path, gff_type: str, protein_ids: Sequence[str],
                         min_match: float = G.GFF_MIN_MATCH_DEFAULT,
                         out_tsv=None) -> ResolutionResult:
    """Resolve, write ``out_tsv`` (optional), log a summary, raise GffResolutionError below min_match."""
    ids = list(protein_ids)
    result = resolve_gff(gff_path, gff_type, ids, out_tsv=out_tsv)
    logger.info(result.summary())
    if result.cross_gene_collisions:
        logger.warning(f"GFF resolution: {result.cross_gene_collisions} identifiers pointed to a protein "
                       "already used by another gene; the first gene keeps it.")
    min_match = float(min_match or 0.0)
    if min_match > 0 and (result.n_coding == 0 or result.match_rate < min_match):
        ranked = rank_gff_types(gff_path, ids)
        raise GffResolutionError(format_resolution_failure(result, min_match, gff_path, ids, ranked))
    if result.n_coding and result.match_rate < 1.0:
        logger.warning(
            f"GFF resolution: {result.n_coding - result.n_matched} protein-coding genes have no protein in the "
            f"FASTA (protein_id=unknown in cgc.gff); see {G.GFF_RESOLUTION_FILE}."
        )
    return result
