"""Shared building blocks of the dbCAN CGC figures.

Every CGC figure (circle / overview, linear CGC map, coverage, synteny) reads the run
outputs and colours genes through this module, so that

* there is ONE palette (``base_constants.CGC_FEATURE_COLORS``, the original dbCAN colours);
* a gene with several signature roles (``CAZyme|GH13_31+TC|8.A.9.1.1``) is filled with ONE
  colour, chosen by the CGCFinder Gene Type priority CAZyme > TC > TF > STP > Sulfatase >
  Peptidase; its other roles appear only in its text label (``GH13_31 + TC 8.A.9.1.1``);
* IDs that contain ``|`` are parsed with the CGC id pattern, never by split position.
"""
from __future__ import annotations

import csv
import logging
import math
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Patch, Polygon  # noqa: E402

import dbcan.constants.plots_constants as P  # noqa: E402

logger = logging.getLogger(__name__)

PALETTE = P.CGC_FEATURE_COLORS
OTHER = P.OTHER_GENE_TYPE

# 'contig|CGC12' or 'contig|CGC12|protein|type|family'; the contig may itself contain '|'
_CGC_ID_RE = re.compile(r"^(?P<contig>.+?)\|(?P<cgc>CGC\d+)(?:\|(?P<rest>.*))?$")
_ECAMI_SUFFIX = re.compile(r"_e\d+$")


class PlotInputError(RuntimeError):
    """The inputs needed for a requested figure are missing or inconsistent."""


def setup_matplotlib() -> None:
    plt.rcParams.update({
        "font.family": P.PLOT_FONT_FAMILY,
        "font.size": P.PLOT_BASE_FONT_SIZE,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })


# ---------------------------------------------------------------------------
# ids and gene types
# ---------------------------------------------------------------------------
def split_cgc_id(text: str) -> Optional[Tuple[str, str, str]]:
    """'NC_1|CGC3|NP_1.1|TC|2.A' -> ('NC_1', 'CGC3', 'NP_1.1|TC|2.A'); None if no CGC token."""
    m = _CGC_ID_RE.match(text or "")
    if not m:
        return None
    return m.group("contig"), m.group("cgc"), m.group("rest") or ""


def normalize_gene_type(name: str) -> str:
    """Map a signature type spelling (SULFATLAS, PEPTIDASE, cazyme ...) to a palette key."""
    return P.CGC_GENE_TYPE_ALIASES.get((name or "").strip().lower(), OTHER)


def parse_roles(annotation: str) -> List[Tuple[str, str]]:
    """'CAZyme|GH13_31+TC|8.A.9.1.1' -> [('CAZyme', 'GH13_31'), ('TC', '8.A.9.1.1')].

    Roles of the same type are merged (details joined by ';'); roles are ordered by the
    Gene Type priority; 'null' and unknown types are dropped.
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
    order = P.CGC_GENE_TYPE_PRIORITY
    return [(t, ";".join(d)) for t, d in sorted(merged.items(), key=lambda kv: order.get(kv[0], 99))]


def primary_gene_type(annotation: str, fallback: str = "") -> str:
    """The single colour type of a gene: highest-priority role, else the fallback, else Other."""
    roles = parse_roles(annotation)
    if roles:
        return roles[0][0]
    return normalize_gene_type(fallback) if fallback else OTHER


def gene_color(gtype: str) -> str:
    return PALETTE.get(gtype, PALETTE[OTHER])


def short_label(gtype: str, detail: str) -> str:
    """Compact family label of one role."""
    if gtype == "CAZyme":
        fams: List[str] = []
        for d in re.split(r"[|;]", detail):
            d = _ECAMI_SUFFIX.sub("", d.strip())
            if d and d not in fams and not any(f.startswith(d + "_") for f in fams):
                fams.append(d)
        return "+".join(fams) if fams else "CAZyme"
    first = detail.split(";")[0].strip()
    if gtype == "TC":
        return f"TC {first}" if first else "TC"
    if gtype == "Peptidase":
        return first.split(".")[0] or "Peptidase"
    if gtype == "TF":
        if first.isdigit():
            return f"TF {first}"
        m = re.match(r"^([A-Z0-9]+)_[A-Z0-9]{3,5}$", first)   # UniProt mnemonic CYNR_ECOLI
        return (m.group(1) if m else first)[:16] or "TF"
    return first[:16] or gtype


# ---------------------------------------------------------------------------
# data model + readers
# ---------------------------------------------------------------------------
@dataclass
class Gene:
    pid: str
    contig: str
    start: int
    end: int
    strand: str
    annotation: str = "null"
    gtype: str = OTHER
    roles: List[Tuple[str, str]] = field(default_factory=list)

    @classmethod
    def build(cls, pid, contig, start, end, strand, annotation, type_hint=""):
        g = cls(str(pid), str(contig), int(start), int(end), strand or ".", annotation or "null")
        g.roles = parse_roles(g.annotation)
        g.gtype = primary_gene_type(g.annotation, type_hint)
        return g

    @property
    def mid(self) -> float:
        return (self.start + self.end) / 2

    @property
    def label(self) -> str:
        """Family label; extra roles only appear here (never as extra colours)."""
        parts: List[str] = []
        for t, d in self.roles:
            lab = short_label(t, d)
            if lab.lower() not in (p.lower() for p in parts):
                parts.append(lab)
        return " + ".join(parts)


def read_cgc_table(path: Path) -> Dict[str, List[Gene]]:
    """cgc_standard_out.tsv -> {'contig|CGCn': [Gene, ...]} sorted by start."""
    path = Path(path)
    if not path.exists():
        raise PlotInputError(f"CGC table not found: {path}")
    out: Dict[str, List[Gene]] = {}
    with path.open(newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            cgc = (row.get(P.CGC_ID_FIELD) or "").strip()
            if not cgc or cgc.startswith("#"):
                continue
            g = Gene.build(row["Protein ID"], row["Contig ID"], row["Gene Start"], row["Gene Stop"],
                           row["Gene Strand"], row.get("Gene Annotation", "null"), row.get("Gene Type", ""))
            out.setdefault(f"{g.contig}|{cgc}", []).append(g)
    for genes in out.values():
        genes.sort(key=lambda g: (g.start, g.end))
    return out


def _gff_attrs(text: str) -> Dict[str, str]:
    return dict(kv.split("=", 1) for kv in text.strip().split(";") if "=" in kv)


def read_gff_genes(path: Path, feature_types: Iterable[str] = ("gene", "CDS")) -> Dict[str, List[Gene]]:
    """cgc.gff (or a dbCAN-PUL cgc.gff) -> {contig: [Gene, ...]} sorted by start.

    When a locus has both a ``gene`` and a ``CDS`` row, only one is kept.
    """
    path = Path(path)
    if not path.exists():
        raise PlotInputError(f"GFF not found: {path}")
    wanted = set(feature_types)
    out: Dict[str, List[Gene]] = {}
    seen = set()
    with path.open() as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in wanted:
                continue
            key = (f[0], f[3], f[4], f[6])
            if key in seen:
                continue
            seen.add(key)
            attrs = _gff_attrs(f[8])
            g = Gene.build(attrs.get(P.GFF_PROTEIN_ID_COL, "NA"), f[0], f[3], f[4], f[6],
                           attrs.get(P.GFF_CGC_ANNOTATION_COL, "null"))
            out.setdefault(f[0], []).append(g)
    for genes in out.values():
        genes.sort(key=lambda g: (g.start, g.end))
    return out


def contig_sizes(path: Path) -> Dict[str, int]:
    """##sequence-region lengths when present, else the last feature end per contig."""
    sizes: Dict[str, int] = {}
    regions: Dict[str, int] = {}
    with Path(path).open() as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if line.startswith("##sequence-region"):
                parts = line.split()
                if len(parts) >= 4:
                    regions[parts[1]] = int(parts[3])
                continue
            if line.startswith("#") or not line.strip():
                continue
            f = line.split("\t")
            if len(f) >= 5 and f[4].strip().isdigit():
                sizes[f[0]] = max(sizes.get(f[0], 0), int(f[4]))
    sizes.update(regions)
    return sizes


@dataclass
class SubstrateCall:
    cgcid: str
    pul: str = ""
    pul_substrate: str = ""
    sub_substrate: str = ""

    def names(self) -> List[str]:
        """Substrate names from dbCAN-PUL then dbCAN-sub, de-duplicated, no scores."""
        out: List[str] = []
        for text in (self.pul_substrate, self.sub_substrate):
            for name in re.split(r"[;,]", text or ""):
                name = name.strip()
                if name and name != "-" and name not in out:
                    out.append(name)
        return out


def read_substrates(path: Path) -> Dict[str, SubstrateCall]:
    """substrate_prediction.tsv -> {cgcid: SubstrateCall}; {} when the file does not exist."""
    path = Path(path)
    out: Dict[str, SubstrateCall] = {}
    if not path.exists():
        return out
    with path.open() as fh:
        next(fh, None)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if not f or not f[0].strip() or f[0].startswith("#"):
                continue
            f += [""] * (7 - len(f))
            out[f[0].strip()] = SubstrateCall(f[0].strip(), f[1].strip(), f[2].strip(), f[5].strip())
    return out


def read_depth(path: Path, contig: str, lo: int, hi: int) -> Dict[int, float]:
    """samtools-depth style file (contig, position, depth) -> {position: depth}.

    Only rows of ``contig`` inside [lo, hi] are kept (a depth file that covers several
    contigs used to be merged by position alone).
    """
    out: Dict[int, float] = {}
    with Path(path).open() as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3 or f[0] != contig:
                continue
            try:
                pos = int(f[1])
                if lo <= pos <= hi:
                    out[pos] = float(f[2])
            except ValueError:
                continue
    return out


# ---------------------------------------------------------------------------
# drawing primitives
# ---------------------------------------------------------------------------
def arrow_vertices(x0: float, x1: float, strand: str, y: float, h: float, head: float):
    hl = min(head, 0.6 * (x1 - x0))
    b, t = y - h / 2, y + h / 2
    if strand == "+":
        return [(x0, b), (x1 - hl, b), (x1, y), (x1 - hl, t), (x0, t)]
    if strand == "-":
        return [(x1, b), (x0 + hl, b), (x0, y), (x0 + hl, t), (x1, t)]
    return [(x0, b), (x1, b), (x1, t), (x0, t)]


def draw_gene(ax, x0: float, x1: float, strand: str, gtype: str, y: float, h: float, head: float,
              zorder: float = 3):
    """Solid, single-colour gene arrow (thin white outline separates abutting genes)."""
    poly = Polygon(arrow_vertices(x0, x1, strand, y, h, head), closed=True,
                   facecolor=gene_color(gtype), edgecolor=P.PLOT_GENE_GAP_COLOR, lw=0.6, zorder=zorder)
    ax.add_patch(poly)
    return poly


def text_widths(fig_w: float, ax_frac: float, xlo: float, xhi: float, texts: Sequence[str],
                fontsize: float, **kw) -> List[float]:
    """Width of each text in data units of an x axis spanning [xlo, xhi] over ax_frac of fig_w."""
    fig = plt.figure(figsize=(fig_w, 2))
    ax = fig.add_axes([0, 0, ax_frac, 1])
    ax.set_xlim(xlo, xhi)
    renderer = fig.canvas.get_renderer()
    inv = ax.transData.inverted()
    out = []
    for t in texts:
        art = ax.text(0, 0, t or " ", fontsize=fontsize, **kw)
        bb = art.get_window_extent(renderer=renderer)
        (a, _), (b, _) = inv.transform([(bb.x0, bb.y0), (bb.x1, bb.y0)])
        out.append(abs(b - a))
    plt.close(fig)
    return out


def assign_tiers(spans: Sequence[Tuple[float, float]], pad: float) -> List[int]:
    """Greedy interval packing: the lowest tier where the span does not overlap."""
    order = sorted(range(len(spans)), key=lambda i: spans[i][0])
    ends: List[float] = []
    out = [0] * len(spans)
    for i in order:
        x0, x1 = spans[i]
        for t, e in enumerate(ends):
            if x0 > e + pad:
                ends[t] = x1
                out[i] = t
                break
        else:
            ends.append(x1)
            out[i] = len(ends) - 1
    return out


@dataclass
class LabelLayout:
    labels: List[str]
    anchors: List[float]      # x of the gene the label belongs to
    xs: List[float]           # x of the label centre (clamped inside the axis)
    tiers: List[int]
    ntiers: int


def layout_labels(labels: Sequence[str], anchors: Sequence[float], xlo: float, xhi: float,
                  fig_w: float, ax_frac: float, fontsize: float) -> LabelLayout:
    widths = text_widths(fig_w, ax_frac, xlo, xhi, labels, fontsize)
    span = xhi - xlo
    xs = [min(max(a, xlo + w / 2), xhi - w / 2) for a, w in zip(anchors, widths)]
    keep = [i for i, lab in enumerate(labels) if lab]
    tiers = [0] * len(labels)
    packed = assign_tiers([(xs[i] - widths[i] / 2, xs[i] + widths[i] / 2) for i in keep], pad=0.012 * span)
    for i, t in zip(keep, packed):
        tiers[i] = t
    return LabelLayout(list(labels), list(anchors), xs, tiers, (max(packed) + 1) if packed else 0)


def draw_labels(ax, lay: LabelLayout, y0: float, tier_h: float, direction: int, fontsize: float,
                gene_edge: float, span: float):
    """Draw stacked labels from y0 outward (direction +1 = up, -1 = down) with leader lines."""
    for lab, anchor, x, t in zip(lay.labels, lay.anchors, lay.xs, lay.tiers):
        if not lab:
            continue
        y = y0 + direction * t * tier_h
        ax.text(x, y, lab, ha="center", va="bottom" if direction > 0 else "top", fontsize=fontsize,
                color=P.PLOT_INK, zorder=6)
        if t > 0 or abs(x - anchor) > 0.01 * span:
            yl = y - direction * 0.06
            ax.plot([anchor, anchor, x], [gene_edge, yl, yl], color=P.PLOT_LEADER_COLOR, lw=0.6, zorder=2)


def nice_step(span: float, target: int = 6) -> float:
    raw = span / max(1, target)
    if raw <= 0:
        return 1.0
    mag = 10 ** math.floor(math.log10(raw))
    for m in (1, 2, 2.5, 5, 10):
        if raw <= m * mag:
            return m * mag
    return 10 * mag


def fmt_bp(v: float) -> str:
    if abs(v) >= 1e6:
        return f"{v / 1e6:.2f}".rstrip("0").rstrip(".") + " Mb"
    if abs(v) >= 1e3:
        return f"{v / 1e3:.1f}".rstrip("0").rstrip(".") + " kb"
    return f"{int(round(v))} bp"


def fmt_kb_tick(v: float, step: float) -> str:
    if v >= 10e6 and step >= 1e5:
        return f"{v / 1e6:,.{0 if step % 1e6 == 0 else 1}f} Mb"
    digits = 0 if step % 1000 == 0 else (1 if step % 100 == 0 else 2)
    return f"{v / 1e3:,.{digits}f} kb"


def draw_ruler(ax, lo: float, hi: float, y: float, fontsize: float = P.PLOT_SMALL_FONT_SIZE,
               tick_h: float = 0.15):
    """Genomic coordinate ruler (round kb ticks) between lo and hi."""
    step = nice_step(hi - lo, 6)
    ax.plot([lo, hi], [y, y], color="#6E6E6E", lw=0.8, solid_capstyle="butt")
    tick = math.ceil(lo / step) * step
    while tick <= hi:
        ax.plot([tick, tick], [y, y - tick_h], color="#6E6E6E", lw=0.8)
        ax.text(tick, y - tick_h - 0.08, fmt_kb_tick(tick, step), ha="center", va="top",
                fontsize=fontsize, color=P.PLOT_INK_MUTED)
        tick += step


def draw_scale_bar(ax, x_right: float, y: float, span: float, fontsize: float = P.PLOT_SMALL_FONT_SIZE):
    length = nice_step(span, 5)
    x0 = x_right - length
    ax.plot([x0, x_right], [y, y], color=P.PLOT_INK, lw=1.3, solid_capstyle="butt", clip_on=False)
    for xx in (x0, x_right):
        ax.plot([xx, xx], [y - 0.08, y + 0.08], color=P.PLOT_INK, lw=1.0, clip_on=False)
    ax.text((x0 + x_right) / 2, y - 0.14, fmt_bp(length), ha="center", va="top", fontsize=fontsize,
            color=P.PLOT_INK)
    return length


def types_present(genes: Iterable[Gene]) -> List[str]:
    found = {g.gtype for g in genes}
    return [t for t in PALETTE if t in found]


def type_legend_handles(types: Sequence[str]):
    handles = [Patch(facecolor=gene_color(t), edgecolor="none") for t in types]
    return handles, list(types)


def add_type_legend(fig, types: Sequence[str], y: float = 0.0, x: float = 0.5, loc: str = "lower center",
                    fontsize: float = P.PLOT_LABEL_FONT_SIZE):
    if not types:
        return None
    hs, ls = type_legend_handles(types)
    return fig.legend(hs, ls, loc=loc, bbox_to_anchor=(x, y), ncol=len(hs), frameon=False,
                      fontsize=fontsize, handlelength=1.4, columnspacing=1.4, handletextpad=0.5)


def safe_name(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", text).strip("_") or "plot"


def normalize_formats(formats: Optional[Iterable[str]]) -> List[str]:
    out = []
    for f in formats or P.PLOT_DEFAULT_FORMATS:
        for part in str(f).split(","):
            part = part.strip().lower().lstrip(".")
            if not part:
                continue
            if part not in P.PLOT_SUPPORTED_FORMATS:
                raise PlotInputError(f"Unsupported figure format '{part}' (use {', '.join(P.PLOT_SUPPORTED_FORMATS)})")
            if part not in out:
                out.append(part)
    return out or list(P.PLOT_DEFAULT_FORMATS)


def save_figure(fig, out_dir: Path, template: str, name_fields: Dict[str, str], formats, dpi: int) -> List[Path]:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    written = []
    for ext in normalize_formats(formats):
        path = out_dir / template.format(ext=ext, **{k: safe_name(v) for k, v in name_fields.items()})
        fig.savefig(path, dpi=dpi, facecolor="white")
        written.append(path)
        logger.info(f"Saved figure: {path}")
    plt.close(fig)
    return written


def text_width_data_in(text: str, fontsize: float, **kw) -> float:
    """Rendered width of ``text`` in inches."""
    fig = plt.figure(figsize=(4, 1))
    art = fig.text(0, 0, text or " ", fontsize=fontsize, **kw)
    bb = art.get_window_extent(renderer=fig.canvas.get_renderer())
    plt.close(fig)
    return bb.width / fig.dpi


def middle_truncate(text: str, n: int) -> str:
    return text if len(text) <= n else text[: n // 2 - 1] + "…" + text[-(n - n // 2):]


def fmt_axis_tick(v: float, step: float, total: float) -> str:
    """Tick label of a genome axis: Mb for sequences of 2 Mb or more, kb otherwise."""
    if total >= 2e6:
        return "0" if v == 0 else f"{v / 1e6:.{0 if step % 1e6 == 0 else 1}f} Mb"
    return "0" if v == 0 else fmt_kb_tick(v, step)


def plural(n: int, word: str) -> str:
    return f"{n} {word}{'' if n == 1 else 's'}"
