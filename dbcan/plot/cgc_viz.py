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

import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Polygon

import dbcan.constants.plots_constants as P
from dbcan.IO.cgc_ids import normalize_gene_type, parse_roles, split_cgc_id, split_cgc_key  # noqa: F401 (re-exported)

logger = logging.getLogger(__name__)

PALETTE = P.CGC_FEATURE_COLORS
OTHER = P.OTHER_GENE_TYPE

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
# split_cgc_id, normalize_gene_type and parse_roles live in dbcan.IO.cgc_ids (shared with
# substrate_prediction and dbcan_utils).


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
        # dbCAN-PUL gffs give 'TF|0036286+TF|HTH_AraC': prefer the named family over a bare number
        named = [d.strip() for d in detail.split(";") if d.strip() and not d.strip().isdigit()]
        if not named:
            return f"TF {first}" if first else "TF"
        first = named[0]
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
    gene_type: str = ""     # 'Gene Type' column of cgc_standard_out.tsv as written (null, CAZyme, ...)

    @classmethod
    def build(cls, pid, contig, start, end, strand, annotation, type_hint=""):
        g = cls(str(pid), str(contig), int(start), int(end), strand or ".", annotation or "null")
        g.gene_type = (type_hint or "").strip()
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
    """samtools-depth style file (contig, position, depth; any whitespace) -> {position: depth}.

    Only rows of ``contig`` inside [lo, hi] are kept (a depth file that covers several
    contigs used to be merged by position alone). Raises :class:`PlotInputError` when the
    file has no row at all for ``contig`` (wrong file, or contig names that do not match the
    cgc.gff), since the figure would otherwise show a flat zero coverage.
    """
    out: Dict[int, float] = {}
    seen_contig = False
    other: List[str] = []
    with Path(path).open() as fh:
        for line in fh:
            f = line.split()
            if len(f) < 3 or f[0].startswith("#"):
                continue
            if f[0] != contig:
                if len(other) < 3 and f[0] not in other:
                    other.append(f[0])
                continue
            seen_contig = True
            try:
                pos = int(f[1])
                if lo <= pos <= hi:
                    out[pos] = float(f[2])
            except ValueError:
                continue
    if not seen_contig:
        found = ", ".join(repr(c) for c in other) or "none"
        raise PlotInputError(
            f"{path}: no read-depth rows for contig '{contig}' (contigs found: {found}). "
            "The contig names of the depth file must match the cgc.gff.")
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
    anchors: List[float]      # x where the leader leaves the gene (inside the gene's extent)
    xs: List[float]           # x of the label centre (clamped inside the axis)
    tiers: List[int]
    ntiers: int
    widths: List[float] = field(default_factory=list)

    def box(self, i: int) -> Tuple[float, float]:
        w = self.widths[i] if self.widths else 0.0
        return self.xs[i] - w / 2, self.xs[i] + w / 2


def _free_points(lo: float, hi: float, blocked: Sequence[Tuple[float, float]]) -> List[Tuple[float, float]]:
    """Sub-intervals of [lo, hi] not covered by any of the ``blocked`` intervals."""
    free = [(lo, hi)]
    for b0, b1 in blocked:
        nxt = []
        for f0, f1 in free:
            if b1 <= f0 or b0 >= f1:
                nxt.append((f0, f1))
                continue
            if b0 > f0:
                nxt.append((f0, b0))
            if b1 < f1:
                nxt.append((b1, f1))
        free = nxt
    return [(a, b) for a, b in free if b >= a]


def layout_labels(labels: Sequence[str], anchors: Sequence[float], xlo: float, xhi: float,
                  fig_w: float, ax_frac: float, fontsize: float,
                  extents: Optional[Sequence[Tuple[float, float]]] = None) -> LabelLayout:
    """Stack gene labels in tiers so that neither labels nor leader lines cross any text.

    A label in tier ``t > 0`` is joined to its gene by a vertical leader through tiers
    ``0 .. t-1``. The leader may leave the gene anywhere inside its extent (``extents``;
    default the anchor point); the label is centred over the leader, or starts / ends at
    it when centring does not fit. Tiers are filled bottom-up, left to right. A label is
    put in the current tier only when (1) it does not overlap a label of that tier, (2) its
    leader passes no label of a lower tier, and (3) every label still waiting keeps some
    point of its gene that no label covers, so its leader can always rise to a later tier.
    """
    widths = text_widths(fig_w, ax_frac, xlo, xhi, labels, fontsize)
    span = xhi - xlo
    if extents is None:
        extents = [(a, a) for a in anchors]
    tiers, leaders, xs, forced = _pack_labels(labels, widths, anchors, extents, xlo, xhi, span)
    if forced:
        logger.debug(f"{forced} gene label(s) could not avoid every leader line")
    used = sorted({tiers[i] for i, lab in enumerate(labels) if lab})
    remap = {t: k for k, t in enumerate(used)}      # drop empty tiers
    tiers = [remap.get(t, 0) for t in tiers]
    return LabelLayout(list(labels), leaders, xs, tiers, len(used), widths)


def _pack_labels(labels, widths, anchors, extents, xlo, xhi, span):
    """Tier packing of :func:`layout_labels`; returns (tiers, leader xs, centres, n forced)."""
    n = len(widths)
    pad = 0.012 * span            # gap between two labels of one tier
    lpad = 0.004 * span           # clearance between a leader and a label box
    xs = [min(max(a, xlo + w / 2), xhi - w / 2) for a, w in zip(anchors, widths)]
    leaders = list(anchors)
    tiers = [0] * n
    blocked: List[Tuple[float, float]] = []       # label boxes of the finished tiers (+ clearance)

    def leader_room(i: int, boxes: Sequence[Tuple[float, float]]) -> List[Tuple[float, float]]:
        a = anchors[i]
        e0, e1 = extents[i]
        inset = 0.15 * (e1 - e0)
        return _free_points(min(e0 + inset, a), max(e1 - inset, a), boxes)

    def leader_points(i: int) -> List[float]:
        a = anchors[i]
        pts = []
        for f0, f1 in leader_room(i, blocked):
            pts += [min(max(a, f0), f1), (f0 + f1) / 2, f0, f1]
        out = []
        for p in sorted(pts, key=lambda p: abs(p - a)):
            if all(not (b0 < p < b1) for b0, b1 in blocked) and p not in out:
                out.append(p)
        return out

    def centre(i: int, lx: float, shift: int) -> float:
        # shift 0: centred over the leader; +1 / -1: the label starts / ends just at the leader
        cx = lx + shift * 0.42 * widths[i]
        return min(max(cx, xlo + widths[i] / 2), xhi - widths[i] / 2)

    waiting = sorted((i for i, lab in enumerate(labels) if lab), key=lambda i: anchors[i])
    forced = 0
    t = 0
    while waiting:
        row: List[Tuple[float, float]] = []                 # boxes of tier t
        placed_any = False
        for i in list(waiting):
            pts = [anchors[i]] if t == 0 else leader_points(i)
            choice = None
            for lx in pts:
                for shift in (0, 1, -1):
                    cx = centre(i, lx, shift)
                    b0, b1 = cx - widths[i] / 2, cx + widths[i] / 2
                    if any(not (b1 + pad <= r0 or b0 >= r1 + pad) for r0, r1 in row):
                        continue                              # (1)
                    newbox = (b0 - lpad, b1 + lpad)
                    if any(j != i and not leader_room(j, blocked + row + [newbox]) for j in waiting):
                        continue                              # (3)
                    choice = (lx, cx, b0, b1)
                    break
                if choice:
                    break
            if choice is None:
                continue
            lx, cx, b0, b1 = choice
            tiers[i], leaders[i], xs[i] = t, lx, cx
            row.append((b0, b1))
            waiting.remove(i)
            placed_any = True
        if not placed_any:              # cannot happen with a free leader point; never loop forever
            i = waiting.pop(0)
            pts = leader_points(i) or [anchors[i]]
            tiers[i], leaders[i], xs[i] = t, pts[0], centre(i, pts[0], 0)
            row.append((xs[i] - widths[i] / 2, xs[i] + widths[i] / 2))
            forced += 1
        blocked += [(r0 - lpad, r1 + lpad) for r0, r1 in row]
        t += 1
    return tiers, leaders, xs, forced


def draw_labels(ax, lay: LabelLayout, y0: float, tier_h: float, direction: int, fontsize: float,
                gene_edge: float, span: float):
    """Draw stacked labels from y0 outward (direction +1 = up, -1 = down) with leader lines."""
    for i, (lab, anchor, x, t) in enumerate(zip(lay.labels, lay.anchors, lay.xs, lay.tiers)):
        if not lab:
            continue
        y = y0 + direction * t * tier_h
        ax.text(x, y, lab, ha="center", va="bottom" if direction > 0 else "top", fontsize=fontsize,
                color=P.PLOT_INK, zorder=6)
        if t > 0 or abs(x - anchor) > 0.01 * span:
            yl = y - direction * 0.06
            b0, b1 = lay.box(i)
            inset = 0.02 * (b1 - b0)
            # vertical leader up to just under the label; a horizontal part only when the
            # leader ends beside the label (a label pushed in from the axis edge)
            x_end = min(max(anchor, b0 + inset), b1 - inset)
            ax.plot([anchor, anchor, x_end], [gene_edge, yl, yl], color=P.PLOT_LEADER_COLOR, lw=0.6, zorder=2)


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
