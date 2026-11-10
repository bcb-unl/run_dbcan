"""CGC model and fixed-width gene-track geometry used by the expression figures.

``dbcan_plot`` CGC figures were redesigned on :mod:`dbcan.plot.cgc_viz` /
:mod:`dbcan.plot.cgc_map` (one palette, one renderer per plot type), which removed the
old helpers from :mod:`dbcan.utils.plots`. The SeMa-Trap style expression panels in
:mod:`dbcan.plot.expression_plots` are drawn on the old 1000-unit pixel geometry
(expression trajectories are aligned to it), so the helpers they need live here.

TODO: move the expression gene panel onto :class:`dbcan.plot.cgc_map.GeneTrack` and drop
this module.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Sequence

from dbcan.utils.utils import cgc_standard_line

# Geometry of the gene track (data units; the x axis is GENOME_PLOT_WIDTH wide).
TRACK_WIDTH = 1000
TRACK_HEIGHT = 160
GENE_HALF_HEIGHT = 10
ARROW_HEAD = 4
SCALE_TICKS = 10


class CGC_Standard_Out:
    """Rows of cgc_standard_out.tsv as :class:`cgc_standard_line` records (file order)."""

    def __init__(self, filename: str):
        self.genes: List[cgc_standard_line] = []
        with open(filename) as fh:
            next(fh, None)
            for line in fh:
                if not line.strip() or line.startswith("CGC#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 8:
                    continue
                self.genes.append(cgc_standard_line(parts))

    def __iter__(self):
        return iter(self.genes)

    def CGCID2genes(self) -> Dict[str, List[cgc_standard_line]]:
        cgcdict: Dict[str, List[cgc_standard_line]] = {}
        for gene in self:
            cgcdict.setdefault(gene.cgcid, []).append(gene)
        return cgcdict


class CGC:
    """The genes of one CGC (cgcid 'contig|CGCn')."""

    def __init__(self, genes: Sequence[cgc_standard_line]):
        self.genes = list(genes)
        self.ID = self.genes[0].cgcid
        self.start = min(g.gene_start for g in self.genes)
        self.end = max(g.gene_end for g in self.genes)
        self.gene_num = len(self.genes)

    def __iter__(self):
        return iter(self.genes)

    def __repr__(self):
        return "\t".join([self.ID, str(self.start), str(self.end), str(self.gene_num)])

    def __len__(self):
        return len(self.genes)

    def get_positions(self):
        return ([g.gene_start for g in self], [g.gene_end for g in self], [g.strand for g in self])

    def get_proteinID(self) -> List[str]:
        return [g.seqid for g in self]

    def get_cgc_CAZyme(self) -> List[str]:
        return [g.gene_type for g in self]

    def get_annotations(self) -> List[str]:
        return [getattr(g, "annotation", "") or g.gene_type for g in self]


class CGC_standard_out_2CGC:
    """cgc_standard_out.tsv grouped into :class:`CGC` objects."""

    def __init__(self, dbcan: CGC_Standard_Out):
        self.CGCs = [CGC(genes) for genes in dbcan.CGCID2genes().values()]

    def __iter__(self):
        return iter(self.CGCs)

    def cgcid2CGC(self) -> Dict[str, CGC]:
        return {cgc.ID: cgc for cgc in self}


def Get_Position(starts, ends, strands, labels, yshift: float = 0):
    """Gene polygons, connector lines, labels and scale of one CGC on the fixed-width track.

    Returns ``(polygons, lines, texts, scale_positions, scale_text)``; polygons and lines are
    space-separated coordinate strings ("x1 y1 x2 y2 ...").
    """
    plot_start_y = TRACK_HEIGHT / 2 - GENE_HALF_HEIGHT - yshift
    s, e, st = list(starts), list(ends), list(strands)
    shift_pos = min(s) if s else 0
    maxbp = (max(e) - min(s)) if s and e else 1
    pix = TRACK_WIDTH / maxbp if maxbp > 0 else 1.0
    s = [v - shift_pos for v in s]
    e = [v - shift_pos for v in e]
    h, tri = GENE_HALF_HEIGHT, ARROW_HEAD

    polygons, lines, texts = [], [], []
    for i in range(len(s)):
        x0, x1 = s[i] * pix, e[i] * pix
        y0, ym, y1 = plot_start_y, plot_start_y + h, plot_start_y + 2 * h
        if st[i] == "+":
            pts = [(x0, y0), (x1 - tri, y0), (x1, ym), (x1 - tri, y1), (x0, y1)]
        else:
            pts = [(x0, ym), (x0 + tri, y0), (x1, y0), (x1, y1), (x0 + tri, y1)]
        polygons.append(" ".join(f"{x} {y}" for x, y in pts))
        if i < len(s) - 1:
            lines.append(f"{x1} {ym} {s[i + 1] * pix} {ym}")
        texts.append(labels[i].split(".")[0] if i < len(labels) else str(i + 1))

    each_bp = maxbp / SCALE_TICKS
    each_pix = each_bp * pix
    y_scale = plot_start_y - 50 + 3 * h
    scale_positions = [f"0 {y_scale} {SCALE_TICKS * each_pix} {y_scale}"]
    scale_text = []
    y_tick = y_scale - 1
    for i in range(SCALE_TICKS + 1):
        scale_positions.append(f"{i * each_pix} {y_tick} {i * each_pix} {y_tick + 0.6 * h}")
        scale_text.append(str(int(each_bp * i) + shift_pos))
    return polygons, lines, texts, scale_positions, scale_text


def _points2(coord: str):
    x1, y1, x2, y2 = coord.split()
    return float(x1), float(y1), float(x2), float(y2)


def plot_genome_line(lines, ax) -> None:
    """Connectors between adjacent genes."""
    for line in lines:
        x1, y1, x2, y2 = _points2(line)
        ax.plot([x1, x2], [y1, y2], color="gray", lw=2, solid_capstyle="round", zorder=5)


def plot_scale_line(lines, label, ax) -> None:
    """Horizontal scale with tick marks and bp labels."""
    for i, line in enumerate(lines):
        x1, y1, x2, y2 = _points2(line)
        ax.plot([x1, x2], [y1, y2], color="gray", lw=2, zorder=5)
        if i >= 1:
            ax.text(x1, y1 - 6, label[i - 1], va="top", ha="center", fontsize=8)


def add_gene_annotations(ax, starts, ends, annotations, font_size: float = 10, rotation: float = 60,
                         y_offset: float = 25) -> Optional[float]:
    """Annotation text above each gene arrow; returns the highest y used (for the y limit)."""
    if not starts or not ends:
        return None
    shift_pos = min(starts)
    maxbp = max(ends) - shift_pos
    if maxbp <= 0:
        return None
    pix = TRACK_WIDTH / maxbp
    label_base_y = TRACK_HEIGHT / 2 - GENE_HALF_HEIGHT + 2 * GENE_HALF_HEIGHT + y_offset
    for s, e, ann in zip(starts, ends, annotations):
        ax.text(((s - shift_pos) + (e - shift_pos)) / 2 * pix, label_base_y, ann,
                ha="center", va="bottom", fontsize=font_size, rotation=rotation)
    return label_base_y + font_size * 0.62
