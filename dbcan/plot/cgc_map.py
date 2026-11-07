"""Linear CGC map (``dbcan_plot CGC_plot``) and its read-depth variant (``CGC_coverage_plot``).

Layout, top to bottom: title (``CGCn  contig:start–end``) and one substrate line (names
only), optional read-depth panels (one per sample, shared y axis), the gene track (arrows
to scale by strand, one colour per gene, stacked family labels), a kb ruler with a scale
bar, and the gene-type legend (only the types present). Flanking genes are not drawn.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

import numpy as np
import matplotlib.pyplot as plt

import dbcan.constants.plots_constants as P
from dbcan.plot import cgc_viz as V

logger = logging.getLogger(__name__)

U = 0.27            # inches per y data unit in the gene axis
GENE_H = 1.0
L_IN, R_IN = 0.45, 0.4
TITLE_IN = 0.62
LEGEND_IN = 0.42
PANEL_IN = 0.85     # height of one read-depth panel
PANEL_GAP_IN = 0.14


def substrate_text(names: Optional[Sequence[str]]) -> str:
    return "Substrate: " + (", ".join(names) if names else "not predicted")


def cgc_title(cgcid: str, genes: Sequence[V.Gene]) -> str:
    parts = V.split_cgc_id(cgcid)
    name = parts[1] if parts else cgcid
    lo, hi = min(g.start for g in genes), max(g.end for g in genes)
    return f"{name}  {genes[0].contig}:{lo:,}–{hi:,}"


def load_depth_tracks(depth_files: Sequence[Tuple[str, Path]], contig: str, lo: int, hi: int):
    """[(name, path)] -> [(name, xs, ys)] for the window [lo, hi] of ``contig``."""
    xs = np.arange(lo, hi + 1)
    tracks = []
    for name, path in depth_files:
        path = Path(path)
        if not path.exists():
            raise V.PlotInputError(f"Read-depth file not found: {path}")
        depth = V.read_depth(path, contig, lo, hi)
        if not depth:
            logger.warning(f"{path}: no depth rows for contig '{contig}' in {lo}-{hi}; drawing zero depth.")
        tracks.append((name, xs, np.array([depth.get(int(p), 0.0) for p in xs], dtype=float)))
    return tracks


def draw_depth_panels(fig, tracks, genes: Sequence[V.Gene], xlo: float, xhi: float, W: float, H: float,
                      bottom_in: float, left_in: float = L_IN, right_in: float = R_IN):
    """One panel per sample stacked upward from ``bottom_in``; first sample on top; shared y."""
    ymax = max([1.0] + [float(t[2].max()) for t in tracks if len(t[2])]) * 1.08
    ax_w = W - left_in - right_in
    axes = []
    for i, (name, xs, ys) in enumerate(tracks):
        j = len(tracks) - 1 - i
        b = bottom_in + j * (PANEL_IN + PANEL_GAP_IN)
        cax = fig.add_axes([left_in / W, b / H, ax_w / W, PANEL_IN / H])
        for g in genes:
            cax.axvspan(g.start, g.end, color=P.COVERAGE_GENE_SPAN_COLOR, zorder=0, lw=0)
        cax.fill_between(xs, ys, 0, step="mid", facecolor=P.COVERAGE_FILL_COLOR,
                         alpha=P.COVERAGE_FILL_ALPHA, lw=0, zorder=1)
        cax.plot(xs, ys, drawstyle="steps-mid", color=P.COVERAGE_LINE_COLOR,
                 alpha=P.COVERAGE_LINE_ALPHA, lw=0.7, zorder=2)
        cax.set_xlim(xlo, xhi)
        cax.set_ylim(0, ymax)
        cax.set_xticks([])
        for s in ("top", "right", "bottom"):
            cax.spines[s].set_visible(False)
        cax.spines["left"].set_color("#8A8A8A")
        cax.tick_params(axis="y", labelsize=P.PLOT_SMALL_FONT_SIZE, colors=P.PLOT_INK_MUTED, length=2)
        cax.set_ylabel(P.COVERAGE_YLABEL, fontsize=P.PLOT_SMALL_FONT_SIZE, color=P.PLOT_INK_MUTED)
        if len(tracks) > 1 or name:
            cax.text(0.005, 0.97, name, transform=cax.transAxes, ha="left", va="top",
                     fontsize=P.PLOT_LABEL_FONT_SIZE, color=P.PLOT_INK)
        axes.append(cax)
    return axes


def plot_cgc_map(genes: Sequence[V.Gene], cgcid: str, out_dir: Path, *, substrates: Optional[Sequence[str]] = None,
                 depth_files: Optional[Sequence[Tuple[str, Path]]] = None, show_protein_ids: bool = False,
                 formats=P.PLOT_DEFAULT_FORMATS, dpi: int = P.PLOT_DEFAULT_DPI,
                 template: Optional[str] = None) -> List[Path]:
    """Draw one CGC; with ``depth_files`` the read-depth panels are added above the genes."""
    if not genes:
        raise V.PlotInputError(f"CGC {cgcid} has no genes")
    V.setup_matplotlib()
    genes = sorted(genes, key=lambda g: g.start)
    contig = genes[0].contig
    cstart, cend = min(g.start for g in genes), max(g.end for g in genes)
    span = max(1, cend - cstart)
    xlo, xhi = cstart - 0.015 * span, cend + 0.015 * span
    head = 0.018 * span

    W = max(8.0, min(16.0, 3.0 + 0.6 * len(genes)))
    ax_w = W - L_IN - R_IN
    fs = P.PLOT_LABEL_FONT_SIZE
    lay = V.layout_labels([g.label for g in genes], [g.mid for g in genes], xlo, xhi, W, ax_w / W, fs)
    tier_h = (fs * 1.25 / 72) / U + 0.08
    lab_base = GENE_H / 2 + 0.3
    y_top = lab_base + max(1, lay.ntiers) * tier_h + 0.1

    y_cursor = -GENE_H / 2 - 0.2
    pid_h = 0.0
    if show_protein_ids:
        longest = max(len(g.pid) for g in genes)
        pid_h = (longest * P.PLOT_SMALL_FONT_SIZE * 0.55 / 72) * 0.5 / U + 0.15
        y_pid = y_cursor
        y_cursor -= pid_h
    y_ruler = y_cursor - 0.25
    y_scale = y_ruler - 1.05
    y_bottom = y_scale - 0.7

    ax_h = (y_top - y_bottom) * U
    tracks = []
    if depth_files:
        tracks = load_depth_tracks(depth_files, contig, int(np.floor(xlo)), int(np.ceil(xhi)))
    depth_in = len(tracks) * (PANEL_IN + PANEL_GAP_IN)
    H = TITLE_IN + depth_in + ax_h + LEGEND_IN
    fig = plt.figure(figsize=(W, H))
    ax = fig.add_axes([L_IN / W, LEGEND_IN / H, ax_w / W, ax_h / H])
    ax.set_xlim(xlo, xhi)
    ax.set_ylim(y_bottom, y_top)
    ax.axis("off")

    ax.plot([cstart, cend], [0, 0], color=P.PLOT_BACKBONE_COLOR, lw=1.0, zorder=1, solid_capstyle="butt")
    for g in genes:
        V.draw_gene(ax, g.start, g.end, g.strand, g.gtype, 0, GENE_H, head)
    V.draw_labels(ax, lay, lab_base, tier_h, +1, fs, GENE_H / 2 + 0.03, span)
    if show_protein_ids:
        for g in genes:
            ax.text(g.mid, y_pid, g.pid, rotation=30, rotation_mode="anchor", ha="right", va="top",
                    fontsize=P.PLOT_SMALL_FONT_SIZE, color=P.PLOT_INK_MUTED)
    V.draw_ruler(ax, cstart, cend, y_ruler)
    V.draw_scale_bar(ax, cend, y_scale, span)

    if tracks:
        draw_depth_panels(fig, tracks, genes, xlo, xhi, W, H, LEGEND_IN + ax_h + PANEL_GAP_IN)

    fig.text(L_IN / W, 1 - 0.18 / H, cgc_title(cgcid, genes), ha="left", va="top",
             fontsize=P.PLOT_TITLE_FONT_SIZE, fontweight="bold", color=P.PLOT_INK)
    fig.text(L_IN / W, 1 - 0.42 / H, substrate_text(substrates), ha="left", va="top",
             fontsize=P.PLOT_SUBTITLE_FONT_SIZE, color=P.PLOT_INK_MUTED)
    V.add_type_legend(fig, V.types_present(genes), y=0.02 / H)

    if template is None:
        template = P.CGC_COVERAGE_FILE_TEMPLATE if tracks else P.CGC_MAP_FILE_TEMPLATE
    return V.save_figure(fig, out_dir, template, {"cgc": cgcid}, formats, dpi)
