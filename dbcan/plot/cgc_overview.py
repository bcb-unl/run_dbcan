"""Assembly overview of CGCs (``run_dbcan cgc_circle_plot`` for multi-contig inputs).

A tidy ideogram: one row per contig that carries at least one CGC (all rows share one bp
scale), the contig name on the left and its length at the end of the bar, CGC blocks on
the bar (CAZyme colour when a substrate was predicted, grey otherwise) and labels only for
the CGCs with a substrate prediction (``CGC6 alpha-glucan``). When no substrate prediction
exists the CGC ids are shown where they fit.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List, Sequence

import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle

import dbcan.constants.plot_cgc_circle_constants as C
import dbcan.constants.plots_constants as P
from dbcan.plot import cgc_viz as V

logger = logging.getLogger(__name__)

ROW_IN = 0.46          # height of one contig row
BAR_IN = 0.09          # contig bar thickness
CGC_IN = 0.2           # CGC block height
TITLE_IN = 0.62
FOOT_IN = 0.5
FIG_W = 11.0
R_IN = 0.75            # room for the length label after the longest bar


def _cgc_span(genes: Sequence[V.Gene]):
    return min(g.start for g in genes), max(g.end for g in genes)


def plot_overview(cgcs: Dict[str, List[V.Gene]], sizes: Dict[str, int], calls: Dict[str, V.SubstrateCall],
                  out_dir: Path, *, formats=("png", "pdf"), dpi: int = P.PLOT_DEFAULT_DPI,
                  max_contigs: int = C.OVERVIEW_MAX_CONTIGS) -> List[Path]:
    V.setup_matplotlib()
    by_contig: Dict[str, List[str]] = {}
    for cid, genes in cgcs.items():
        by_contig.setdefault(genes[0].contig, []).append(cid)
    for c, ids in by_contig.items():
        sizes.setdefault(c, max(g.end for cid in ids for g in cgcs[cid]))
    contigs = sorted(by_contig, key=lambda c: (-sizes[c], c))
    hidden = contigs[max_contigs:]
    contigs = contigs[:max_contigs]
    if hidden:
        logger.warning(f"Overview shows the {max_contigs} longest of {len(by_contig)} contigs with CGCs.")
    n_pred = sum(1 for cid in cgcs if calls.get(cid) and calls[cid].names())
    label_all = not any(c.names() for c in calls.values())

    fs = P.PLOT_LABEL_FONT_SIZE
    name_w = max([V.text_width_data_in(V.middle_truncate(c, 28), fs) for c in contigs] + [0.5])
    L_IN = 0.35 + name_w + 0.15
    ax_w = FIG_W - L_IN - R_IN
    wmax = max([sizes[c] for c in contigs] + [1])
    xlo, xhi = -0.005 * wmax, 1.005 * wmax

    # labels: substrate-predicted CGCs first (up to 2 tiers), then ids where they fit (1 tier)
    rows = []
    for c in contigs:
        items = []
        for cid in by_contig[c]:
            a, b = _cgc_span(cgcs[cid])
            names = calls[cid].names() if cid in calls else []
            name = V.split_cgc_id(cid)[1] if V.split_cgc_id(cid) else cid
            items.append(dict(cid=cid, a=a, b=b, pred=bool(names),
                              text=f"{name} {', '.join(names)}" if names else name))
        rows.append(items)
    all_items = [it for items in rows for it in items]
    widths = V.text_widths(FIG_W, ax_w / FIG_W, xlo, xhi, [it["text"] for it in all_items], fs - 0.5)
    for it, w in zip(all_items, widths):
        it["w"] = w
        it["x"] = min(max((it["a"] + it["b"]) / 2, xlo + w / 2), xhi - w / 2)
    pad = 0.01 * wmax
    row_tiers = []
    for items in rows:
        for it in sorted(items, key=lambda i: (not i["pred"], i["x"])):
            it["tier"] = None
            if not (it["pred"] or label_all):
                continue
            l, r = it["x"] - it["w"] / 2, it["x"] + it["w"] / 2
            limit = 2 if it["pred"] else 1
            for t in range(limit):
                placed_spans = [o for o in items if o.get("tier") == t and o is not it]
                if all(r + pad < o["x"] - o["w"] / 2 or l - pad > o["x"] + o["w"] / 2 for o in placed_spans):
                    it["tier"] = t
                    break
        row_tiers.append(max([it["tier"] for it in items if it["tier"] is not None] + [0]) + 1)

    tier_in = fs * 1.3 / 72
    heights = [ROW_IN + (t - 1) * tier_in for t in row_tiers]
    body_in = sum(heights)
    H = TITLE_IN + body_in + FOOT_IN
    fig = plt.figure(figsize=(FIG_W, H))
    ax = fig.add_axes([L_IN / FIG_W, FOOT_IN / H, ax_w / FIG_W, body_in / H])
    ax.set_xlim(xlo, xhi)
    ax.set_ylim(-body_in, 0)       # y in inches, downwards
    ax.axis("off")

    y_top = 0.0
    for c, items, h, nt in zip(contigs, rows, heights, row_tiers):
        yb = y_top - h + 0.12                       # bar centre near the bottom of the row
        L = sizes[c]
        ax.add_patch(Rectangle((0, yb - BAR_IN / 2), L, BAR_IN, facecolor=C.OVERVIEW_CONTIG_COLOR,
                               edgecolor="none", zorder=1))
        fig.text((L_IN - 0.15) / FIG_W, (FOOT_IN + body_in + yb) / H, V.middle_truncate(c, 28), ha="right",
                 va="center", fontsize=fs, color=P.PLOT_INK)
        ax.text(L + 0.006 * wmax, yb, V.fmt_bp(L), ha="left", va="center", fontsize=fs - 0.8,
                color=P.PLOT_INK_MUTED)
        for it in items:
            wbox = max(it["b"] - it["a"], 0.002 * wmax)
            xb = (it["a"] + it["b"]) / 2 - wbox / 2
            color = V.gene_color("CAZyme") if it["pred"] else C.OVERVIEW_CGC_COLOR
            ax.add_patch(Rectangle((xb, yb - CGC_IN / 2), wbox, CGC_IN, facecolor=color, edgecolor="none",
                                   zorder=3 if it["pred"] else 2))
            if it["tier"] is None:
                continue
            ly = yb + CGC_IN / 2 + 0.04 + it["tier"] * tier_in
            ax.text(it["x"], ly, it["text"], ha="center", va="bottom", fontsize=fs - 0.5,
                    color=P.PLOT_INK if it["pred"] else P.PLOT_INK_MUTED, zorder=5)
            mid = (it["a"] + it["b"]) / 2
            if it["tier"] > 0 or abs(it["x"] - mid) > 0.01 * wmax:
                ax.plot([mid, mid, it["x"]], [yb + CGC_IN / 2, ly - 0.02, ly - 0.02],
                        color=P.PLOT_LEADER_COLOR, lw=0.5, zorder=2)
        y_top -= h

    n_contigs = len(by_contig)
    fig.text(0.35 / FIG_W, 1 - 0.18 / H, "CGC overview", ha="left", va="top", fontsize=P.PLOT_TITLE_FONT_SIZE,
             fontweight="bold", color=P.PLOT_INK)
    sub = (f"{len(cgcs)} CGC{'s' if len(cgcs) != 1 else ''} on {n_contigs} contig{'s' if n_contigs != 1 else ''}"
           f" · {n_pred} with substrate prediction")
    if hidden:
        sub += f" · {len(hidden)} shorter contig{'s' if len(hidden) != 1 else ''} with CGCs not shown"
    fig.text(0.35 / FIG_W, 1 - 0.42 / H, sub, ha="left", va="top", fontsize=P.PLOT_SUBTITLE_FONT_SIZE,
             color=P.PLOT_INK_MUTED)

    hs = [Patch(facecolor=V.gene_color("CAZyme")), Patch(facecolor=C.OVERVIEW_CGC_COLOR)]
    ls = ["CGC with substrate prediction", "CGC"]
    fig.legend(hs, ls, loc="lower left", bbox_to_anchor=(L_IN / FIG_W, 0.04 / H), ncol=2, frameon=False,
               fontsize=fs, handlelength=1.4, columnspacing=1.6)
    bar = V.nice_step(wmax, 5)
    frac = bar / (xhi - xlo)
    sax = fig.add_axes([(L_IN + ax_w * (1 - frac) - 0.05) / FIG_W, 0.12 / H, ax_w * frac / FIG_W, 0.22 / H])
    sax.set_xlim(0, 1)
    sax.set_ylim(0, 1)
    sax.axis("off")
    sax.plot([0, 1], [0.8, 0.8], color=P.PLOT_INK, lw=1.3, solid_capstyle="butt")
    for xx in (0, 1):
        sax.plot([xx, xx], [0.6, 1.0], color=P.PLOT_INK, lw=1.0)
    sax.text(0.5, 0.45, V.fmt_bp(bar), ha="center", va="top", fontsize=P.PLOT_SMALL_FONT_SIZE, color=P.PLOT_INK)
    return V.save_figure(fig, out_dir, C.CGC_OVERVIEW_FILE_TEMPLATE, {}, formats, dpi)
