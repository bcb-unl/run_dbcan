"""Draw cgc_rules.png: how cgc_finder joins signature genes into CGCs (default settings), that each gene counts
once, by its final Gene Type, and how gff_process chooses the labels.

Same style as pipeline_overview.png (see _style.py). The gene arrows use the gene-type colours of the run_dbcan
plots (CGC_FEATURE_COLORS); everything else is neutral.

Run from the repository root:  python docs/source/_static/img/figures/src/cgc_rules.py
"""
from matplotlib.patches import Polygon, Rectangle

import _style as st
from dbcan.constants.base_constants import CGC_FEATURE_COLORS as COL

W = 8.8
GENE_W, GENE_H, TIP = 0.44, 0.19, 0.12
MONO = {"mono": True, "color": st.GREY, "size": st.S}


def gene(ax, x, y, kind, label="", stripe=None):
    """A forward gene arrow with its left end at x and its centre line at y."""
    h = GENE_H / 2
    ax.add_patch(Polygon([(x, y - h), (x + GENE_W - TIP, y - h), (x + GENE_W, y), (x + GENE_W - TIP, y + h),
                          (x, y + h)], closed=True, fc=COL[kind], ec="none", zorder=2))
    if stripe:  # a second hit on the same gene (a domain of another type)
        ax.add_patch(Rectangle((x + 0.2, y - h), 0.075, GENE_H, fc=COL[stripe], ec="none", zorder=3))
    if label:
        st.text(ax, x + GENE_W / 2, y - h - 0.12, label, size=st.S, ha="center")


def gene_row(ax, x, y, genes, step):
    st.hline(ax, x - 0.08, x + step * (len(genes) - 1) + GENE_W + 0.08, y, color=st.RULE, zorder=1)
    for k, g in enumerate(genes):
        gene(ax, x + k * step, y, *g)


def bracket(ax, xa, xb, y, label, ok=True):
    """Bracket under genes: solid ink = a CGC, dashed grey = rejected candidate."""
    color = st.INK if ok else st.GREY
    ax.plot([xa, xa, xb, xb], [y - 0.06, y, y, y - 0.06], color=color, lw=st.LW,
            linestyle="-" if ok else (0, (3, 2)), solid_joinstyle="miter")
    st.text(ax, (xa + xb) / 2, y + 0.17, label, size=st.S, color=color, ha="center",
            weight="bold" if ok else "normal")


def panel_title(ax, x, y, letter, title):
    st.text(ax, x, y, letter, size=st.L, weight="bold")
    st.text(ax, x + 0.24, y, title, size=st.L, weight="bold")


def legend(ax, x_right, y):
    items = [("CAZyme", "CAZyme"), ("TC", "TC"), ("TF", "TF"), ("STP", "STP"), ("Other", "Null gene")]
    renderer = ax.figure.canvas.get_renderer()
    widths = []
    for _, name in items:  # measure, then place right-aligned
        t = st.text(ax, 0, 0, name, size=st.S)
        widths.append(t.get_window_extent(renderer).width / ax.figure.dpi)
        t.remove()
    sw, gap = 0.22, 0.24
    x = x_right - sum(sw + 0.07 + w for w in widths) - gap * (len(items) - 1)
    for (kind, name), w in zip(items, widths):
        ax.add_patch(Rectangle((x, y - 0.06), sw, 0.12, fc=COL[kind], ec="none"))
        st.text(ax, x + sw + 0.07, y, name, size=st.S)
        x += sw + 0.07 + w + gap


def draw():
    H = 5.1
    fig, ax = st.canvas(W, H)
    N = ("Other",)
    note = {"size": st.S, "color": st.GREY, "ha": "center"}   # notes under a gene row, as bracket()

    # ---- A: joining rule
    panel_title(ax, 0, 0.1, "A", "Joining rule (defaults)")
    legend(ax, W, 0.1)
    st.text(ax, 0.24, 0.4, "Signature genes are joined while at most 2 other genes lie between them "
            "(--num_null_gene 2).")
    st.text(ax, 0.24, 0.4 + 0.2, "A CGC needs at least 1 CAZyme and 1 TC (--additional_genes TC). "
            "A TF gene counts as a null gene by default.")
    step, x0, ya = 0.62, 0.32, 1.18
    gene_row(ax, x0, ya, [("CAZyme", "GH13"), N, ("TC", "TC"), N, N, ("CAZyme", "GH31"), N, N, N,
                          ("CAZyme", "GH2"), ("CAZyme", "GH42"), N, ("TF", "TF")], step)
    gx = lambda k: x0 + k * step  # noqa: E731
    yb = ya + 0.22
    bracket(ax, gx(0), gx(5) + GENE_W, yb, "CGC: GH13 + TC + GH31")
    st.text(ax, (gx(6) + gx(8) + GENE_W) / 2, yb + 0.17, "3 null genes: cluster ends", **note)
    bracket(ax, gx(9), gx(10) + GENE_W, yb, "Not a CGC: no TC", ok=False)
    st.text(ax, gx(12) + GENE_W / 2, yb + 0.17, "TF: null gene", **note)

    y_rule = 1.94
    st.hline(ax, 0, W, y_rule, lw=0.5)

    # ---- B: each gene counts once
    yB = y_rule + 0.28
    panel_title(ax, 0, yB, "B", "Each gene counts once, by its final Gene Type")
    st.text(ax, 0.24, yB + 0.3, "With --additional_genes STP, a GH13 gene with an STP domain")
    st.text(ax, 0.24, yB + 0.5, "has Gene Type CAZyme: it does not fill the STP role.")
    stepB, xB = 0.9, 0.32
    hit = ("CAZyme", "GH13 + STP", "STP")
    for k, (genes, (ia, ib), label, ok) in enumerate([
        ([hit, ("CAZyme", "GH31"), N], (0, 1), "Not a CGC: no STP gene", False),
        ([hit, ("CAZyme", "GH31"), N, ("STP", "STP")], (0, 3), "CGC: GH13 + GH31 + STP", True),
    ]):
        y = yB + 1.1 + k * 0.98
        gene_row(ax, xB, y, genes, stepB)
        bracket(ax, xB + ia * stepB, xB + ib * stepB + GENE_W, y + 0.22, label, ok=ok)

    # ---- C: how labels are chosen
    xC = 4.75
    st.vline(ax, xC - 0.3, y_rule + 0.18, H - 0.05, color=st.RULE, lw=0.5)
    panel_title(ax, xC, yB, "C", "How labels are chosen (gff_process)")
    xr = xC + 2.55              # second column of the table
    t0 = yB + 0.36
    st.hline(ax, xC, W, t0, color=st.INK)
    st.text(ax, xC, t0 + 0.17, "Hits on one gene", size=st.S, color=st.GREY)
    st.runs(ax, xr, t0 + 0.17, [("Label in ", {}), ("cgc.gff", {"mono": True})], size=st.S, color=st.GREY)
    st.hline(ax, xC, W, t0 + 0.34, lw=0.5)
    rows = [(["CAZyme (#ofTools ≥ 2) and a TC,", "TF (DIAMOND) or sulfatase hit"], "CAZyme only"),
            (["Not a CAZyme; several full-length", "hits (TC, TF DIAMOND, sulfatase)"], "Best bit score only")]
    y = t0 + 0.53
    for left, right in rows:
        for k, s in enumerate(left):
            st.text(ax, xC, y + k * 0.2, s)
        st.text(ax, xr, y, right, weight="bold")
        y += 0.2 * len(left) + 0.16
    st.hline(ax, xC, W, y - 0.06, color=st.INK)
    y += 0.22
    st.text(ax, xC, y, "Gene Type column: one type per gene, by priority")
    st.text(ax, xC, y + 0.21, "CAZyme > TC > TF > STP > Sulfatase > Peptidase", weight="bold")
    st.runs(ax, xC, y + 0.52, [("All hits stay in ", {}), ("total_cgc_info.tsv", MONO), (".", {})],
            size=st.S, color=st.GREY)
    return fig


if __name__ == "__main__":
    st.save(draw(), "cgc_rules.png")
