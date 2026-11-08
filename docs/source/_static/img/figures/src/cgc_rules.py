"""Draw cgc_rules.png: how cgc_finder joins signature genes into CGCs (default settings),
one gene = one role, and how the labels of one protein are used (gff_process / cgc_finder).

Run from the repository root:  python docs/source/_static/img/figures/src/cgc_rules.py
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

from dbcan.constants.base_constants import CGC_FEATURE_COLORS as COL

OUT = Path(__file__).resolve().parents[1] / "cgc_rules.png"
INK, MUTED, OK, BAD = "#1F2A36", "#566573", "#1E8449", "#C0392B"

fig, ax = plt.subplots(figsize=(10.4, 5.4), dpi=150)
ax.set_xlim(0, 104)
ax.set_ylim(0, 54)
ax.axis("off")


def gene(x, y, kind, label="", w=5.6, h=2.6, fwd=True, stripe=None):
    tip = 1.2
    if fwd:
        pts = [(x, y - h / 2), (x + w - tip, y - h / 2), (x + w, y), (x + w - tip, y + h / 2), (x, y + h / 2)]
    else:
        pts = [(x + tip, y - h / 2), (x + w, y - h / 2), (x + w, y + h / 2), (x + tip, y + h / 2), (x, y)]
    ax.add_patch(Polygon(pts, closed=True, fc=COL[kind], ec="white", lw=0.8))
    if stripe:  # a second label on the same gene
        ax.add_patch(Polygon([(x + w * 0.62, y - h / 2), (x + w * 0.62, y + h / 2),
                              (x + w * 0.80, y + h / 2), (x + w * 0.80, y - h / 2)],
                             closed=True, fc=COL[stripe], ec="white", lw=0.6))
    if label:
        ax.text(x + w / 2, y + h / 2 + 1.0, label, ha="center", va="bottom", fontsize=7.4, color=INK)
    return x + w


def row(x, y, genes, step=6.2):
    ax.plot([x - 0.5, x + step * len(genes) - 0.4], [y, y], color="#BDC3C7", lw=1, zorder=0)
    for k, g in enumerate(genes):
        gene(x + k * step, y, *g)


def bracket(xa, xb, y, text, color):
    ax.plot([xa, xa, xb, xb], [y + 0.8, y, y, y + 0.8], color=color, lw=1.3)
    ax.text((xa + xb) / 2, y - 1.0, text, ha="center", va="top", fontsize=8, color=color, weight="bold")


def title(y, t, sub=""):
    ax.text(1, y, t, fontsize=9.5, color=INK, weight="bold", va="center")
    if sub:
        ax.text(1, y - 2.5, sub, fontsize=7.8, color=MUTED, va="center")


N = ("Other", "")
# ---- A: joining rule
title(52, "A  Joining rule (defaults)",
      "signature genes are joined while at most --num_null_gene 2 other genes lie between them; "
      "a CGC needs >= 1 CAZyme and >= 1 TC (--additional_genes TC)")
ya, x = 41, 3
row(x, ya, [("CAZyme", "GH13"), N, ("TC", "TC"), N, N, ("CAZyme", "GH31"),
            N, N, N, ("CAZyme", "GH2"), ("CAZyme", "GH42"), N, ("TF", "TF")], step=6.2)
s = 6.2
bracket(x, x + 5 * s + 5.6, ya - 3.2, "CGC: GH13 + TC + GH31 (<= 2 null genes apart)", OK)
ax.text(x + 7 * s + 2.8, ya - 4.2, "3 null genes:\ncluster ends", ha="center", va="top", fontsize=7.4, color=MUTED)
bracket(x + 9 * s, x + 10 * s + 5.6, ya - 3.2, "no TC: not a CGC", BAD)
ax.text(x + 12 * s + 2.8, ya - 4.2, "TF counts as a null\ngene by default", ha="center", va="top", fontsize=7.4,
        color=MUTED)

# ---- B: one gene = one role
title(29, "B  One gene = one role",
      "CAZyme and each required extra type come from different genes")
ax.text(1, 24, "(example: --additional_genes STP; STP domain outside the CAZyme domain)", fontsize=7.8,
        color=MUTED, va="center")
yb = 15.5
row(3, yb, [N, ("CAZyme", "GH13 + STP", 5.6, 2.6, True, "STP"), N])
bracket(3 + s, 3 + s + 5.6, yb - 3.2, "one gene only: not a CGC", BAD)
row(30, yb, [("CAZyme", "GH13 + STP", 5.6, 2.6, True, "STP"), ("CAZyme", "GH31"), N])
bracket(30, 30 + s + 5.6, yb - 3.2, "CGC: GH31 = CAZyme, GH13 gene = STP", OK)

# ---- C: how the labels of one protein are used
xc = 62
ax.text(xc, 29, "C  How the labels of one protein are used", fontsize=9.5, color=INK, weight="bold", va="center")
lines = [
    ("CAZyme (#ofTools >= 2) + TC / TF(DIAMOND) / sulfatase hit", "cgc.gff: CAZyme only", "gff_process"),
    ("not a CAZyme, several full-length hits", "cgc.gff: best bit score only", "gff_process"),
    ("STP / fungal TF / peptidase domain on the CAZyme domain", "kept in cgc.gff; counts as CAZyme only",
     "cgc_finder"),
]
for k, (a, b, step) in enumerate(lines):
    y = 24.5 - k * 5.4
    ax.text(xc, y, a, fontsize=7.6, color=INK, va="center")
    ax.text(xc + 2, y - 2.2, "→ " + b, fontsize=7.8, color=OK, va="center", weight="bold")
    ax.text(103.5, y - 2.2, step, fontsize=7.0, color=MUTED, va="center", ha="right", style="italic")
ax.text(xc, 7.5, "all hits stay in total_cgc_info.tsv", fontsize=7.4, color=MUTED, va="center", style="italic")

# legend
lx = 3
for k, kind in enumerate(["CAZyme", "TC", "TF", "STP", "Other"]):
    ax.add_patch(Polygon([(lx, 2.2), (lx + 2.4, 2.2), (lx + 2.4, 3.8), (lx, 3.8)], fc=COL[kind], ec="none"))
    ax.text(lx + 3, 3, "null gene" if kind == "Other" else kind, fontsize=7.6, color=INK, va="center")
    lx += 11

fig.savefig(OUT, dpi=150, bbox_inches="tight", pad_inches=0.08, facecolor="white")
print(OUT, OUT.stat().st_size)
