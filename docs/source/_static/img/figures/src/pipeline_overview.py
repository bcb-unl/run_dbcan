"""Draw pipeline_overview.png: the four run_dbcan steps, their inputs/outputs and the easy_* wrappers.

Run from the repository root:  python docs/source/_static/img/figures/src/pipeline_overview.py
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

from dbcan.constants.base_constants import CGC_FEATURE_COLORS as COL

OUT = Path(__file__).resolve().parents[1] / "pipeline_overview.png"
SLATE, GREY, INK = "#34495E", "#ABB2B9", "#1F2A36"

STEPS = [
    ("1  CAZyme_annotation", COL["CAZyme"],
     ["pyrodigal (prok / meta)", "DIAMOND: CAZy.dmnd", "pyHMMER: dbCAN.hmm,", "   dbCAN-sub.hmm"],
     ["overview.tsv", "uniInput.faa / .gff"]),
    ("2  gff_process", COL["TC"],
     ["DIAMOND: TCDB, TF, SulfAtlas,", "   MEROPS", "pyHMMER: STP.hmm", "   (TF.hmm with --fungi)"],
     ["total_cgc_info.tsv", "cgc.gff"]),
    ("3  cgc_finder", SLATE,
     ["CAZyme + signature genes", "within --num_null_gene 2", "(default: CAZyme + TC)"],
     ["cgc_standard_out.tsv", "cgc_standard_out_summary.tsv"]),
    ("4  substrate_prediction", COL["TF"],
     ["DIAMOND vs dbCAN-PUL", "   (PUL.dmnd, dbCAN-PUL.xlsx)", "dbCAN-sub substrate scores"],
     ["substrate_prediction.tsv", "synteny_pdf/"]),
]

fig, ax = plt.subplots(figsize=(10.4, 4.6), dpi=150)
ax.set_xlim(0, 104)  # data units
ax.set_ylim(0, 46)
ax.axis("off")

# inputs
ax.text(1, 44, "Input", fontsize=9, color=INK, weight="bold", va="center")
ax.text(9, 44, "genome / metagenome DNA (--mode prok | meta)   or   proteins (--mode protein)",
        fontsize=8.5, color=INK, va="center")
ax.text(9, 40.8, "protein mode: CAZyme annotation needs only the FASTA; the CGC steps 2-4 also need its GFF "
        "(--input_gff, --gff_type)", fontsize=8.0,
        color="#566573", va="center")

bw, gap, x0, ytop, h = 23.5, 2.2, 1, 37, 24
for i, (title, colr, tools, outs) in enumerate(STEPS):
    x = x0 + i * (bw + gap)
    ax.add_patch(FancyBboxPatch((x, ytop - h), bw, h, boxstyle="round,pad=0,rounding_size=1.2",
                                fc="white", ec=colr, lw=1.6))
    ax.add_patch(FancyBboxPatch((x, ytop - 4.2), bw, 4.2, boxstyle="round,pad=0,rounding_size=1.2",
                                fc=colr, ec=colr, lw=1.6))
    ax.text(x + bw / 2, ytop - 2.1, title, ha="center", va="center", fontsize=9, color="white", weight="bold")
    for k, t in enumerate(tools):
        ax.text(x + 1.2, ytop - 6.6 - k * 2.4, t, fontsize=7.6, color=INK, va="center")
    ax.plot([x + 1, x + bw - 1], [ytop - 16.2, ytop - 16.2], color=GREY, lw=0.6)
    ax.text(x + 1.2, ytop - 17.8, "writes", fontsize=7, color="#7F8C8D", va="center", style="italic")
    for k, t in enumerate(outs):
        ax.text(x + 1.2, ytop - 20 - k * 2.4, t, fontsize=7.1, color=INK, va="center", family="monospace")
    if i:
        ax.add_patch(FancyArrowPatch((x - gap + 0.1, ytop - h / 2), (x - 0.1, ytop - h / 2),
                                     arrowstyle="-|>", mutation_scale=11, color=INK, lw=1.2))

# wrappers
def bracket(i0, i1, y, label):
    xa = x0 + i0 * (bw + gap)
    xb = x0 + i1 * (bw + gap) + bw
    ax.plot([xa, xa, xb, xb], [y + 1, y, y, y + 1], color=INK, lw=1)
    ax.text((xa + xb) / 2, y - 1.4, label, ha="center", va="center", fontsize=8.5, color=INK)

bracket(0, 2, 10.5, "run_dbcan easy_CGC  = steps 1-3")
bracket(0, 3, 5.0, "run_dbcan easy_substrate  = steps 1-4")

ax.text(1, 1, "Afterwards: figures (cgc_circle_plot, dbcan_plot)  ·  read abundance (dbcan_utils)  ·  "
        "databases: run_dbcan database --db_dir db", fontsize=7.8, color="#566573", va="center")

fig.savefig(OUT, dpi=150, bbox_inches="tight", pad_inches=0.08, facecolor="white")
print(OUT, OUT.stat().st_size)
