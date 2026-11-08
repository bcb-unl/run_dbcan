"""Draw protocol_flow.png: the four modules of the metagenomics protocol, the command of each step
and the files it writes (used by docs/source/metagenomics_pipeline/run_from_raw_reads.rst).

Run from the repository root:  python docs/source/_static/img/figures/src/protocol_flow.py
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

from dbcan.constants.base_constants import CGC_FEATURE_COLORS as C

OUT = Path(__file__).resolve().parents[1] / "protocol_flow.png"

# module title, colour, title colour, [(command, [output files])]
MODULES = [
    ("1  Reads to genes  (P1-P4)", C["Other"], "white", [
        ("kraken2, trim_galore", ["S_1_val_1.fq.gz", "S_2_val_2.fq.gz"]),
        ("megahit", ["megahit_S/S.contigs.fa"]),
        ("pyrodigal -p meta", ["S.faa  S.gff"]),
    ]),
    ("2  run_dbcan  (P5-P7)", C["CAZyme"], "white", [
        ("run_dbcan easy_substrate", ["S.dbCAN/", "  overview.tsv", "  cgc.gff", "  cgc_standard_out.tsv",
                                      "  substrate_prediction.tsv", "  synteny_pdf/"]),
    ]),
    ("3  Abundance  (P8-P12)", C["TF"], "white", [
        ("bwa mem, samtools", ["S.bam  S.bam.bai"]),
        ("dbcan_utils cal_coverage", ["S.depth.txt"]),
        ("samtools depth -r", ["S.cgc.depth.txt"]),
        ("dbcan_utils *_abund", ["S_abund/  (7 tables)"]),
    ]),
    ("4  Plots  (P13-P15)", C["STP"], "#333333", [
        ("dbcan_plot heatmap_plot", ["heatmap.pdf"]),
        ("dbcan_plot bar_plot", ["fam/subfam/EC .pdf"]),
        ("dbcan_plot CGC_*_plot", ["CGC map, coverage,", "synteny figures"]),
    ]),
]

COL_W, GAP = 3.25, 0.45
CMD_H, FILE_H, STEP_GAP = 0.34, 0.27, 0.12
TOP = 0.0

def body_height(steps):
    return sum(CMD_H + FILE_H * len(f) + STEP_GAP for _, f in steps) + 0.1

H = max(body_height(s) for *_, s in MODULES)
fig, ax = plt.subplots(figsize=(4 * COL_W + 3 * GAP, H + 1.1))
ax.set_xlim(-0.05, 4 * COL_W + 3 * GAP + 0.05)
ax.set_ylim(-H - 0.7, 0.6)
ax.axis("off")

for i, (title, colour, tcol, steps) in enumerate(MODULES):
    x0 = i * (COL_W + GAP)
    ax.add_patch(FancyBboxPatch((x0, 0.05), COL_W, 0.45, boxstyle="round,pad=0.02,rounding_size=0.08",
                                fc=colour, ec="none"))
    ax.text(x0 + 0.12, 0.275, title, color=tcol, fontsize=11, fontweight="bold", va="center")
    ax.add_patch(FancyBboxPatch((x0, -H), COL_W, H - 0.02, boxstyle="round,pad=0.02,rounding_size=0.08",
                                fc="#F7F7F7", ec=colour, lw=1.3))
    y = -0.1
    for cmd, files in steps:
        ax.text(x0 + 0.12, y - CMD_H / 2, cmd, fontsize=9, family="monospace", fontweight="bold",
                color="#222222", va="center")
        y -= CMD_H
        for k, f in enumerate(files):
            prefix = "→ " if k == 0 else "  "
            ax.text(x0 + 0.2, y - FILE_H / 2, prefix + f, fontsize=8.5, family="monospace", color="#555555",
                    va="center")
            y -= FILE_H
        y -= STEP_GAP
    if i < 3:
        xa = x0 + COL_W + 0.04
        ax.add_patch(FancyArrowPatch((xa, 0.275), (xa + GAP - 0.08, 0.275), arrowstyle="-|>",
                                     mutation_scale=15, color="#444444", lw=1.5))

ax.text(0, -H - 0.45, "S = sample (e.g. Wet2014). Module 3 maps the clean reads of module 1 to the contigs of the "
        "same sample and counts reads per gene of S.dbCAN/cgc.gff.", fontsize=9, color="#555555", va="center")

fig.savefig(OUT, dpi=150, bbox_inches="tight", facecolor="white")
print(OUT)
