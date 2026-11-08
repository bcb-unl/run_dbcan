"""Draw protocol_flow.png: the four modules of the metagenomics protocol (2 x 2 grid), the command of each
step and the files it writes (used by docs/source/metagenomics_pipeline/run_from_raw_reads.rst).

All title bars use one neutral dark grey: the coloured CGC palette (orange = CAZyme, purple = TF, ...) carries
gene-type meaning in the CGC plots, so it is not reused for the modules here.

Run from the repository root:  python docs/source/_static/img/figures/src/protocol_flow.py
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

from dbcan.constants.base_constants import CGC_FEATURE_COLORS as C

OUT = Path(__file__).resolve().parents[1] / "protocol_flow.png"

TITLE_BG = "#4D5656"      # darker shade of the palette's neutral grey (C["Other"]); white text on it is AA
FRAME = C["Other"]
TEXT, FILE_TEXT = "#1B1B1B", "#3D3D3D"

# module title, [(command, [output files])]
MODULES = [
    ("1  Reads to genes  (P1-P4)", [
        ("kraken2, trim_galore", ["S_1_val_1.fq.gz", "S_2_val_2.fq.gz"]),
        ("megahit", ["megahit_S/S.contigs.fa"]),
        ("pyrodigal -p meta", ["S.faa, S.gff"]),
    ]),
    ("2  run_dbcan  (P5-P7)", [
        ("run_dbcan easy_substrate", ["S.dbCAN/overview.tsv", "S.dbCAN/cgc.gff",
                                      "S.dbCAN/cgc_standard_out.tsv",
                                      "S.dbCAN/substrate_prediction.tsv", "S.dbCAN/synteny_pdf/"]),
    ]),
    ("3  Abundance  (P8-P12)", [
        ("bwa mem | samtools sort", ["S.bam, S.bam.bai"]),
        ("dbcan_utils cal_coverage", ["S.depth.txt"]),
        ("samtools depth -r", ["S.cgc.depth.txt"]),
        ("dbcan_utils *_abund", ["S_abund/ (7 tables)"]),
    ]),
    ("4  Plots  (P13-P15)", [
        ("dbcan_plot heatmap_plot", ["substrate_heatmap.pdf"]),
        ("dbcan_plot bar_plot", ["fam_abund.pdf, subfam_abund.pdf,", "EC_abund.pdf"]),
        ("dbcan_plot CGC_*_plot", ["CGC map, coverage and", "synteny PDFs"]),
    ]),
]

COL_W, COL_GAP, ROW_GAP = 4.2, 0.7, 0.55
TITLE_H, CMD_H, FILE_H, STEP_GAP = 0.42, 0.36, 0.29, 0.12


def body_height(steps):
    return sum(CMD_H + FILE_H * len(f) + STEP_GAP for _, f in steps) + 0.12


ROW_H = [max(body_height(MODULES[i][1]), body_height(MODULES[i + 1][1])) for i in (0, 2)]
W = 2 * COL_W + COL_GAP
H = 2 * TITLE_H + sum(ROW_H) + ROW_GAP + 0.75

fig, ax = plt.subplots(figsize=(W + 0.1, H))
ax.set_xlim(-0.05, W + 0.05)
ax.set_ylim(-H + 0.05, 0.05)
ax.axis("off")

anchors = {}
for i, (title, steps) in enumerate(MODULES):
    row, col = divmod(i, 2)
    x0 = col * (COL_W + COL_GAP)
    y_top = 0 if row == 0 else -(TITLE_H + ROW_H[0] + ROW_GAP)
    bh = ROW_H[row]
    ax.add_patch(FancyBboxPatch((x0, y_top - TITLE_H), COL_W, TITLE_H,
                                boxstyle="round,pad=0.0,rounding_size=0.08", fc=TITLE_BG, ec="none"))
    ax.text(x0 + 0.15, y_top - TITLE_H / 2, title, color="white", fontsize=12, fontweight="bold", va="center")
    ax.add_patch(FancyBboxPatch((x0, y_top - TITLE_H - bh), COL_W, bh,
                                boxstyle="round,pad=0.0,rounding_size=0.08", fc="#F7F7F7", ec=FRAME, lw=1.2))
    y = y_top - TITLE_H - 0.1
    for cmd, files in steps:
        ax.text(x0 + 0.15, y - CMD_H / 2, cmd, fontsize=10.5, family="monospace", fontweight="bold",
                color=TEXT, va="center")
        y -= CMD_H
        for k, f in enumerate(files):
            ax.text(x0 + 0.25, y - FILE_H / 2, ("→ " if k == 0 else "  ") + f, fontsize=10,
                    family="monospace", color=FILE_TEXT, va="center")
            y -= FILE_H
        y -= STEP_GAP
    anchors[i] = (x0, y_top, bh)

arrow = dict(arrowstyle="-|>", mutation_scale=18, color="#444444", lw=1.6)
# 1 -> 2 and 3 -> 4: horizontal, at the title bars
for a, b in ((0, 1), (2, 3)):
    xa, ya, _ = anchors[a]
    xb, _, _ = anchors[b]
    ax.add_patch(FancyArrowPatch((xa + COL_W + 0.06, ya - TITLE_H / 2), (xb - 0.06, ya - TITLE_H / 2), **arrow))
# 2 -> 3: from the bottom of module 2 to the title of module 3
x2, y2, bh2 = anchors[1]
x3, y3, _ = anchors[2]
ax.add_patch(FancyArrowPatch((x2 + COL_W / 2, y2 - TITLE_H - bh2 - 0.04), (x3 + COL_W / 2, y3 + 0.04),
                             connectionstyle="arc3,rad=-0.12", **arrow))

ax.text(0, -H + 0.3, "S = sample (e.g. Wet2014). Module 3 maps the clean reads of module 1 to the contigs of the same\n"
        "sample and counts reads per gene of S.dbCAN/cgc.gff.", fontsize=10, color=FILE_TEXT, va="center")

fig.savefig(OUT, dpi=150, bbox_inches="tight", facecolor="white")
print(OUT)
