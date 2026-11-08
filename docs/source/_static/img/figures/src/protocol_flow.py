"""Draw protocol_flow.png: the four modules of the metagenomics protocol (2 x 2 grid), the command of each
step and the files it writes (used by docs/source/metagenomics_pipeline/run_from_raw_reads.rst).

Same style as pipeline_overview.png (see _style.py). Commands and file names are monospace; one file per line.

Run from the repository root:  python docs/source/_static/img/figures/src/protocol_flow.py
"""
import _style as st

# module title, protocol steps, [(command, [output files])]; an output starting with "~" is prose, not a file
MODULES = [
    ("Reads to genes", "P1–P4", [
        ("kraken2, trim_galore", ["S_1_val_1.fq.gz", "S_2_val_2.fq.gz"]),
        ("megahit", ["megahit_S/S.contigs.fa"]),
        ("pyrodigal -p meta", ["S.faa", "S.gff"]),
    ]),
    ("run_dbcan", "P5–P7", [
        ("run_dbcan easy_substrate", ["S.dbCAN/overview.tsv", "S.dbCAN/cgc.gff", "S.dbCAN/cgc_standard_out.tsv",
                                      "S.dbCAN/substrate_prediction.tsv", "S.dbCAN/synteny_pdf/"]),
    ]),
    ("Abundance", "P8–P12", [
        ("bwa mem | samtools sort", ["S.bam", "S.bam.bai"]),
        ("dbcan_utils cal_coverage", ["S.depth.txt"]),
        ("samtools depth -r", ["S.cgc.depth.txt"]),
        ("dbcan_utils *_abund", ["S_abund/", "~7 abundance tables"]),
    ]),
    ("Plots", "P13–P15", [
        ("dbcan_plot heatmap_plot", ["substrate_heatmap.pdf"]),
        ("dbcan_plot bar_plot", ["fam_abund.pdf", "subfam_abund.pdf", "EC_abund.pdf"]),
        ("dbcan_plot CGC_*_plot", ["~CGC map, coverage and synteny PDFs"]),
    ]),
]

W = 6.9
COL_GAP, ROW_GAP = 0.42, 0.5
COL_W = (W - COL_GAP) / 2
HEAD_H = 0.36             # module title band (inside the frame, above a light rule)
CMD, FILE, GAP = 0.19, 0.165, 0.09
PAD = 0.14


def body_height(steps):
    return sum(CMD + FILE * len(files) + GAP for _, files in steps) - GAP + 2 * PAD


def draw():
    acc = st.accent()
    row_h = [HEAD_H + max(body_height(MODULES[i][2]), body_height(MODULES[i + 1][2])) for i in (0, 2)]
    tops = [0.0, row_h[0] + ROW_GAP]
    H = tops[1] + row_h[1] + 0.62
    fig, ax = st.canvas(W, H)

    frames = []
    for i, (title, psteps, steps) in enumerate(MODULES):
        r, c = divmod(i, 2)
        x0, y0, h = c * (COL_W + COL_GAP), tops[r], row_h[r]
        frames.append((x0, y0, h))
        st.box(ax, x0, y0, COL_W, h)
        st.text(ax, x0 + PAD, y0 + HEAD_H / 2, str(i + 1), size=st.L, weight="bold", color=acc)
        st.text(ax, x0 + PAD + 0.2, y0 + HEAD_H / 2, title, size=st.L, weight="bold")
        st.text(ax, x0 + COL_W - PAD, y0 + HEAD_H / 2, psteps, size=st.S, color=st.GREY, ha="right")
        st.hline(ax, x0 + PAD, x0 + COL_W - PAD, y0 + HEAD_H, lw=0.5)
        y = y0 + HEAD_H + PAD
        for cmd, files in steps:
            st.text(ax, x0 + PAD, y + CMD / 2, cmd, mono=True)
            y += CMD
            for f in files:
                prose = f.startswith("~")
                st.text(ax, x0 + PAD + 0.2, y + FILE / 2, f.lstrip("~"), size=st.S, color=st.GREY,
                        mono=not prose)
                y += FILE
            y += GAP

    # 1 -> 2 and 3 -> 4 at the title band; 2 -> 3 as an elbow through the gap between the rows
    for a, b in ((0, 1), (2, 3)):
        xa, ya, _ = frames[a]
        xb, _, _ = frames[b]
        st.arrow(ax, (xa + COL_W + 0.06, ya + HEAD_H / 2), (xb - 0.06, ya + HEAD_H / 2), color=acc, head=8)
    x2, y2, h2 = frames[1]
    x3, y3, _ = frames[2]
    ym = y2 + h2 + ROW_GAP / 2
    xa, xb = x2 + COL_W / 2, x3 + COL_W / 2
    ax.plot([xa, xa, xb], [y2 + h2, ym, ym], color=acc, lw=st.LW, solid_joinstyle="miter")
    st.arrow(ax, (xb, ym), (xb, y3), color=acc, head=8)

    st.text(ax, 0, H - 0.3, "S = sample (e.g. Wet2014). Module 3 maps the clean reads of module 1 to the contigs "
            "of the same sample", size=st.S, color=st.GREY)
    st.runs(ax, 0, H - 0.3 + FILE, [("and counts the reads per gene of ", {}), ("S.dbCAN/cgc.gff", {"mono": True}),
                                     (".", {})], size=st.S, color=st.GREY)
    return fig


if __name__ == "__main__":
    st.save(draw(), "protocol_flow.png")
