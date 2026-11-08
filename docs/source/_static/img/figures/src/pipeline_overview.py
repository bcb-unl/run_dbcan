"""Draw pipeline_overview.png: the four run_dbcan steps, their tools and output files, and the easy_* wrappers.

Writes two variants: pipeline_overview.png (used by the docs) and pipeline_overview_alt.png. Which one is which
follows _style.USE_ACCENT: True -> the docs figure has one muted accent (step numbers, arrows) and the alt is
strictly monochrome; False -> the other way round.

Run from the repository root:  python docs/source/_static/img/figures/src/pipeline_overview.py
"""
import _style as st

STEPS = [
    ("CAZyme_annotation",
     ["Gene calling: pyrodigal (prok, meta)", "DIAMOND: CAZy.dmnd", "pyHMMER: dbCAN.hmm, dbCAN-sub.hmm"],
     ["overview.tsv", "uniInput.faa", "uniInput.gff"]),
    ("gff_process",
     ["DIAMOND: TCDB, SulfAtlas, MEROPS", "TF: DIAMOND, or TF.hmm with --fungi", "pyHMMER: STP.hmm"],
     ["total_cgc_info.tsv", "cgc.gff"]),
    ("cgc_finder",
     ["Joins CAZymes and signature genes with", "at most 2 other genes between them",
      "Default: at least 1 CAZyme and 1 TC"],
     ["cgc_standard_out.tsv", "cgc_standard_out_summary.tsv"]),
    ("substrate_prediction",
     ["DIAMOND: CGCs vs dbCAN-PUL (PUL.dmnd)", "dbCAN-sub: substrates of CAZyme families"],
     ["substrate_prediction.tsv", "synteny_pdf/"]),
]
OTHER = [("run_dbcan database", "downloads the databases (once, before step 1)"),
         ("cgc_circle_plot, dbcan_plot", "figures of the results"),
         ("dbcan_utils", "read abundance (metagenomes)")]
WRAPPERS = [("easy_CGC", 3), ("easy_substrate", 4)]

W = 9.4
X_BOX, BOX_W, BOX_H = 0.0, 2.0, 0.5
X_TOOLS, X_FILES = 2.32, 5.42
X_WRAP = [7.85, 8.85]
LINE = 0.195          # line pitch of 10 pt text
ROW = 0.92            # row pitch


def draw(use_accent):
    acc = st.accent(use_accent)
    y_input = 0.12
    y_head = 0.76
    y0 = y_head + 0.36                     # top of the first step box
    y_end = y0 + 3 * ROW + BOX_H           # bottom of the last step box
    H = y_end + 0.82
    fig, ax = st.canvas(W, H)
    xc = X_BOX + BOX_W / 2

    # input row
    st.text(ax, xc, y_input + 0.06, "Input", size=st.L, weight="bold", ha="center")
    st.text(ax, X_TOOLS, y_input + 0.06,
            "Genome (--mode prok), metagenome contigs (--mode meta) or proteins (--mode protein)")
    st.text(ax, X_TOOLS, y_input + 0.06 + LINE + 0.02,
            "Protein mode: step 1 needs only the FASTA; the CGC steps 2–4 also need its GFF "
            "(--input_gff, --gff_type).", size=st.S, color=st.GREY)
    st.arrow(ax, (xc, y_input + 0.22), (xc, y0), color=acc)

    # column heads
    for x, t, ha, mono in ((X_TOOLS, "Tools and databases", "left", False), (X_FILES, "Output files", "left", False),
                           (X_WRAP[0], "easy_CGC", "center", True), (X_WRAP[1], "easy_substrate", "center", True)):
        st.text(ax, x, y_head, t, size=st.S, color=st.GREY, ha=ha, mono=mono)
    st.hline(ax, X_TOOLS, W, y_head + 0.15, color=st.INK)

    for i, (name, tools, files) in enumerate(STEPS):
        top = y0 + i * ROW
        mid = top + BOX_H / 2
        st.box(ax, X_BOX, top, BOX_W, BOX_H)
        st.text(ax, X_BOX + 0.14, mid, str(i + 1), size=st.L, weight="bold", color=acc)
        st.text(ax, X_BOX + 0.36, mid, name, size=st.L, weight="bold")
        if i:
            st.arrow(ax, (xc, top - ROW + BOX_H), (xc, top), color=acc)
        block_top = mid - (max(len(tools), len(files)) - 1) * LINE / 2
        for k, t in enumerate(tools):
            st.text(ax, X_TOOLS, block_top + k * LINE, t)
        for k, f in enumerate(files):
            st.text(ax, X_FILES, block_top + k * LINE, f, mono=True, color=st.GREY, size=st.S)
        if i < 3:
            st.hline(ax, X_TOOLS, W, top + BOX_H + (ROW - BOX_H) / 2, lw=0.5)

    # wrappers: a line through the steps each one runs
    for x, (_, n) in zip(X_WRAP, WRAPPERS):
        mids = [y0 + i * ROW + BOX_H / 2 for i in range(n)]
        st.vline(ax, x, mids[0], mids[-1], lw=1.2)
        for m in mids:
            ax.plot([x], [m], marker="o", ms=5.2, mfc=st.INK, mec="white", mew=1.0, zorder=3)

    # other commands: name in the tools column, what it does in the files column
    y_foot = y_end + 0.36
    st.hline(ax, X_TOOLS, W, y_foot - 0.2, color=st.INK)
    st.text(ax, xc, y_foot, "Other commands", size=st.L, weight="bold", ha="center")
    for k, (cmd, what) in enumerate(OTHER):
        st.text(ax, X_TOOLS, y_foot + k * LINE, cmd, mono=True)
        st.text(ax, X_FILES, y_foot + k * LINE, what, color=st.GREY)
    return fig


if __name__ == "__main__":
    st.save(draw(st.USE_ACCENT), "pipeline_overview.png")
    st.save(draw(not st.USE_ACCENT), "pipeline_overview_alt.png")
