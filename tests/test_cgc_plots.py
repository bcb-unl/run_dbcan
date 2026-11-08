"""CGC figures: one palette, Gene Type priority colours, readers, layouts and the plot CLIs."""
import re
import shutil
import warnings
from pathlib import Path
from types import SimpleNamespace

import pytest
from click.testing import CliRunner

import dbcan.constants.base_constants as B
import dbcan.constants.cgcfinder_constants as CF
import dbcan.constants.plot_cgc_circle_constants as CC
import dbcan.constants.plots_constants as P
from dbcan.plot import cgc_viz as V

DATA = Path(__file__).parent / "_data" / "cgc_plot"
RUN = DATA / "run"
DB = DATA / "pul_db"
ORIGINAL_COLORS = {
    "CAZyme": "#E67E22", "TC": "#2ECC71", "TF": "#9B59B6", "STP": "#F1C40F",
    "Peptidase": "#16A085", "Sulfatase": "#010E1B", "Other": "#95A5A6",
}


@pytest.fixture
def run_dir(tmp_path):
    dst = tmp_path / "run"
    shutil.copytree(RUN, dst)
    return dst


def plot_cli(*args):
    from dbcan.utils.plots import cli
    return CliRunner().invoke(cli, [str(a) for a in args])


# ---------------------------------------------------------------- palette / priority
def test_palette_is_the_original_single_source():
    assert B.CGC_FEATURE_COLORS == ORIGINAL_COLORS
    assert P.CGC_FEATURE_COLORS is B.CGC_FEATURE_COLORS
    assert CC.CGC_FEATURE_COLORS is B.CGC_FEATURE_COLORS
    assert V.PALETTE is B.CGC_FEATURE_COLORS
    # no other module carries its own copy of the gene-type colours
    pkg = Path(B.__file__).resolve().parents[1]
    offenders = []
    for src in pkg.rglob("*.py"):
        if src.name == "base_constants.py":
            continue
        text = src.read_text()
        if any(re.search(re.escape(c), text, re.I) for c in ORIGINAL_COLORS.values()):
            offenders.append(src.name)
    assert offenders == []


@pytest.mark.parametrize("annotation,expected", [
    ("CAZyme|GH13_31+TC|8.A.9.1.1", "CAZyme"),
    ("TC|2.A.1+CAZyme|GH5", "CAZyme"),
    ("STP|ACT+TF|LacI", "TF"),
    ("TF|LacI+STP|Peripla_BP_3", "TF"),
    ("PEPTIDASE|M38+SULFATLAS|S1_8", "Sulfatase"),
    ("PEPTIDASE|M38", "Peptidase"),
    ("SULFATLAS|S1_19", "Sulfatase"),
    ("null", "Other"),
    ("", "Other"),
])
def test_primary_type_follows_gene_type_priority(annotation, expected):
    assert V.primary_gene_type(annotation) == expected


def test_priority_is_shared_with_cgcfinder():
    from dbcan.annotation.CGCFinder import CGCFinder
    assert CF.PRIORITY is B.CGC_GENE_TYPE_PRIORITY
    for ann in ("CAZyme|GH13_31+TC|8.A.9.1.1", "STP|ACT+TF|LacI", "Peptidase|M38+Sulfatase|S1_8"):
        assert CGCFinder.get_gene_type(ann) == V.primary_gene_type(ann)


def test_multi_role_gene_has_one_colour_and_extra_roles_in_label():
    g = V.Gene.build("p", "c", 1, 900, "+", "CAZyme|GH13_31+TC|8.A.9.1.1")
    assert g.gtype == "CAZyme"
    assert g.label == "GH13_31 + TC 8.A.9.1.1"
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    poly = V.draw_gene(ax, g.start, g.end, g.strand, g.gtype, 0, 1, 50)
    assert poly.get_facecolor()[:3] == pytest.approx(tuple(int(ORIGINAL_COLORS["CAZyme"][i:i + 2], 16) / 255
                                                           for i in (1, 3, 5)))
    assert len(ax.patches) == 1        # no stripes
    plt.close(fig)


def test_types_present_lists_only_drawn_types():
    cgcs = V.read_cgc_table(RUN / "cgc_standard_out.tsv")
    assert V.types_present(cgcs["ctgA|CGC1"]) == ["CAZyme", "TC", "TF", "Other"]


# ---------------------------------------------------------------- ids / readers
def test_split_cgc_id_handles_pipes_in_contig_and_protein():
    assert V.split_cgc_id("ctgB|x|CGC2|b2|CAZyme|GH5_2") == ("ctgB|x", "CGC2", "b2|CAZyme|GH5_2")
    assert V.split_cgc_id("NC_1.1|CGC12") == ("NC_1.1", "CGC12", "")
    assert V.split_cgc_id("no_cgc_here") is None


def test_read_cgc_table_keys_and_substrates_deduplicated():
    cgcs = V.read_cgc_table(RUN / "cgc_standard_out.tsv")
    assert set(cgcs) == {"ctgA|CGC1", "ctgB|x|CGC2"}
    calls = V.read_substrates(RUN / "substrate_prediction.tsv")
    assert calls["ctgA|CGC1"].names() == ["xylan", "arabinan"]
    assert calls["ctgB|x|CGC2"].names() == ["cellulose"]


def test_read_depth_filters_contig_and_window():
    a = V.read_depth(RUN / "two_contigs.depth.txt", "ctgA", 1000, 1100)
    b = V.read_depth(RUN / "two_contigs.depth.txt", "ctgB|x", 1000, 1100)
    assert a and b
    assert min(a) >= 1000 and max(a) <= 1100
    assert set(b.values()) == {500.0}
    assert 500.0 not in a.values()


def test_coverage_panels_use_only_the_cgc_contig():
    from dbcan.plot.cgc_map import load_depth_tracks
    [(name, xs, ys)] = load_depth_tracks([("s", RUN / "two_contigs.depth.txt")], "ctgA", 950, 5250)
    assert ys.max() < 500


# ---------------------------------------------------------------- synteny
def _synteny_inputs(cutoffs=None):
    from dbcan.plot.syntenic_plot import HomologyCutoffs, load_pul_genes, read_pul_links
    cgcs = V.read_cgc_table(RUN / "cgc_standard_out.tsv")
    pul_genes = load_pul_genes(DB, "PULT001")
    links = read_pul_links(RUN / "PUL_blast.out", "ctgA|CGC1", "PULT001", cgcs["ctgA|CGC1"], pul_genes,
                           cutoffs or HomologyCutoffs())
    return cgcs, pul_genes, links


def test_synteny_links_pass_substrate_prediction_cutoffs():
    from dbcan.plot.syntenic_plot import HomologyCutoffs
    _, pul_genes, links = _synteny_inputs()
    assert sorted((l.q.pid, l.s.pid) for l in links) == [("a1", "r3"), ("a2", "r2"), ("a5", "r1")]  # e=0.5 dropped
    _, _, strict = _synteny_inputs(HomologyCutoffs(evalue=1e-2, identity=50.0, coverage=0.0))
    assert sorted(l.q.pid for l in strict) == ["a1", "a2"]
    assert {g.gtype for g in pul_genes} == {"CAZyme", "TC", "Peptidase"}   # PEPTIDASE alias normalised


def test_synteny_links_parse_contig_with_pipe():
    from dbcan.plot.syntenic_plot import HomologyCutoffs, load_pul_genes, read_pul_links
    cgcs = V.read_cgc_table(RUN / "cgc_standard_out.tsv")
    links = read_pul_links(RUN / "PUL_blast.out", "ctgB|x|CGC2", "PULT001", cgcs["ctgB|x|CGC2"],
                           load_pul_genes(DB, "PULT001"), HomologyCutoffs())
    assert [(l.q.pid, l.s.pid) for l in links] == [("b2", "r1")]


def test_revcomp_option_controls_reference_orientation():
    from dbcan.plot.syntenic_plot import reference_transform
    cgcs, pul_genes, links = _synteny_inputs()
    pmap_on, flipped = reference_transform(cgcs["ctgA|CGC1"], pul_genes, links, revcomp=True)
    assert flipped
    r1 = next(g for g in pul_genes if g.pid == "r1")
    assert pmap_on(r1)[2] == "-"            # '+' gene drawn on the reverse strand
    pmap_off, flipped_off = reference_transform(cgcs["ctgA|CGC1"], pul_genes, links, revcomp=False)
    assert not flipped_off and pmap_off(r1)[2] == "+"


def test_substrate_options_carry_synteny_revcomp():
    import click
    from dbcan.configs.base_config import create_config
    from dbcan.configs.cgc_substrate_config import SynPlotConfig
    from dbcan.parameter import cgc_substrate_base_options
    seen = {}

    @click.command()
    @cgc_substrate_base_options
    def cmd(**kw):
        seen.update(kw)

    CliRunner().invoke(cmd, ["--output_dir", "o", "--db_dir", "d"], catch_exceptions=False)
    assert seen["synteny_revcomp"] is True
    CliRunner().invoke(cmd, ["--output_dir", "o", "--db_dir", "d", "--no-synteny-revcomp"], catch_exceptions=False)
    assert seen["synteny_revcomp"] is False
    assert create_config(SynPlotConfig, **seen).synteny_revcomp is False


def test_substrate_prediction_writes_synteny_pdfs(run_dir):
    from dbcan.core import run_dbcan_syn_plot
    from dbcan.configs.cgc_substrate_config import SynPlotConfig
    run_dbcan_syn_plot(SynPlotConfig(db_dir=str(DB), output_dir=str(run_dir), synteny_revcomp=False))
    assert (run_dir / "synteny_pdf" / "ctgA_CGC1_PULT001-syntenic.pdf").stat().st_size > 0


def test_syntenic_plot_import_does_not_touch_global_warning_filters():
    import importlib
    import dbcan.plot.syntenic_plot as S
    before = list(warnings.filters)
    importlib.reload(S)
    assert warnings.filters == before


# ---------------------------------------------------------------- circle / overview
def _single_contig(run_dir):
    for name in ("cgc.gff", "cgc_standard_out.tsv"):
        p = run_dir / name
        lines = p.read_text().splitlines(keepends=True)
        p.write_text("".join(l for l in lines if "ctgB" not in l))


def test_layout_auto_picks_circos_for_one_contig_and_overview_otherwise(run_dir):
    from dbcan.plot.plot_cgc_circle import CGCCircosPlot
    assert CGCCircosPlot(SimpleNamespace(output_dir=str(run_dir), layout="auto")).resolved_layout() == "overview"
    assert CGCCircosPlot(SimpleNamespace(output_dir=str(run_dir), layout="circos")).resolved_layout() == "circos"
    _single_contig(run_dir)
    assert CGCCircosPlot(SimpleNamespace(output_dir=str(run_dir), layout="auto")).resolved_layout() == "circos"


def test_circos_ring_selects_cgc_genes_by_coordinates_not_protein_id(run_dir):
    from dbcan.plot.plot_cgc_circle import CGCCircosPlot
    plot = CGCCircosPlot(SimpleNamespace(output_dir=str(run_dir), layout="auto"))
    ring = plot.cgc_ring_genes("ctgA")
    assert [(g.pid, g.start) for g in ring] == [("a1", 1000), ("a2", 2100), ("a3", 3000), ("NA", 3700), ("a5", 4000)]
    assert all(g.start != 8000 for g in ring)      # the RNA outside the CGC (protein_id=NA) is not on the ring


def test_deg_na_is_not_down_regulated(run_dir):
    from dbcan.plot.plot_cgc_circle import read_deg
    assert read_deg(run_dir / "DEG.tsv") == {"a1": 2.5, "b2": -1.2}


@pytest.mark.parametrize("layout,stem", [("auto", "cgc_overview"), ("circos", "cgc_circos_plot")])
def test_cgc_circle_plot_command_writes_png_and_pdf(run_dir, layout, stem):
    pytest.importorskip("pycirclize")
    from dbcan.main import cli
    res = CliRunner().invoke(cli, ["cgc_circle_plot", "--output_dir", str(run_dir), "--layout", layout])
    assert res.exit_code == 0, res.output
    for ext in ("png", "pdf"):
        assert (run_dir / "cgc_circos" / f"{stem}.{ext}").stat().st_size > 0
    assert not list((run_dir / "cgc_circos").glob("*.svg"))       # no per-contig SVG any more


def test_cgc_circle_plot_single_contig_auto_is_circos(run_dir):
    pytest.importorskip("pycirclize")
    from dbcan.main import cli
    _single_contig(run_dir)
    res = CliRunner().invoke(cli, ["cgc_circle_plot", "--output_dir", str(run_dir)])
    assert res.exit_code == 0, res.output
    assert (run_dir / "cgc_circos" / "cgc_circos_plot.png").exists()


def test_cgc_circle_plot_fails_without_inputs(tmp_path):
    from dbcan.main import cli
    res = CliRunner().invoke(cli, ["cgc_circle_plot", "--output_dir", str(tmp_path)])
    assert res.exit_code != 0


# ---------------------------------------------------------------- dbcan_plot CLI
def test_cgc_plot_writes_into_output_dir(run_dir, tmp_path):
    out = tmp_path / "figs"
    res = plot_cli("CGC_plot", "-i", run_dir, "--cgcid", "ctgA|CGC1", "-o", out, "--format", "png,pdf")
    assert res.exit_code == 0, res.output
    assert (out / "ctgA_CGC1.cgc.png").exists() and (out / "ctgA_CGC1.cgc.pdf").exists()


def test_cgc_coverage_plot_several_samples(run_dir, tmp_path):
    res = plot_cli("CGC_coverage_plot", "-i", run_dir, "--cgcid", "ctgA|CGC1", "-o", tmp_path, "--format", "png",
                   "--reads-count", run_dir / "two_contigs.depth.txt", "--reads-count", run_dir / "sample2.depth.txt",
                   "--samples", "s1,s2")
    assert res.exit_code == 0, res.output
    assert (tmp_path / "ctgA_CGC1.cgc-coverage.png").exists()


@pytest.mark.parametrize("args", [
    ["CGC_synteny_plot"],
    ["CGC_synteny_plot", "--no-synteny-revcomp"],
    ["CGC_synteny_coverage_plot", "--reads-count", RUN / "two_contigs.depth.txt"],
])
def test_synteny_cli(run_dir, tmp_path, args):
    res = plot_cli(*args, "-i", run_dir, "--cgcid", "ctgA|CGC1", "--db-dir", DB, "-o", tmp_path, "--format", "png")
    assert res.exit_code == 0, res.output
    assert len(list(tmp_path.glob("ctgA_CGC1*.png"))) == 1


def test_synteny_cli_pipe_contig_with_explicit_pul(run_dir, tmp_path):
    res = plot_cli("CGC_synteny_plot", "-i", run_dir, "--cgcid", "ctgB|x|CGC2", "--pul", "PULT001",
                   "--db-dir", DB, "-o", tmp_path)
    assert res.exit_code == 0, res.output
    assert (tmp_path / "ctgB_x_CGC2_PULT001-syntenic.pdf").exists()


@pytest.mark.parametrize("args", [
    ["CGC_plot"],                                                     # no --cgcid
    ["CGC_plot", "--cgcid", "ctgA|CGC99"],                            # unknown CGC
    ["CGC_coverage_plot", "--cgcid", "ctgA|CGC1"],                    # no --reads-count
    ["CGC_coverage_plot", "--cgcid", "ctgA|CGC1", "--reads-count", "missing.depth"],
    ["CGC_synteny_plot", "--cgcid", "ctgB|x|CGC2"],                   # no PUL hit
])
def test_cgc_plot_failures_exit_non_zero(run_dir, tmp_path, args):
    res = plot_cli(*args, "-i", run_dir, "-o", tmp_path, "--db-dir", DB)
    assert res.exit_code == 1


def _abund_args():
    return ["-i", f"{DATA / 'abund' / 'A_fam_abund.out'},{DATA / 'abund' / 'B_fam_abund.out'}", "--samples", "A,B"]


@pytest.mark.parametrize("vertical", [False, True])
def test_bar_plot_orientation(tmp_path, monkeypatch, vertical):
    import matplotlib.pyplot as plt
    seen = {}
    real_savefig = plt.savefig

    def spy(path, *a, **k):
        bars = [p for p in plt.gca().patches]
        seen["vertical"] = all(b.get_height() >= b.get_width() for b in bars if b.get_height() > 0)
        seen["horizontal"] = all(b.get_width() >= b.get_height() for b in bars if b.get_width() > 0)
        return real_savefig(path, *a, **k)

    monkeypatch.setattr(plt, "savefig", spy)
    args = ["bar_plot", *_abund_args(), "-o", tmp_path, "--pdf", "fam.pdf"] + (["--vertical-bar"] if vertical else [])
    res = plot_cli(*args)
    assert res.exit_code == 0, res.output
    assert (tmp_path / "fam.pdf").exists()
    assert seen["vertical" if vertical else "horizontal"]


def test_heatmap_honours_pdf_and_output_dir(tmp_path):
    res = plot_cli("heatmap_plot", *_abund_args(), "-o", tmp_path, "--pdf", "fam_heat")
    assert res.exit_code == 0, res.output
    assert (tmp_path / "fam_heat.pdf").exists()


@pytest.mark.parametrize("style", ["seaborn", "seaborn-poster", "ggplot"])
def test_old_seaborn_style_names_work(tmp_path, style):
    res = plot_cli("bar_plot", *_abund_args(), "-o", tmp_path, "--plot-style", style)
    assert res.exit_code == 0, res.output
    assert (tmp_path / "bar_plot.pdf").exists()


def test_bar_plot_missing_file_exits_non_zero(tmp_path):
    res = plot_cli("bar_plot", "-i", tmp_path / "nope.out", "--samples", "A", "-o", tmp_path)
    assert res.exit_code == 1


# ---------------------------------------------------------------- review fixes
def test_plot_fixtures_are_not_git_ignored():
    import subprocess
    root = Path(__file__).resolve().parent.parent
    if subprocess.run(["git", "-C", str(root), "rev-parse"], capture_output=True).returncode != 0:
        pytest.skip("not a git checkout")
    files = [str(p.relative_to(root)) for p in DATA.rglob("*") if p.is_file()]
    res = subprocess.run(["git", "-C", str(root), "check-ignore", "--no-index", *files],
                         capture_output=True, text=True)
    assert res.stdout.strip() == "", f"fixtures ignored by .gitignore: {res.stdout}"


def _boxes(lay):
    return [(i, *lay.box(i), lay.tiers[i]) for i, lab in enumerate(lay.labels) if lab]


@pytest.mark.parametrize("seed", range(12))
def test_label_leaders_never_cross_other_labels(seed):
    import random
    rnd = random.Random(seed)
    words = ["TC 2.A.53.3.12", "Cyanate_lyase", "GH13_31 + TC 8.A.9.1.1", "TC 3.A.1.147.9", "CynR", "GH16_3 + S1_19",
             "TC 3.A.2.1.1", "GT2", "CE4", "HTH_AraC"]
    x, extents = 0, []
    for _ in range(14):
        w = rnd.choice([150, 300, 600, 900, 1500])
        extents.append((x, x + w))
        x += w + rnd.choice([0, 20, 80])
    labels = [rnd.choice(words) for _ in extents]
    anchors = [(a + b) / 2 for a, b in extents]
    lay = V.layout_labels(labels, anchors, -100, x + 100, 8.0, 0.9, P.PLOT_LABEL_FONT_SIZE, extents=extents)
    boxes = _boxes(lay)
    for i, b0, b1, t in boxes:
        lx = lay.anchors[i]
        assert extents[i][0] <= lx <= extents[i][1]              # the leader starts on its own gene
        for j, c0, c1, u in boxes:
            if j == i:
                continue
            if u == t:
                assert b1 <= c0 or c1 <= b0, (labels[i], labels[j])   # no overlap within a tier
            if t > 0 and u < t:
                assert not (c0 < lx < c1), (labels[i], labels[j])    # leader passes no lower label


def test_read_depth_accepts_any_whitespace(tmp_path):
    p = tmp_path / "space.depth.txt"
    p.write_text("".join(l.replace("\t", "  ") for l in (RUN / "two_contigs.depth.txt").open()))
    assert V.read_depth(p, "ctgA", 1000, 1100) == V.read_depth(RUN / "two_contigs.depth.txt", "ctgA", 1000, 1100)


def test_read_depth_unknown_contig_is_an_error():
    with pytest.raises(V.PlotInputError, match="ctgZ"):
        V.read_depth(RUN / "two_contigs.depth.txt", "ctgZ", 1000, 1100)


def test_coverage_plot_contig_mismatch_exits_non_zero(run_dir, tmp_path):
    bad = tmp_path / "renamed.depth.txt"
    bad.write_text((run_dir / "two_contigs.depth.txt").read_text().replace("ctgA", "ctgA_v2"))
    res = plot_cli("CGC_coverage_plot", "-i", run_dir, "--cgcid", "ctgA|CGC1", "-o", tmp_path, "--reads-count", bad)
    assert res.exit_code == 1
    assert "ctgA" in res.output
    assert not list(tmp_path.glob("*.cgc-coverage.*"))


def test_cgc_extent_colours_are_not_gene_type_colours():
    gene_colours = {c.upper() for c in B.CGC_FEATURE_COLORS.values()}
    for c in (CC.CGC_RANGE_COLOR, CC.CGC_RANGE_PREDICTED_COLOR, CC.OVERVIEW_CGC_COLOR, CC.OVERVIEW_CGC_PREDICTED_COLOR):
        assert c.upper() not in gene_colours
    assert CC.CGC_RANGE_COLOR != CC.CGC_RANGE_PREDICTED_COLOR


def test_tf_label_prefers_named_family():
    g = V.Gene.build("p", "c", 1, 100, "+", "TF|0036286+TF|HTH_AraC")
    assert g.label == "HTH_AraC"
    assert V.Gene.build("p", "c", 1, 100, "+", "TF|0036286").label == "TF 0036286"


def test_importing_plot_modules_keeps_the_matplotlib_backend():
    import subprocess
    import sys
    code = ("import matplotlib; matplotlib.use('svg'); import dbcan.plot.cgc_viz, dbcan.utils.plots, "
            "dbcan.plot.plot_cgc_circle, dbcan.plot.syntenic_plot; "
            "assert matplotlib.get_backend().lower() == 'svg', matplotlib.get_backend()")
    res = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True,
                         cwd=Path(__file__).resolve().parent.parent)
    assert res.returncode == 0, res.stderr


def test_heatmap_show_fig_without_display_saves_and_fails(tmp_path, monkeypatch):
    import dbcan.utils.plots as plots
    monkeypatch.setattr(plots, "_interactive_backend", lambda: False)
    res = plot_cli("heatmap_plot", *_abund_args(), "-o", tmp_path, "--show-fig")
    assert res.exit_code == 1
    assert "--show-fig" in res.output
    assert (tmp_path / "heatmap.pdf").exists()


def _circos_figure(run_dir, layout="circos"):
    """Draw the circos and return the matplotlib figure (save_figure captured)."""
    from unittest.mock import patch
    from dbcan.plot.plot_cgc_circle import CGCCircosPlot
    got = {}

    def keep(fig, *a, **k):
        got["fig"] = fig
        return []
    with patch.object(V, "save_figure", side_effect=keep):
        CGCCircosPlot(SimpleNamespace(output_dir=str(run_dir), layout=layout)).plot()
    return got["fig"]


def test_circos_figure_is_sized_to_its_labels(run_dir):
    """Single small contig: no empty band between the subtitle and the rings."""
    pytest.importorskip("pycirclize")
    _single_contig(run_dir)
    fig = _circos_figure(run_dir)
    W, H = fig.get_size_inches()
    polar = [a for a in fig.axes if a.name == "polar"][0]
    box = polar.get_position()
    ring_top_in = (box.y0 + box.height / 2) * H + CC.CGC_CIRCLE_DIAMETER_IN / 2
    subtitle_in = H - 0.47
    gap = subtitle_in - ring_top_in
    assert 0.2 < gap < 1.0, gap                 # was ~1.2 in with the fixed layout
    assert H < W + 0.5


def test_forced_circos_names_contigs_and_labels_ticks_once(run_dir):
    pytest.importorskip("pycirclize")
    fig = _circos_figure(run_dir)
    polar = [a for a in fig.axes if a.name == "polar"][0]
    texts = [t.get_text() for t in polar.texts if t.get_text().strip()]
    assert "ctgA" in texts and "ctgB|x" in texts
    tick_labels = [t for t in texts if t.endswith("kb") or t.endswith("Mb")]
    assert len(tick_labels) == len(set(tick_labels))
    assert any("ticks every" in t.get_text() for t in fig.texts)


def _text_r_range(ax, t, renderer):
    """Radius range (ring units) covered by a drawn polar-axes text: the corners of its
    unrotated box, rotated about the anchor (rotation_mode 'anchor', as pycirclize sets)."""
    import numpy as np
    from matplotlib.transforms import Affine2D
    rot = t.get_rotation()
    t.set_rotation(0)
    bb = t.get_window_extent(renderer)
    t.set_rotation(rot)
    x, y = ax.transData.transform(t.get_position())
    corners = np.array([[bb.x0, bb.y0], [bb.x1, bb.y0], [bb.x0, bb.y1], [bb.x1, bb.y1],
                        [(bb.x0 + bb.x1) / 2, bb.y0], [(bb.x0 + bb.x1) / 2, bb.y1]])
    pts = Affine2D().rotate_deg_around(x, y, rot).transform(corners)
    r = np.hypot(*(pts - ax.transData.transform((0, 0))).T) / (
        np.hypot(*(ax.transData.transform((0, 100)) - ax.transData.transform((0, 0)))) / 100)
    return r.min(), r.max()


@pytest.mark.parametrize("with_deg", [True, False])
def test_circos_text_never_crosses_the_rings(run_dir, with_deg):
    """Contig names stay inside the innermost ring (their outer end is anchored, the text runs
    inwards) and tick labels stay above the CGC ring, with and without the DEG rings."""
    pytest.importorskip("pycirclize")
    if not with_deg:
        (run_dir / "DEG.tsv").unlink()
    long_name = "NZ_JAAXYZ010000001.1_a_very_long_contig_name_from_an_assembler"
    for name in ("cgc.gff", "cgc_standard_out.tsv", "substrate_prediction.tsv"):
        p = run_dir / name
        p.write_text(p.read_text().replace("ctgA", long_name))
    fig = _circos_figure(run_dir)
    ax = [a for a in fig.axes if a.name == "polar"][0]
    renderer = fig.canvas.get_renderer()
    inner_ring = CC.DEG_LOG2FC_RANGE[0] if with_deg else CC.CGC_CAZYME_TRACK_RANGE[0]
    names = [t for t in ax.texts if t.get_text().startswith(("NZ_JAAXYZ", "ctgB"))]
    assert {t.get_text()[:4] for t in names} == {"NZ_J", "ctgB"}
    assert any(t.get_text().endswith("…") for t in names)          # the long name is shortened
    for t in names:
        lo, hi = _text_r_range(ax, t, renderer)
        assert hi <= inner_ring - 1 and lo >= CC.CGC_CONTIG_LABEL_MIN_R - 1, (t.get_text(), lo, hi)
    ticks = [t for t in ax.texts if re.fullmatch(r"[\d.]+ (bp|kb|Mb)", t.get_text())]
    assert ticks
    for t in ticks:
        lo, hi = _text_r_range(ax, t, renderer)
        assert lo >= CC.CGC_RANGE_TRACK_RANGE[1] + 0.3 and hi <= CC.CGC_AXIS_TRACK_RANGE[0], (t.get_text(), lo, hi)


def test_gene_track_helper_draws_the_cgc_map_track_into_any_axis():
    """gene_track() is the CGC map's gene track, reusable by other figures (e.g. expression panels)."""
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon
    from dbcan.plot.cgc_map import gene_track, map_width_in
    cgcs = V.read_cgc_table(RUN / P.CGC_RESULT_FILE)
    cid, genes = max(cgcs.items(), key=lambda kv: len(kv[1]))
    W = map_width_in(len(genes))
    with_bar = gene_track(genes, W)
    track = gene_track(genes, W, scale_bar=False)
    assert track.y_scale is None and with_bar.y_scale is not None
    assert 0 < track.height_in < with_bar.height_in             # no scale bar -> shorter track
    assert track.y_top == with_bar.y_top and track.lay.tiers == with_bar.lay.tiers
    assert track.cstart == min(g.start for g in genes) and track.cend == max(g.end for g in genes)
    fig, ax = plt.subplots()
    try:
        track.draw(ax)
        arrows = [p for p in ax.patches if isinstance(p, Polygon)]
        assert len(arrows) == len(genes)
        assert [p.get_facecolor()[:3] for p in arrows] == [
            tuple(int(V.gene_color(g.gtype)[k:k + 2], 16) / 255 for k in (1, 3, 5))
            for g in sorted(genes, key=lambda g: g.start)]
        assert ax.get_xlim() == (track.xlo, track.xhi) and ax.get_ylim() == (track.y_bottom, track.y_top)
        texts = {t.get_text() for t in ax.texts}
        assert {g.label for g in genes if g.label} <= texts       # same stacked family labels
        assert not any(t == g.pid for t in texts for g in genes)  # no protein ids by default
    finally:
        plt.close(fig)
    with pytest.raises(V.PlotInputError):
        gene_track([], W)


def _hourglass_case():
    import matplotlib
    matplotlib.use("Agg", force=True)
    from dbcan.plot.syntenic_plot import Link

    def gene(pid, s, e):
        return V.Gene(contig="c", start=s, end=e, strand="+", pid=pid, gtype="CAZyme", roles=[("CAZyme", "GH32")])

    q1, q2 = gene("q1", 0, 1000), gene("q2", 1000, 2000)
    s1, s2 = gene("s1", 0, 1000), gene("s2", 1000, 2000)
    # two paralogous genes: best pairs q1-s1 / q2-s2 plus the cross hits
    return [Link(q1, s1, 87.0, 900), Link(q2, s2, 88.0, 950), Link(q1, s2, 50.0, 300), Link(q2, s1, 50.0, 290)]


def test_synteny_draws_every_link_as_a_narrow_ribbon_with_its_identity():
    """Multiple matches are drawn too: every link is a ribbon narrower than its genes, a gene with
    several links gives each its own slot (no two ribbons share an end), and every link gets its
    % identity on it without overlapping tags."""
    import matplotlib.colors as mcolors
    import matplotlib.pyplot as plt
    from dbcan.plot.syntenic_plot import draw_links, link_ports

    links = _hourglass_case()
    pmap = lambda g: (g.start, g.end, g.strand)
    ports = link_ports(links, pmap)
    for l in links:
        qx0, qx1, sx0, sx1 = ports[id(l)]
        assert l.q.start <= qx0 < qx1 <= l.q.end and qx1 - qx0 <= P.SYNTENY_LINK_WIDTH_FRACTION * 1000 + 1e-9
        assert l.s.start <= sx0 < sx1 <= l.s.end
    # the two links leaving q1 sit side by side, ordered by partner (s1 left of s2)
    q1_links = sorted((l for l in links if l.q.pid == "q1"), key=lambda l: ports[id(l)][0])
    assert [l.s.pid for l in q1_links] == ["s1", "s2"]
    assert ports[id(q1_links[0])][1] <= ports[id(q1_links[1])][0]

    fig, ax = plt.subplots(figsize=(6, 3))
    ax.set_xlim(-100, 2100)
    ax.set_ylim(-1, 4)
    drawn = draw_links(ax, links, pmap, mcolors.Normalize(20, 100), 2200 / 6.0)
    fig.canvas.draw()
    boxes = [t.get_window_extent() for t in ax.texts]
    texts = sorted(t.get_text() for t in ax.texts)
    n_ribbons = len(ax.patches)
    plt.close(fig)
    assert len(drawn) == 4 and n_ribbons == 4
    assert texts == ["50%", "50%", "87%", "88%"]
    assert not any(a.overlaps(b) for i, a in enumerate(boxes) for b in boxes[i + 1:])


def test_synteny_single_link_is_centred_and_tagged_on_its_ribbon():
    import matplotlib.colors as mcolors
    import matplotlib.pyplot as plt
    from dbcan.plot.syntenic_plot import Link, draw_links, link_ports

    def gene(pid, s, e):
        return V.Gene(contig="c", start=s, end=e, strand="+", pid=pid, gtype="TC", roles=[("TC", "1.A")])

    link = Link(gene("q", 0, 1000), gene("s", 0, 1000), 93.0, 800)
    pmap = lambda g: (g.start, g.end, g.strand)
    qx0, qx1, sx0, sx1 = link_ports([link], pmap)[id(link)]
    assert (qx0 + qx1) / 2 == pytest.approx(500) and (sx0 + sx1) / 2 == pytest.approx(500)
    fig, ax = plt.subplots(figsize=(6, 3))
    ax.set_xlim(-100, 1100)
    draw_links(ax, [link], pmap, mcolors.Normalize(20, 100), 1200 / 6.0)
    t = ax.texts[0]
    edge = t.get_bbox_patch().get_edgecolor()
    ribbon = ax.patches[0].get_facecolor()
    plt.close(fig)
    # the tag sits ON its ribbon (centre line, half-way) and is outlined in the ribbon's grey
    assert t.get_text() == "93%" and t.get_position()[0] == pytest.approx(500)
    assert tuple(edge[:3]) == pytest.approx(tuple(ribbon[:3]))
