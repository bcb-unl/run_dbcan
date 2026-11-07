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
DB = DATA / "db"
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
