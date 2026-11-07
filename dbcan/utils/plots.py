#!/usr/bin/env python3
"""``dbcan_plot``: CGC figures (map, coverage, synteny) and abundance charts (bar, heatmap)."""
import logging
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib as mpl  # noqa: E402
import matplotlib.colors as colors  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import rich_click as click  # noqa: E402

import dbcan.constants.cgc_substrate_prediction_constants as CSP  # noqa: E402
import dbcan.constants.plots_constants as P  # noqa: E402
from dbcan.configs.plots_config import PlotsConfig  # noqa: E402
from dbcan.parameter import logging_options  # noqa: E402
from dbcan.main import setup_logging  # noqa: E402
from dbcan.plot import cgc_viz as V  # noqa: E402

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42

logger = logging.getLogger(__name__)


class PlotError(click.ClickException):
    """Shown as 'Error: ...' and makes dbcan_plot exit with status 1."""


# ---------------------------------------------------------------------------
# CGC figures
# ---------------------------------------------------------------------------
def _load_cgc(cfg: PlotsConfig):
    if not cfg.cgcid:
        raise PlotError("--cgcid is required for CGC plots (example: 'contig1|CGC1').")
    run = Path(cfg.input_dir)
    cgcs = V.read_cgc_table(run / P.CGC_RESULT_FILE)
    if cfg.cgcid not in cgcs:
        examples = ", ".join(list(cgcs)[:3])
        raise PlotError(f"CGC '{cfg.cgcid}' not found in {run / P.CGC_RESULT_FILE} (ids look like: {examples}).")
    return run, cgcs


def _depth_files(cfg: PlotsConfig):
    files = list(cfg.reads_count or [])
    if not files:
        raise PlotError("--reads-count is required (samtools depth output; repeat it for several samples).")
    names = list(cfg.samples or [])
    if names and len(names) != len(files):
        raise PlotError(f"--samples has {len(names)} names for {len(files)} --reads-count files.")
    if not names:
        names = [Path(f).name.split(".")[0] for f in files] if len(files) > 1 else [""]
    return list(zip(names, files))


def _substrates(run: Path, cgcid: str):
    call = V.read_substrates(run / P.CGC_SUB_PREDICTION_FILE).get(cgcid)
    return call.names() if call else None


def CGC_plot(cfg: PlotsConfig):
    from dbcan.plot.cgc_map import plot_cgc_map
    run, cgcs = _load_cgc(cfg)
    return plot_cgc_map(cgcs[cfg.cgcid], cfg.cgcid, Path(cfg.output_dir), substrates=_substrates(run, cfg.cgcid),
                        show_protein_ids=cfg.show_protein_ids, formats=cfg.formats, dpi=cfg.dpi)


def CGC_plot_reads_count(cfg: PlotsConfig):
    from dbcan.plot.cgc_map import plot_cgc_map
    run, cgcs = _load_cgc(cfg)
    return plot_cgc_map(cgcs[cfg.cgcid], cfg.cgcid, Path(cfg.output_dir), substrates=_substrates(run, cfg.cgcid),
                        depth_files=_depth_files(cfg), show_protein_ids=cfg.show_protein_ids,
                        formats=cfg.formats, dpi=cfg.dpi)


def _synteny(cfg: PlotsConfig, with_depth: bool):
    from dbcan.plot.syntenic_plot import HomologyCutoffs, draw_cgc_vs_pul
    run, cgcs = _load_cgc(cfg)
    cutoffs = HomologyCutoffs(cfg.evalue_cutoff, cfg.identity_cutoff, cfg.coverage_cutoff)
    return draw_cgc_vs_pul(run, Path(cfg.db_dir), cfg.cgcid, Path(cfg.output_dir), pul=cfg.pul, cutoffs=cutoffs,
                           revcomp=cfg.synteny_revcomp, depth_files=_depth_files(cfg) if with_depth else None,
                           formats=cfg.formats, dpi=cfg.dpi, cgc_table=cgcs,
                           template=P.CGC_SYNTENY_COVERAGE_FILE_TEMPLATE if with_depth else P.SYNTENY_FILE_TEMPLATE)


def CGC_syntenic_with_PUL(cfg: PlotsConfig):
    return _synteny(cfg, with_depth=False)


def CGC_syntenic_with_PUL_abund(cfg: PlotsConfig):
    return _synteny(cfg, with_depth=True)


# ---------------------------------------------------------------------------
# abundance charts
# ---------------------------------------------------------------------------
def combined_datafram_based_on_first_col(pd_lists, samples):
    if len(pd_lists) <= 1:
        return pd_lists[0]
    col_name = pd_lists[0].columns
    on_merge_col = col_name[0]
    merged_table = pd.merge(pd_lists[0], pd_lists[1], on=[on_merge_col], how="outer")
    for i in range(len(pd_lists)):
        ori_names = pd_lists[i].columns
        mod_names = [ori_names[0]] + [ori_names[j] + "_" + samples[i] for j in range(1, len(ori_names))]
        pd_lists[i].columns = mod_names
    for i in range(2, len(pd_lists)):
        merged_table = pd.merge(merged_table, pd_lists[i], on=[on_merge_col], how="outer")

    abundance_col = col_name[1]
    merged_table = merged_table.fillna(0)
    merged_table["diff_abs"] = np.abs(merged_table[abundance_col + "_x"] - merged_table[abundance_col + "_y"])
    merged_table["diff"] = merged_table[abundance_col + "_x"] - merged_table[abundance_col + "_y"]

    rename_columns = []
    abund_index = 0
    for column in merged_table.columns:
        if abundance_col in column and column != on_merge_col:
            rename_columns.append(samples[abund_index])
            abund_index += 1
        else:
            rename_columns.append(column)
    merged_table.columns = rename_columns
    merged_table = merged_table.sort_values("diff_abs", ascending=False)
    return merged_table, on_merge_col


def filter_out_enzyme_number(table):
    bools = []
    for i in table.iloc[:, 0]:
        i = str(i)
        bools.append(not (i in ("PL0", "GH0", "GT0", "CBM0", "AA0", "CE0")))
    return table[bools]


def add_column_type(table):
    table["fam"] = [re.sub(r"[0-9]+", "", str(i)) for i in table["CAZy"]]
    return table


def check_input_files_min_lines(input_files, min_lines=2):
    for f in input_files.split(","):
        if not os.path.isfile(f):
            raise PlotError(f"Abundance file not found: {f}")
        with open(f) as fh:
            lines = sum(1 for _ in fh)
        if lines < min_lines:
            raise PlotError(f"Abundance file {f} has fewer than {min_lines} lines; nothing to plot.")
    return True


def resolve_plot_style(style: str) -> str:
    """Map pre-3.6 seaborn style names to the 'seaborn-v0_8*' names and check availability."""
    style = P.PLOT_STYLE_ALIASES.get(style, style)
    if style not in plt.style.available and style != "default":
        raise PlotError(f"Matplotlib style '{style}' is not available (choose from {', '.join(P.PLOT_STYLE_CHOICES)}).")
    return style


def _read_abundance(input_files: str, samples: str, plot_style: str):
    tables = [filter_out_enzyme_number(pd.read_csv(f, sep="\t")) for f in input_files.split(",")]
    samples_list = samples.split(",")
    if len(tables) != len(samples_list):
        raise PlotError(f"--samples has {len(samples_list)} names for {len(tables)} input files.")
    plt.style.use(resolve_plot_style(plot_style))
    for t, s in zip(tables, samples_list):
        t["sample"] = s
    return tables, samples_list


def _select_rows(data, col, value, top):
    if col and value:
        if col not in data.columns:
            raise PlotError(f"--col '{col}' is not a column of the input ({', '.join(map(str, data.columns))}).")
        return data.loc[data[col].isin(value.split(","))]
    return data.iloc[0:int(top), :]


def _output_path(cfg: PlotsConfig, name: str, default: str) -> Path:
    name = name or default
    if not name.lower().endswith(".pdf"):
        name += ".pdf"
    path = Path(name)
    if not path.is_absolute() and path.parent == Path("."):
        path = Path(cfg.output_dir) / path
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def heatmap_plot(cfg: PlotsConfig, input_files: str, samples: str, show_abund: bool, cluster_map: bool,
                 palette: str, col: str, value: str, plot_style: str, top: int, show_fig: bool):
    import seaborn as sns
    tables, samples_list = _read_abundance(input_files, samples, plot_style)
    if len(tables) == 1:
        data = tables[0].rename(columns={tables[0].columns[1]: samples_list[0]})
    else:
        data, _ = combined_datafram_based_on_first_col(tables, samples_list)
    data = _select_rows(data, col, value, top)
    data = data.set_index(data.iloc[:, 0])[samples_list]
    sns.set_style("whitegrid")
    sns.set_context("paper")

    n_rows, n_cols = data.shape
    fig_width = min(max(n_cols * 0.7, 6), 30)
    fig_height = min(max(n_rows * 0.5, 4), 30)
    if palette:
        cmap = palette
    else:
        cmap = colors.LinearSegmentedColormap.from_list(
            "my_list", ["aliceblue", "skyblue", "deepskyblue", "orange", "tomato", "red"])
        cmap.set_under("white")

    if cluster_map:
        g = sns.clustermap(data, cmap=cmap, cbar=True, vmin=0.1, dendrogram_ratio=0.03, cbar_pos=(0.1, 1, 0.1, 0.1),
                           col_cluster=False, cbar_kws={"shrink": 0.3}, figsize=(fig_width, fig_height))
        plt.setp(g.ax_heatmap.get_xticklabels(), rotation=30, ha="right", fontsize=10)
        plt.setp(g.ax_heatmap.get_yticklabels(), fontsize=10)
        default = P.HEATMAP_CLUSTER_FILE
    else:
        plt.figure(figsize=(fig_width, fig_height))
        ax = sns.heatmap(data, cmap=cmap, yticklabels=True, annot=show_abund, fmt=".0f",
                         linewidths=0 if not show_abund else 0.5, cbar=True, vmin=0.1,
                         cbar_kws={"shrink": 0.3, "anchor": (0, 0.0)})
        ax.collections[0].colorbar.ax.tick_params(labelsize=8)
        plt.xticks(rotation=30, ha="right", fontsize=10)
        plt.yticks(fontsize=10)
        plt.tight_layout(pad=0.5)
        default = P.HEATMAP_FILE
    if show_fig:
        plt.show()
        return []
    out = _output_path(cfg, cfg.pdf, default)
    plt.savefig(out, bbox_inches="tight")
    plt.close("all")
    logger.info(f"Saved plot: {out}")
    return [out]


def bar_plot(cfg: PlotsConfig, input_files: str, samples: str, plot_style: str, col: str, value: str, top: int,
             vertical_bar: bool, pdf: str):
    tables, samples_list = _read_abundance(input_files, samples, plot_style)
    result = combined_datafram_based_on_first_col(tables, samples_list)
    if isinstance(result, tuple):
        data, x = result
    else:
        data = result.rename(columns={result.columns[1]: samples_list[0]})
        x = data.columns[0]
    data = _select_rows(data, col, value, top)

    if vertical_bar:
        data.plot(x=x, y=samples_list, kind="bar")
        plt.xticks(rotation=90)
        plt.xlabel("")
        plt.ylabel("Abundance")
    else:
        data.plot.barh(x=x, y=samples_list)
        plt.gca().invert_yaxis()          # largest difference on top
        plt.ylabel("")
        plt.xlabel("Abundance")

    lower = input_files.lower()
    title_word = "subfamilies" if "subfam" in lower else ("ECs" if "ec" in os.path.basename(lower) else "families")
    plt.title(f"The most top {top} different {title_word}")
    out = _output_path(cfg, pdf, P.BAR_PLOT_FILE)
    plt.savefig(out, bbox_inches="tight")
    plt.close("all")
    logger.info(f"Saved plot: {out}")
    return [out]


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
click.rich_click.USE_RICH_MARKUP = True
click.rich_click.USE_MARKDOWN = True
click.rich_click.SHOW_ARGUMENTS = True
click.rich_click.GROUP_ARGUMENTS_OPTIONS = True
click.rich_click.STYLE_ERRORS_SUGGESTION = "magenta italic"
click.rich_click.SHOW_METAVARS_COLUMN = False
click.rich_click.APPEND_METAVARS_HELP = True
click.rich_click.MAX_WIDTH = 100


@click.command()
@logging_options
@click.argument("function", type=click.Choice([
    "CGC_plot", "CGC_coverage_plot", "CGC_synteny_plot", "CGC_synteny_coverage_plot", "bar_plot", "heatmap_plot"]))
@click.option("-i", "--input", "input_path", required=True,
              help="run_dbcan output directory (CGC plots), or comma-separated abundance files (bar/heatmap).")
@click.option("-o", "--output-dir", "--output_dir", "output_dir", default=".", show_default=True,
              help="Directory for the figures (created if needed).")
@click.option("--db-dir", "--db_dir", default="db", show_default=True,
              help="dbCAN database directory (dbCAN-PUL for the synteny plots).")
@click.option("--cgcid", "cgcid", default=None, help="CGC identifier contig|CGCn (example: 'contig1|CGC1').")
@click.option("--pul", "pul", default=None,
              help="dbCAN-PUL id to compare with (synteny plots; default: the CGC's best hit in substrate_prediction.tsv).")
@click.option("--reads-count", "--reads_count", "--readscount", "reads_count", multiple=True,
              help="samtools depth file (contig, position, depth); repeat for several samples (one panel each).")
@click.option("--samples", "samples", default=None,
              help="Comma-separated sample names: one per bar/heatmap input file or per --reads-count file.")
@click.option("--format", "--formats", "formats", default="pdf", show_default=True,
              help="Figure format(s) of the CGC plots, comma-separated: pdf, png, svg.")
@click.option("--dpi", default=P.PLOT_DEFAULT_DPI, show_default=True, type=int, help="Resolution of PNG output.")
@click.option("--show-protein-ids", "--show_protein_ids", is_flag=True,
              help="CGC_plot / CGC_coverage_plot: also print the protein IDs under the genes.")
@click.option("--synteny-revcomp/--no-synteny-revcomp", "--synteny_revcomp/--no-synteny_revcomp",
              "synteny_revcomp", default=True, show_default=True,
              help="Synteny plots: draw the dbCAN-PUL reverse-complemented when most homologous genes lie on the opposite strand.")
@click.option("--evalue-cutoff", "--evalue_cutoff", default=CSP.DEFAULT_EVALUE_CUTOFF, show_default=True, type=float,
              help="Synteny links: maximum DIAMOND E-value (as in substrate_prediction).")
@click.option("--identity-cutoff", "--identity_cutoff", default=CSP.DEFAULT_IDENTITY_CUTOFF, show_default=True,
              type=float, help="Synteny links: minimum DIAMOND percent identity, 0-100 (as in substrate_prediction).")
@click.option("--coverage-cutoff", "--coverage_cutoff", default=CSP.DEFAULT_COVERAGE_CUTOFF, show_default=True,
              type=float, help="Synteny links: minimum query coverage 0-1 (as in substrate_prediction).")
@click.option("--top", "top", default=20, show_default=True, type=int, help="Maximum number of categories to display.")
@click.option("--plot-style", "--plot_style", default="ggplot", show_default=True,
              type=click.Choice(list(P.PLOT_STYLE_CHOICES) + list(P.PLOT_STYLE_ALIASES)),
              help="Matplotlib style for bar/heatmap ('seaborn' is mapped to 'seaborn-v0_8').")
@click.option("--vertical-bar", "--vertical_bar", is_flag=True,
              help="bar_plot: vertical bars (default: horizontal bars).")
@click.option("--show-fig", "--show_fig", is_flag=True, help="Display the heatmap window instead of saving it.")
@click.option("--show-abund", "--show_abund", is_flag=True, help="Annotate heatmap cells with the abundance values.")
@click.option("--palette", "palette", default=None, help="Heatmap colormap name (e.g. viridis, Reds).")
@click.option("--cluster-map", "--cluster_map", is_flag=True, help="Cluster heatmap rows.")
@click.option("--col", "col", default=None, help="Column name used with --value for abundance filtering.")
@click.option("--value", "value", default=None, help="Comma-separated values to keep from column --col.")
@click.option("--pdf", "pdf", default=None,
              help="PDF file name of bar_plot (default bar_plot.pdf) or heatmap_plot (default heatmap.pdf); "
                   "a bare name is written into --output-dir.")
@click.option("--show-annotation", "--show_annotation", is_flag=True, hidden=True,
              help="Deprecated: family labels are always drawn.")
def cli(function, log_level, log_file, verbose, input_path, output_dir, db_dir, cgcid, pul, reads_count, samples,
        formats, dpi, show_protein_ids, synteny_revcomp, evalue_cutoff, identity_cutoff, coverage_cutoff, top,
        plot_style, vertical_bar, show_fig, show_abund, palette, cluster_map, col, value, pdf, show_annotation):
    """# dbCAN plotting utilities

    ## CGC figures (input: a run_dbcan output directory)

    **CGC_plot**: gene map of one CGC
    ```
    dbcan_plot CGC_plot -i sample.dbCAN --cgcid 'contigX|CGC1' -o figures --format pdf,png
    ```

    **CGC_coverage_plot**: gene map with read depth (one panel per sample)
    ```
    dbcan_plot CGC_coverage_plot -i sample.dbCAN --cgcid 'contigX|CGC1' --reads-count a.depth.txt --reads-count b.depth.txt --samples A,B
    ```

    **CGC_synteny_plot**: CGC vs its best dbCAN-PUL hit
    ```
    dbcan_plot CGC_synteny_plot -i sample.dbCAN --cgcid 'contigX|CGC1' --db-dir db
    ```

    **CGC_synteny_coverage_plot**: synteny with read depth
    ```
    dbcan_plot CGC_synteny_coverage_plot -i sample.dbCAN --cgcid 'contigX|CGC1' --db-dir db --reads-count a.depth.txt
    ```

    ## Abundance charts

    ```
    dbcan_plot bar_plot -i a/fam_abund.out,b/fam_abund.out --samples A,B --vertical-bar -o figures
    dbcan_plot heatmap_plot -i a/fam_abund.out,b/fam_abund.out --samples A,B --show-abund --pdf fam_heatmap.pdf
    ```
    """
    setup_logging(log_level, log_file, verbose)
    if show_annotation:
        logger.warning("--show-annotation is deprecated: family labels are always drawn on CGC plots.")
    try:
        cfg = PlotsConfig(
            input_dir=input_path if os.path.isdir(input_path) else (os.path.dirname(os.path.abspath(input_path)) or "."),
            db_dir=db_dir, output_dir=output_dir, formats=tuple(V.normalize_formats([formats])), dpi=dpi,
            cgcid=cgcid, pul=pul, reads_count=list(reads_count), show_protein_ids=show_protein_ids,
            synteny_revcomp=synteny_revcomp, evalue_cutoff=evalue_cutoff, identity_cutoff=identity_cutoff,
            coverage_cutoff=coverage_cutoff, samples=(samples.split(",") if samples else None), top=top,
            plot_style=plot_style, vertical_bar=vertical_bar, show_fig=show_fig, show_abund=show_abund,
            palette=palette, cluster_map=cluster_map, filter_col=col, filter_value=value, pdf=pdf,
        )
        if function in ("CGC_plot", "CGC_coverage_plot", "CGC_synteny_plot", "CGC_synteny_coverage_plot") \
                and not os.path.isdir(input_path):
            raise PlotError(f"-i must be a run_dbcan output directory for {function}: {input_path}")
        if function == "CGC_plot":
            CGC_plot(cfg)
        elif function == "CGC_coverage_plot":
            CGC_plot_reads_count(cfg)
        elif function == "CGC_synteny_plot":
            CGC_syntenic_with_PUL(cfg)
        elif function == "CGC_synteny_coverage_plot":
            CGC_syntenic_with_PUL_abund(cfg)
        else:
            check_input_files_min_lines(input_path, min_lines=2)
            if not samples:
                raise PlotError(f"--samples is required for {function}")
            if function == "bar_plot":
                bar_plot(cfg, input_files=input_path, samples=samples, plot_style=plot_style, col=col, value=value,
                         top=top, vertical_bar=vertical_bar, pdf=pdf)
            else:
                heatmap_plot(cfg, input_files=input_path, samples=samples, show_abund=show_abund,
                             cluster_map=cluster_map, palette=palette, col=col, value=value, plot_style=plot_style,
                             top=top, show_fig=show_fig)
    except click.ClickException:
        raise
    except V.PlotInputError as e:
        raise PlotError(str(e))
    except (OSError, ValueError, KeyError, pd.errors.ParserError) as e:
        logger.debug("dbcan_plot failed", exc_info=True)
        raise PlotError(f"{function} failed: {e}")


def main():
    cli()


if __name__ == "__main__":
    cli()
