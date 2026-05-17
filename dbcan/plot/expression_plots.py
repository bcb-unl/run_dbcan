"""CGC expression and DEG visualization."""

from __future__ import annotations

import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Set

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Polygon

from dbcan.constants.plots_constants import CGC_RESULT_FILE
from dbcan.utils.plots import (
    CGC,
    CGC_Standard_Out,
    CGC_standard_out_2CGC,
    Get_Position,
    plot_Polygon,
    plot_genome_line,
    plot_scale_line,
    points2,
)

logger = logging.getLogger(__name__)

DEG_UP_COLOR = "#D62728"
DEG_DOWN_COLOR = "#1F77B4"
GENE_COLORS = {
    "CAZyme": "#E67E22",
    "TC": "#2ECC71",
    "TF": "#9B59B6",
    "STP": "#F1C40F",
    "Peptidase": "#16A085",
    "Sulfatase": "#34495E",
}
DEFAULT_GENE_COLOR = "#95A5A6"


def _load_gene_deg(expression_dir: Path) -> pd.DataFrame:
    path = expression_dir / "gene_deg.tsv"
    if not path.is_file():
        return pd.DataFrame(columns=["protein_id", "log2FoldChange", "padj", "is_deg"])
    return pd.read_csv(path, sep="\t")


def _load_cgc_deg(expression_dir: Path) -> pd.DataFrame:
    path = expression_dir / "cgc_deg.tsv"
    if not path.is_file():
        return pd.DataFrame(columns=["cgcid", "is_cgc_de"])
    return pd.read_csv(path, sep="\t")


def _load_cgc_gene_counts(expression_dir: Path) -> pd.DataFrame:
    path = expression_dir / "cgc_gene_counts.tsv"
    if not path.is_file():
        raise FileNotFoundError(f"cgc_gene_counts.tsv not found in {expression_dir}")
    return pd.read_csv(path, sep="\t")


def _gene_x_positions(starts: List[int], ends: List[int]) -> List[float]:
    """Pixel x centers matching Get_Position layout."""
    shift = min(starts)
    maxbp = max(ends) - shift
    width = 1000
    pix = width / maxbp if maxbp > 0 else 1.0
    return [((s - shift) + (e - shift)) / 2 * pix for s, e in zip(starts, ends)]


def _plot_expression_panel(
    ax,
    cgc: CGC,
    counts_df: pd.DataFrame,
    metric: str,
    sample_ids: List[str],
    show_error_bars: bool,
    replicate_groups: Optional[Dict[str, str]],
):
    """Top panel: multi-sample gene expression lines."""
    starts, ends, _ = cgc.get_positions()
    protein_ids = cgc.get_proteinID()
    xs = _gene_x_positions(starts, ends)

    val_col = "tpm" if metric == "tpm" else "count"
    palette = sns.color_palette("colorblind", n_colors=max(len(sample_ids), 3))

    for si, sample_id in enumerate(sample_ids):
        sub = counts_df[
            (counts_df["cgcid"] == cgc.ID)
            & (counts_df["sample_id"] == sample_id)
        ]
        pid2val = dict(zip(sub["protein_id"], sub[val_col]))
        ys = [float(pid2val.get(pid, 0)) for pid in protein_ids]
        color = palette[si % len(palette)]
        ax.plot(xs, ys, "-o", color=color, label=sample_id, linewidth=1.5, markersize=5)

    if show_error_bars and replicate_groups:
        groups: Dict[str, List[str]] = {}
        for sid, grp in replicate_groups.items():
            groups.setdefault(grp, []).append(sid)
        for gi, (grp, members) in enumerate(groups.items()):
            if len(members) < 2:
                continue
            mean_y = []
            sem_y = []
            for pid in protein_ids:
                vals = []
                for sid in members:
                    sub = counts_df[
                        (counts_df["cgcid"] == cgc.ID)
                        & (counts_df["sample_id"] == sid)
                        & (counts_df["protein_id"] == pid)
                    ]
                    if len(sub):
                        vals.append(float(sub[val_col].iloc[0]))
                mean_y.append(np.mean(vals) if vals else 0)
                sem_y.append(np.std(vals, ddof=1) / np.sqrt(len(vals)) if len(vals) > 1 else 0)
            ax.errorbar(
                xs, mean_y, yerr=sem_y, fmt="none", ecolor="gray", alpha=0.5, capsize=2
            )

    ax.set_ylabel("Counts" if metric == "counts" else "TPM")
    ax.legend(loc="upper right", fontsize=7, frameon=False)
    for pos in ["top", "right"]:
        ax.spines[pos].set_visible(False)
    ax.set_xticks([])


def _plot_gene_panel_with_deg(
    ax,
    cgc: CGC,
    gene_deg: pd.DataFrame,
    deg_marker: str,
):
    """Bottom panel: gene arrows with DEG edge highlighting."""
    starts, ends, strands = cgc.get_positions()
    types = cgc.get_cgc_CAZyme()
    labels = cgc.get_proteinID()

    polygens, lines, texts, scale_positions, scale_text = Get_Position(
        starts, ends, strands, labels
    )

    deg_idx = gene_deg.set_index("protein_id") if len(gene_deg) else pd.DataFrame()

    for j, poly_str in enumerate(polygens):
        polygen = poly_str.split()
        points = []
        gtype = types[j] if j < len(types) else "Other"
        face = GENE_COLORS.get(gtype, DEFAULT_GENE_COLOR)
        edge = "none"
        lw = 0
        pid = labels[j] if j < len(labels) else ""
        if pid in deg_idx.index and bool(deg_idx.loc[pid].get("is_deg", False)):
            lfc = float(deg_idx.loc[pid].get("log2FoldChange", 0) or 0)
            edge = DEG_UP_COLOR if lfc >= 0 else DEG_DOWN_COLOR
            lw = 2.5
        for i in range(len(polygen) // 2):
            points.append([float(polygen[2 * i]), float(polygen[2 * i + 1])])
        ax.add_patch(
            Polygon(points, facecolor=face, edgecolor=edge, alpha=0.55, linewidth=lw)
        )
        if deg_marker in ("arrow", "star") and lw > 0 and j < len(labels):
            cx = (points[0][0] + points[-1][0]) / 2
            cy = max(p[1] for p in points) + 8
            sym = "▲" if edge == DEG_UP_COLOR else "▼"
            if deg_marker == "star":
                sym = "*"
            ax.text(cx, cy, sym, ha="center", va="bottom", fontsize=9, color=edge)

    plot_genome_line(lines, ax)
    plot_scale_line(scale_positions, scale_text, ax)

    genelabelcolor = list(GENE_COLORS.values()) + [DEFAULT_GENE_COLOR]
    geneslabels = list(GENE_COLORS.keys()) + ["Other"]
    handles = [Patch(color=c, alpha=0.5) for c in genelabelcolor]
    handles.append(Line2D([0], [0], color=DEG_UP_COLOR, lw=2, label="DEG up"))
    handles.append(Line2D([0], [0], color=DEG_DOWN_COLOR, lw=2, label="DEG down"))
    labels_leg = geneslabels + ["DEG up", "DEG down"]
    ax.legend(handles, labels_leg, frameon=False, loc="upper right", fontsize=7)
    ax.set_ylim(0, 150)
    ax.set_xlim(-10, 1100)
    ax.axis("off")
    ax.set_xlabel("Gene Position", fontsize=9)


def plot_cgc_expression_panel(
    cgc: CGC,
    counts_df: pd.DataFrame,
    gene_deg: pd.DataFrame,
    cgc_deg_row: Optional[pd.Series],
    out_pdf: Path,
    metric: str = "counts",
    deg_marker: str = "edge",
    show_error_bars: bool = False,
    replicate_groups: Optional[Dict[str, str]] = None,
) -> None:
    """Render dual-panel CGC expression figure."""
    sample_ids = sorted(counts_df["sample_id"].unique())
    px = 1 / plt.rcParams["figure.dpi"]
    fig = plt.figure(figsize=(14 * 1.2, 8))
    ax_top = fig.add_subplot(211)
    ax_bot = fig.add_subplot(212)

    _plot_expression_panel(
        ax_top, cgc, counts_df, metric, sample_ids, show_error_bars, replicate_groups
    )

    badge = ""
    if cgc_deg_row is not None:
        if bool(cgc_deg_row.get("is_cgc_de", False)):
            padj = cgc_deg_row.get("cgc_padj", "")
            lfc = cgc_deg_row.get("cgc_log2FoldChange", "")
            badge = f"  [DE CGC  padj={padj:.3g}  log2FC={lfc:.2f}]" if pd.notna(padj) else "  [DE CGC]"
        else:
            badge = "  [Not DE]"
    fig.text(0.5, 0.48, f"CGCID: {cgc.ID}{badge}", ha="center", fontsize=11, fontweight="bold")

    _plot_gene_panel_with_deg(ax_bot, cgc, gene_deg, deg_marker)
    fig.subplots_adjust(hspace=0.35)
    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_pdf, bbox_inches="tight", dpi=150)
    plt.close(fig)
    logger.info(f"Saved {out_pdf}")


def _list_cgc_ids(dbcan_dir: str, expression_dir: Path, cgcid: Optional[str], only_de: bool, top: Optional[int]) -> List[str]:
    pul_ann = Path(dbcan_dir) / CGC_RESULT_FILE
    dbcan = CGC_Standard_Out(str(pul_ann))
    all_ids = list(CGC_standard_out_2CGC(dbcan).cgcid2CGC().keys())
    if cgcid:
        return [cgcid] if cgcid in all_ids else []
    ids = all_ids
    if only_de:
        cgc_deg = _load_cgc_deg(expression_dir)
        if len(cgc_deg):
            de_set = set(cgc_deg[cgc_deg["is_cgc_de"].astype(bool)]["cgcid"])
            ids = [i for i in ids if i in de_set]
    if top is not None and top > 0:
        cgc_deg = _load_cgc_deg(expression_dir)
        if len(cgc_deg) and "cgc_log2FoldChange" in cgc_deg.columns:
            ranked = cgc_deg.set_index("cgcid").reindex(ids)
            ranked["_abs"] = ranked["cgc_log2FoldChange"].abs().fillna(0)
            ids = list(ranked.sort_values("_abs", ascending=False).head(top).index)
        else:
            ids = ids[:top]
    return ids


def _plot_one_cgc(args_tuple):
    cgcid, dbcan_dir, expression_dir, out_dir, metric, deg_marker, show_error_bars, repl = args_tuple
    pul_ann = Path(dbcan_dir) / CGC_RESULT_FILE
    cgcs = CGC_standard_out_2CGC(CGC_Standard_Out(str(pul_ann)))
    cgc = cgcs.cgcid2CGC().get(cgcid)
    if not cgc:
        return None
    counts_df = _load_cgc_gene_counts(Path(expression_dir))
    gene_deg = _load_gene_deg(Path(expression_dir))
    cgc_deg = _load_cgc_deg(Path(expression_dir))
    row = None
    if len(cgc_deg):
        m = cgc_deg[cgc_deg["cgcid"] == cgcid]
        row = m.iloc[0] if len(m) else None
    out_pdf = Path(out_dir) / f"{cgcid.replace('|', '_')}.expression.pdf"
    plot_cgc_expression_panel(
        cgc, counts_df, gene_deg, row, out_pdf, metric, deg_marker, show_error_bars, repl
    )
    return str(out_pdf)


def run_cgc_expression_plots(
    dbcan_dir: str,
    expression_dir: str,
    cgcid: Optional[str] = None,
    only_de: bool = False,
    top: Optional[int] = None,
    metric: str = "counts",
    deg_marker: str = "edge",
    output_dir: Optional[str] = None,
    show_error_bars: bool = False,
    replicate_groups: Optional[Dict[str, str]] = None,
    max_cgc: int = 500,
    force: bool = False,
    threads: int = 1,
) -> List[str]:
    """Plot one or all CGCs; returns list of output PDF paths."""
    expr = Path(expression_dir)
    out = Path(output_dir or expr / "plots" / "CGC_expression")
    out.mkdir(parents=True, exist_ok=True)

    ids = _list_cgc_ids(dbcan_dir, expr, cgcid, only_de, top)
    if not ids:
        logger.error("No CGC ids to plot")
        return []

    if len(ids) > max_cgc and not force:
        raise ValueError(
            f"Would plot {len(ids)} CGCs (max {max_cgc}). Use --top, --only-de, --cgcid, or --force"
        )

    outputs = []
    if threads <= 1 or len(ids) == 1:
        for cid in ids:
            r = _plot_one_cgc(
                (cid, dbcan_dir, str(expr), str(out), metric, deg_marker, show_error_bars, replicate_groups)
            )
            if r:
                outputs.append(r)
    else:
        tasks = [
            (cid, dbcan_dir, str(expr), str(out), metric, deg_marker, show_error_bars, replicate_groups)
            for cid in ids
        ]
        with ThreadPoolExecutor(max_workers=threads) as pool:
            futs = {pool.submit(_plot_one_cgc, t): t[0] for t in tasks}
            for fut in as_completed(futs):
                r = fut.result()
                if r:
                    outputs.append(r)
    logger.info(f"Generated {len(outputs)} CGC expression plots in {out}")
    return outputs


def plot_deg_volcano(gene_deg_path: Path, out_pdf: Path) -> None:
    df = pd.read_csv(gene_deg_path, sep="\t")
    if df.empty:
        return
    fig, ax = plt.subplots(figsize=(8, 6))
    df = df.copy()
    df["neglog10p"] = -np.log10(df["pvalue"].clip(lower=1e-300))
    colors = ["#CCCCCC"] * len(df)
    if "is_deg" in df.columns:
        up = df["is_deg"] & (df["log2FoldChange"] > 0)
        down = df["is_deg"] & (df["log2FoldChange"] < 0)
        colors = np.where(up, DEG_UP_COLOR, np.where(down, DEG_DOWN_COLOR, "#CCCCCC"))
    ax.scatter(df["log2FoldChange"], df["neglog10p"], c=colors, s=12, alpha=0.7)
    ax.axhline(-np.log10(0.05), color="gray", ls="--", lw=0.8)
    ax.set_xlabel("log2 Fold Change")
    ax.set_ylabel("-log10(pvalue)")
    ax.set_title("Gene-level DEG volcano")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Saved {out_pdf}")


def plot_deg_cgc_bar(cgc_deg_path: Path, out_pdf: Path) -> None:
    df = pd.read_csv(cgc_deg_path, sep="\t")
    if df.empty:
        return
    n_de = int(df["is_cgc_de"].sum()) if "is_cgc_de" in df.columns else 0
    n_not = len(df) - n_de
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.bar(["DE CGC", "Non-DE CGC"], [n_de, n_not], color=[DEG_UP_COLOR, "#CCCCCC"])
    ax.set_ylabel("Count")
    ax.set_title("CGC differential expression summary")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Saved {out_pdf}")


def plot_cgc_expression_heatmap(
    expression_dir: Path,
    out_pdf: Path,
    top: int = 50,
    cluster_map: bool = False,
) -> None:
    """Optional heatmap of CGC mean expression across samples."""
    matrix_path = expression_dir / "cgc_expression_matrix.tsv"
    if not matrix_path.is_file():
        logger.warning(f"Missing {matrix_path}")
        return
    df = pd.read_csv(matrix_path, sep="\t", index_col=0)
    if top and len(df) > top:
        df = df.loc[df.mean(axis=1).sort_values(ascending=False).head(top).index]
    fig, ax = plt.subplots(figsize=(max(6, df.shape[1] * 0.5), max(8, df.shape[0] * 0.15)))
    if cluster_map:
        sns.clustermap(df, cmap="viridis", figsize=fig.get_size_inches())
        plt.savefig(out_pdf, bbox_inches="tight")
    else:
        sns.heatmap(df, ax=ax, cmap="viridis")
        fig.savefig(out_pdf, bbox_inches="tight")
    plt.close("all")
    logger.info(f"Saved {out_pdf}")
