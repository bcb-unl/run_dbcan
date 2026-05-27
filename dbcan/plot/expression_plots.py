"""CGC expression and DEG visualization (SeMa-Trap style)."""

from __future__ import annotations

import logging
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import to_rgb
from matplotlib.patches import Patch, Polygon, Rectangle

from dbcan.constants.expression_constants import (
    COUNT_MATRIX_FILE,
    GENE_NORM_COUNTS_FILE,
    PLOT_METRIC_CHOICES,
    SAMPLE_METADATA_FILE,
)
from dbcan.constants.plots_constants import CGC_RESULT_FILE
from dbcan.constants.utils_constants import OVERVIEW_FILE
from dbcan.utils.plots import (
    CGC,
    CGC_Standard_Out,
    CGC_standard_out_2CGC,
    Get_Position,
    add_gene_annotations,
    plot_genome_line,
    plot_scale_line,
)

logger = logging.getLogger(__name__)

GENOME_PLOT_WIDTH = 1000
X_LIM = (-10, 1100)
GENE_LOCUS_LINE_COLOR = "#E8E8E8"
GENE_TRACK_RECT_HEIGHT_FRAC = 0.03
GENE_TRACK_RECT_HEIGHT_MIN = 0.06
GENE_TRACK_RECT_HEIGHT_MAX = 0.15
GENE_TRACK_FILL_ALPHA = 0.18
GENE_TRACK_EDGE_WIDTH = 1.35
GENE_TRACK_LINK_WIDTH = 0.9
GENE_TRACK_LINK_ALPHA = 0.75
LEGEND_LINE_WIDTH = 1.6

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

_CAZYME_FAMILY_RE = re.compile(r"(GH|GT|PL|CE|AA|CBM|SLH|COH|dockerin)\d+", re.I)
_PFAM_RE = re.compile(r"PF\d{5}", re.I)


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


def _load_sample_metadata(expression_dir: Path) -> pd.DataFrame:
    path = expression_dir / SAMPLE_METADATA_FILE
    if not path.is_file():
        raise FileNotFoundError(f"{SAMPLE_METADATA_FILE} not found in {expression_dir}")
    df = pd.read_csv(path, sep="\t")
    if "sample_id" in df.columns:
        df = df.set_index("sample_id")
    return df


def _load_overview(dbcan_dir: str) -> pd.DataFrame:
    path = Path(dbcan_dir) / OVERVIEW_FILE
    if not path.is_file():
        return pd.DataFrame()
    df = pd.read_csv(path, sep="\t")
    if "Gene ID" in df.columns:
        df = df.set_index("Gene ID")
    return df


def _load_count_matrix(expression_dir: Path) -> Optional[pd.DataFrame]:
    path = expression_dir / COUNT_MATRIX_FILE
    if not path.is_file():
        return None
    df = pd.read_csv(path, sep="\t", index_col=0)
    return df


def _load_norm_counts(expression_dir: Path) -> Optional[pd.DataFrame]:
    path = expression_dir / GENE_NORM_COUNTS_FILE
    if not path.is_file():
        return None
    return pd.read_csv(path, sep="\t", index_col=0)


def gene_bp_to_pixel(
    starts: List[int],
    ends: List[int],
    width: int = GENOME_PLOT_WIDTH,
) -> Tuple[List[float], List[float], float, float]:
    """Map gene start/end bp to pixel x (same transform as Get_Position)."""
    shift = float(min(starts))
    maxbp = float(max(ends) - shift)
    pix = width / maxbp if maxbp > 0 else 1.0
    x_starts = [(s - shift) * pix for s in starts]
    x_ends = [(e - shift) * pix for e in ends]
    return x_starts, x_ends, shift, pix


def build_sem_trap_trace(
    x_starts: List[float],
    x_ends: List[float],
    y_values: List[float],
    baseline: float = 0.0,
) -> Tuple[np.ndarray, np.ndarray]:
    """Build step polyline: horizontal per gene, drop to baseline between genes."""
    xs: List[float] = []
    ys: List[float] = []
    for i, (x0, x1, y) in enumerate(zip(x_starts, x_ends, y_values)):
        if i > 0:
            xs.extend([xs[-1], x0])
            ys.extend([baseline, baseline])
        xs.extend([x0, x1])
        ys.extend([y, y])
    return np.asarray(xs, dtype=float), np.asarray(ys, dtype=float)


def _metric_ylabel(metric: str) -> str:
    labels = {
        "log2_tpm": "log2(TPM + 1)",
        "log2_fpkm": "log2(FPKM + 1)",
        "log2_rpkm": "log2(RPKM + 1)",
        "log2_rpm": "log2(RPM + 1)",
        "log2_cpm": "log2(CPM + 1)",
        "log2_norm": "log2(norm. counts + 1)",
        "log2fc": "log2 fold change",
        "gene_zscore": "gene-wise z-score",
        "tpm": "TPM",
        "fpkm": "FPKM",
        "rpkm": "RPKM",
        "rpm": "RPM",
        "counts": "raw counts",
    }
    return labels.get(metric, metric)


def _library_sizes(count_matrix: pd.DataFrame) -> Dict[str, float]:
    return {sid: float(count_matrix.loc[sid].sum()) for sid in count_matrix.index}


def _gene_lengths(cgc: CGC) -> Dict[str, int]:
    return {
        g.seqid: max(int(g.gene_end) - int(g.gene_start) + 1, 1)
        for g in cgc.genes
    }


def _abundance_from_row(
    row: pd.Series,
    pid: str,
    pid2len: Dict[str, int],
    lib_sizes: Dict[str, float],
    sample_id: str,
) -> Dict[str, float]:
    """Resolve TPM/FPKM/RPM/RPKM for one gene-sample (table columns or fallback)."""
    count = float(row.get("count", 0) or 0)
    length = pid2len.get(pid, 1)
    length_kb = max(length, 1) / 1000.0
    lib = max(lib_sizes.get(sample_id, 1.0), 1.0)
    norm_reads = lib / 1e6

    def _col(name: str, fallback: float) -> float:
        if name in row.index and pd.notna(row[name]):
            return float(row[name])
        return fallback

    fpkm_fb = count / norm_reads / length_kb if norm_reads > 0 else 0.0
    rpm_fb = count / norm_reads if norm_reads > 0 else 0.0
    tpm_fb = float(row.get("tpm", 0) or 0)
    if tpm_fb <= 0 and count > 0:
        tpm_fb = fpkm_fb

    return {
        "count": count,
        "tpm": _col("tpm", tpm_fb),
        "fpkm": _col("fpkm", fpkm_fb),
        "rpm": _col("rpm", rpm_fb),
        "rpkm": _col("rpkm", fpkm_fb),
    }


def _transform_expression(
    metric: str,
    row: pd.Series,
    sample_id: str,
    protein_id: str,
    pid2len: Dict[str, int],
    lib_sizes: Dict[str, float],
    norm_matrix: Optional[pd.DataFrame],
    gene_deg: pd.DataFrame,
) -> float:
    """Return y-value for a gene/sample under the chosen metric."""
    if metric == "log2fc":
        if len(gene_deg) and protein_id in gene_deg.set_index("protein_id").index:
            return float(gene_deg.set_index("protein_id").loc[protein_id, "log2FoldChange"])
        return float("nan")

    abund = _abundance_from_row(row, protein_id, pid2len, lib_sizes, sample_id)
    count = abund["count"]

    if metric == "log2_tpm":
        return float(np.log2(abund["tpm"] + 1.0))
    if metric == "log2_fpkm":
        return float(np.log2(abund["fpkm"] + 1.0))
    if metric == "log2_rpkm":
        return float(np.log2(abund["rpkm"] + 1.0))
    if metric == "log2_rpm":
        return float(np.log2(abund["rpm"] + 1.0))
    if metric == "tpm":
        return abund["tpm"]
    if metric == "fpkm":
        return abund["fpkm"]
    if metric == "rpkm":
        return abund["rpkm"]
    if metric == "rpm":
        return abund["rpm"]
    if metric == "counts":
        return count
    if metric == "log2_cpm":
        lib = lib_sizes.get(sample_id, 1.0) or 1.0
        cpm = (count / lib) * 1e6
        return float(np.log2(cpm + 1.0))
    if metric == "log2_norm":
        if norm_matrix is not None and sample_id in norm_matrix.index and protein_id in norm_matrix.columns:
            return float(np.log2(float(norm_matrix.loc[sample_id, protein_id]) + 1.0))
        lib = lib_sizes.get(sample_id, 1.0) or 1.0
        cpm = (count / lib) * 1e6
        return float(np.log2(cpm + 1.0))
    return float(np.log2(abund["tpm"] + 1.0))


def _condition_order(conditions: List[str]) -> List[str]:
    def sort_key(c: str):
        cl = c.lower()
        if cl in ("control", "ctrl", "wt", "wildtype", "wild_type"):
            return (0, c)
        return (1, c)

    return sorted(set(conditions), key=sort_key)


def format_gene_label(
    gene_type: str,
    annotation: Optional[str],
    protein_id: str,
    overview_row: Optional[pd.Series] = None,
) -> str:
    """Short functional label for gene arrows and heatmap columns."""
    ann = (annotation or "").strip()
    gtype = (gene_type or "Other").strip()
    if ann and ann.lower() != "null" and "|" in ann:
        prefix, detail = ann.split("|", 1)
        prefix = prefix.strip()
        detail = detail.strip()
        if prefix == "CAZyme":
            m = _CAZYME_FAMILY_RE.search(detail)
            return m.group(0) if m else detail.split("_")[0][:12]
        if prefix == "TC":
            parts = detail.split(".")
            short = ".".join(parts[:3]) if len(parts) >= 3 else detail
            return f"TC:{short}"[:16]
        if detail:
            return detail[:16]
    if gtype == "CAZyme" and overview_row is not None:
        for col in ("Recommend Results", "dbCAN_sub", "DIAMOND"):
            val = str(overview_row.get(col, "") or "")
            m = _CAZYME_FAMILY_RE.search(val)
            if m:
                return m.group(0)
    if overview_row is not None:
        for col in ("Recommend Results", "dbCAN_sub", "DIAMOND", "dbCAN_hmm"):
            val = str(overview_row.get(col, "") or "")
            if val in ("-", "", "nan"):
                continue
            m = _PFAM_RE.search(val)
            if m:
                return m.group(0)
            m = _CAZYME_FAMILY_RE.search(val)
            if m:
                return m.group(0)
            if val and val != "-":
                return val.split("(")[0].split("_")[0][:14]
    if gtype and gtype.lower() not in ("null", "other", ""):
        return gtype[:14]
    return "Other"


def _gene_labels_for_cgc(cgc: CGC, overview: pd.DataFrame) -> List[str]:
    types = cgc.get_cgc_CAZyme()
    annos = cgc.get_annotations()
    pids = cgc.get_proteinID()
    labels = []
    for i, pid in enumerate(pids):
        row = overview.loc[pid] if len(overview) and pid in overview.index else None
        labels.append(
            format_gene_label(
                types[i] if i < len(types) else "Other",
                annos[i] if i < len(annos) else None,
                pid,
                row,
            )
        )
    return labels


def _aggregate_by_condition(
    cgc: CGC,
    counts_df: pd.DataFrame,
    metadata: pd.DataFrame,
    metric: str,
    expression_dir: Path,
    gene_deg: pd.DataFrame,
) -> Dict[str, List[float]]:
    """Per-condition mean y per gene (genomic order)."""
    protein_ids = cgc.get_proteinID()
    pid2len = _gene_lengths(cgc)
    count_matrix = _load_count_matrix(expression_dir)
    norm_matrix = _load_norm_counts(expression_dir)
    lib_sizes = _library_sizes(count_matrix) if count_matrix is not None else {}

    if metric == "log2fc" and len(gene_deg):
        deg_idx = gene_deg.set_index("protein_id")
        return {
            "_single": [
                float(deg_idx.loc[pid, "log2FoldChange"]) if pid in deg_idx.index else 0.0
                for pid in protein_ids
            ]
        }

    sample_to_cond = metadata["condition"].to_dict()
    conditions = _condition_order(list(metadata["condition"].unique()))
    result: Dict[str, List[float]] = {c: [] for c in conditions}

    sub = counts_df[counts_df["cgcid"] == cgc.ID]
    for cond in conditions:
        sample_ids = [s for s, c in sample_to_cond.items() if c == cond]
        for pid in protein_ids:
            vals = []
            for sid in sample_ids:
                rows = sub[(sub["sample_id"] == sid) & (sub["protein_id"] == pid)]
                if len(rows):
                    r = rows.iloc[0]
                    v = _transform_expression(
                        metric,
                        r,
                        sid,
                        pid,
                        pid2len,
                        lib_sizes,
                        norm_matrix,
                        gene_deg,
                    )
                    if not np.isnan(v):
                        vals.append(v)
            result[cond].append(float(np.mean(vals)) if vals else 0.0)
    return result


def _apply_gene_zscore(cond_values: Dict[str, List[float]]) -> Dict[str, List[float]]:
    """Within each condition, z-score expression across genes in the CGC."""
    out: Dict[str, List[float]] = {}
    for cond, ys in cond_values.items():
        if cond == "_single":
            out[cond] = ys
            continue
        arr = np.asarray(ys, dtype=float)
        if len(arr) < 2:
            out[cond] = list(arr)
            continue
        mu = float(np.mean(arr))
        sd = float(np.std(arr, ddof=0))
        out[cond] = list((arr - mu) / sd) if sd > 1e-12 else [0.0] * len(arr)
    return out


def _darken_color(color, factor: float = 0.75):
    r, g, b = to_rgb(color)[:3]
    return (r * factor, g * factor, b * factor)


def _draw_gene_locus_guides(ax, x_starts: List[float], x_ends: List[float]) -> None:
    """Keep the expression panel clean; the lower gene map already shows spans."""
    return


def _track_rect_height(cond_values: Dict[str, List[float]], metric: str) -> float:
    """Fixed narrow height for gene-span rectangles (data coordinates)."""
    flat = [float(v) for ys in cond_values.values() for v in ys if not np.isnan(v)]
    if not flat:
        return 0.08
    if metric in ("log2fc", "gene_zscore"):
        span = 2.0 * max(abs(min(flat)), abs(max(flat)), 0.35)
    else:
        span = max(max(flat) * 1.12 + 0.05, 0.5)
    return min(
        max(span * GENE_TRACK_RECT_HEIGHT_FRAC, GENE_TRACK_RECT_HEIGHT_MIN),
        GENE_TRACK_RECT_HEIGHT_MAX,
    )


def _condition_visual_offsets(
    conditions: List[str],
    cond_values: Dict[str, List[float]],
    metric: str,
) -> Dict[str, float]:
    """Tiny display-only offsets so overlapping condition tracks remain legible."""
    if len(conditions) <= 1:
        return {conditions[0]: 0.0} if conditions else {}
    flat = [float(v) for ys in cond_values.values() for v in ys if not np.isnan(v)]
    if not flat:
        max_abs = 0.06
    elif metric in ("log2fc", "gene_zscore"):
        span = 2.0 * max(abs(min(flat)), abs(max(flat)), 0.35)
        max_abs = min(max(span * 0.02, 0.03), 0.09)
    else:
        span = max(max(flat) - min(flat), max(flat), 0.5)
        max_abs = min(max(span * 0.02, 0.03), 0.09)
    mid = (len(conditions) - 1) / 2.0
    return {
        cond: (mid - i) / mid * max_abs if mid else 0.0
        for i, cond in enumerate(conditions)
    }


def _draw_condition_trajectory(
    ax,
    x_starts: List[float],
    x_ends: List[float],
    y_values: List[float],
    color,
    rect_height: float,
    zorder: int = 2,
) -> None:
    """
    Gene-level expression track: fixed-height narrow rects at y=expression per gene span,
    linked by thin lines between adjacent edge centers (continuous trajectory).
    """
    ys = [0.0 if np.isnan(y) else float(y) for y in y_values]
    edge = _darken_color(color)
    r, g, b = to_rgb(color)[:3]
    half_h = rect_height / 2.0
    n = len(ys)

    for i, (x0, x1, y) in enumerate(zip(x_starts, x_ends, ys)):
        width = max(float(x1 - x0), 1.0)
        ax.add_patch(
            Rectangle(
                (x0, y - half_h),
                width,
                rect_height,
                facecolor=(r, g, b, GENE_TRACK_FILL_ALPHA),
                edgecolor=edge,
                linewidth=GENE_TRACK_EDGE_WIDTH,
                zorder=zorder + 1,
                clip_on=True,
            )
        )

    for i in range(n - 1):
        ax.plot(
            [x_ends[i], x_starts[i + 1]],
            [ys[i], ys[i + 1]],
            color=edge,
            linewidth=GENE_TRACK_LINK_WIDTH,
            alpha=GENE_TRACK_LINK_ALPHA,
            zorder=zorder,
            solid_capstyle="round",
        )


def _set_expression_ylim(ax, cond_values: Dict[str, List[float]], metric: str) -> None:
    flat = [v for ys in cond_values.values() for v in ys if not np.isnan(v)]
    if not flat:
        ax.set_ylim(0, 1)
        return
    if metric in ("log2fc", "gene_zscore"):
        ymax = max(abs(min(flat)), abs(max(flat)), 0.5)
        pad = ymax * 0.15
        ax.set_ylim(-ymax - pad, ymax + pad)
    else:
        ymax = max(flat)
        ax.set_ylim(0, ymax * 1.12 + 0.05)


def plot_sem_trap_expression_panel(
    ax,
    cgc: CGC,
    counts_df: pd.DataFrame,
    metadata: pd.DataFrame,
    metric: str,
    expression_dir: Path,
    gene_deg: pd.DataFrame,
) -> None:
    """
    Gene-level expression track panel (SeMa-Trap style).

    Per condition: narrow horizontal rects at each gene's genomic span (y = expression,
    fixed rect height); thin lines connect right-edge center to next gene's left-edge center.
    """
    starts, ends, _ = cgc.get_positions()
    x_starts, x_ends, _, _ = gene_bp_to_pixel(starts, ends)
    ax.set_facecolor("#FFFFFF")
    _draw_gene_locus_guides(ax, x_starts, x_ends)
    baseline = 0.0
    y_metric = metric

    if metric == "log2fc":
        if not len(gene_deg):
            logger.warning("gene_deg.tsv missing; falling back to log2_tpm for expression panel")
            y_metric = "log2_tpm"
        else:
            cond_values = _aggregate_by_condition(
                cgc, counts_df, metadata, metric, expression_dir, gene_deg
            )
            if "_single" in cond_values:
                ys = cond_values["_single"]
                color = sns.color_palette("colorblind")[0]
                rh = _track_rect_height({"_fc": ys}, "log2fc")
                _draw_condition_trajectory(ax, x_starts, x_ends, ys, color, rh, zorder=2)
                ax.axhline(0, color="#BBBBBB", lw=0.9, zorder=1)
                _set_expression_ylim(ax, {"_fc": ys}, "log2fc")
                ax.set_ylabel(_metric_ylabel("log2fc"))
                ax.plot([], [], color=_darken_color(color), lw=LEGEND_LINE_WIDTH, label="log2FC")
                ax.legend(loc="upper right", fontsize=7, frameon=False)
                ax.yaxis.grid(True, color="#EEEEEE", lw=0.6, zorder=0)
                ax.set_axisbelow(True)
                for pos in ["top", "right"]:
                    ax.spines[pos].set_visible(False)
                return

    cond_values = _aggregate_by_condition(
        cgc, counts_df, metadata, y_metric, expression_dir, gene_deg
    )
    if y_metric == "gene_zscore":
        cond_values = _apply_gene_zscore(cond_values)

    conditions = _condition_order(list(cond_values.keys()))
    palette = sns.color_palette("colorblind", n_colors=max(len(conditions), 3))
    rect_h = _track_rect_height(cond_values, y_metric)
    offsets = _condition_visual_offsets(conditions, cond_values, y_metric)
    legend_lines = []
    display_values: Dict[str, List[float]] = {}

    for ci, cond in enumerate(conditions):
        color = palette[ci % len(palette)]
        offset = offsets.get(cond, 0.0)
        display_values[cond] = [
            (0.0 if np.isnan(y) else float(y)) + offset
            for y in cond_values[cond]
        ]
        _draw_condition_trajectory(
            ax, x_starts, x_ends, display_values[cond], color, rect_h, zorder=2 + ci,
        )
        n = int((metadata["condition"] == cond).sum())
        label = f"{cond} (n={n})" if n > 1 else cond
        legend_lines.append(
            ax.plot(
                [], [], color=_darken_color(color), lw=LEGEND_LINE_WIDTH, label=label,
            )[0]
        )

    if y_metric in ("log2fc", "gene_zscore"):
        ax.axhline(baseline, color="#BBBBBB", lw=0.9, zorder=1)
    _set_expression_ylim(ax, display_values, y_metric)
    ax.set_ylabel(_metric_ylabel(y_metric))
    ax.legend(handles=legend_lines, loc="upper right", fontsize=7, frameon=False)
    ax.yaxis.grid(True, color="#EEEEEE", lw=0.6, zorder=0)
    ax.set_axisbelow(True)
    for pos in ["top", "right"]:
        ax.spines[pos].set_visible(False)


def _plot_gene_panel_with_deg(
    ax,
    cgc: CGC,
    gene_deg: pd.DataFrame,
    deg_marker: str,
    dbcan_dir: str,
    alpha: float = 0.05,
) -> None:
    """Bottom panel: gene arrows, annotations, DEG edge + markers."""
    starts, ends, strands = cgc.get_positions()
    types = cgc.get_cgc_CAZyme()
    labels = cgc.get_proteinID()
    overview = _load_overview(dbcan_dir)
    gene_labels = _gene_labels_for_cgc(cgc, overview)

    polygens, lines, texts, scale_positions, scale_text = Get_Position(
        starts, ends, strands, labels
    )

    deg_idx = gene_deg.set_index("protein_id") if len(gene_deg) else pd.DataFrame()
    show_marker = deg_marker in ("arrow", "star", "both")
    use_edge = deg_marker in ("edge", "both", "arrow", "star")

    for j, poly_str in enumerate(polygens):
        polygen = poly_str.split()
        points = []
        gtype = types[j] if j < len(types) else "Other"
        face = GENE_COLORS.get(gtype, DEFAULT_GENE_COLOR)
        edge = "none"
        lw = 0
        pid = labels[j] if j < len(labels) else ""
        is_deg = pid in deg_idx.index and bool(deg_idx.loc[pid].get("is_deg", False))
        if is_deg and use_edge:
            lfc = float(deg_idx.loc[pid].get("log2FoldChange", 0) or 0)
            edge = DEG_UP_COLOR if lfc >= 0 else DEG_DOWN_COLOR
            lw = 2.5
        for i in range(len(polygen) // 2):
            points.append([float(polygen[2 * i]), float(polygen[2 * i + 1])])
        ax.add_patch(
            Polygon(points, facecolor=face, edgecolor=edge, alpha=0.55, linewidth=lw)
        )
        if is_deg and show_marker and j < len(labels):
            cx = (points[0][0] + points[-1][0]) / 2
            cy = max(p[1] for p in points) + 8
            sym = "▲" if edge == DEG_UP_COLOR else "▼"
            if deg_marker == "star":
                sym = "*"
            padj = deg_idx.loc[pid].get("padj", np.nan)
            if pd.notna(padj) and float(padj) < alpha:
                sym += "*"
            ax.text(cx, cy, sym, ha="center", va="bottom", fontsize=9, color=edge)

    plot_genome_line(lines, ax)
    plot_scale_line(scale_positions, scale_text, ax)

    ylim_top = 150
    top_y = add_gene_annotations(
        ax, starts, ends, gene_labels, font_size=11, rotation=22, y_offset=8
    )
    if top_y:
        ylim_top = max(ylim_top, top_y + 5)

    genelabelcolor = list(GENE_COLORS.values()) + [DEFAULT_GENE_COLOR]
    geneslabels = list(GENE_COLORS.keys()) + ["Other"]
    handles = [Patch(color=c, alpha=0.5) for c in genelabelcolor]
    ax.legend(handles, geneslabels, frameon=False, loc="upper right", fontsize=7)
    ax.set_ylim(0, ylim_top)
    ax.set_xlim(*X_LIM)
    ax.axis("off")


def build_cgc_gene_matrix_with_labels(
    cgc: CGC,
    counts_df: pd.DataFrame,
    metadata: pd.DataFrame,
    metric: str,
    expression_dir: Path,
    gene_deg: pd.DataFrame,
    dbcan_dir: str,
    row_by: str = "sample",
) -> Tuple[pd.DataFrame, List[str]]:
    """Build heatmap matrix with short column labels."""
    protein_ids = cgc.get_proteinID()
    pid2len = _gene_lengths(cgc)
    overview = _load_overview(dbcan_dir)
    col_labels = _gene_labels_for_cgc(cgc, overview)

    count_matrix = _load_count_matrix(expression_dir)
    norm_matrix = _load_norm_counts(expression_dir)
    lib_sizes = _library_sizes(count_matrix) if count_matrix is not None else {}
    sub = counts_df[counts_df["cgcid"] == cgc.ID]
    rows: Dict[str, List[float]] = {}

    plot_metric = metric
    if row_by == "condition":
        cond_values = _aggregate_by_condition(
            cgc, counts_df, metadata, plot_metric, expression_dir, gene_deg
        )
        if plot_metric == "gene_zscore":
            cond_values = _apply_gene_zscore(cond_values)
        for cond in _condition_order(list(cond_values.keys())):
            if cond == "_single":
                continue
            rows[cond] = cond_values[cond]
    else:
        for sid in sorted(sub["sample_id"].unique()):
            ys = []
            for pid in protein_ids:
                rows_df = sub[(sub["sample_id"] == sid) & (sub["protein_id"] == pid)]
                if len(rows_df):
                    r = rows_df.iloc[0]
                    v = _transform_expression(
                        metric,
                        r,
                        sid,
                        pid,
                        pid2len,
                        lib_sizes,
                        norm_matrix,
                        gene_deg,
                    )
                    ys.append(0.0 if np.isnan(v) else v)
                else:
                    ys.append(0.0)
            rows[sid] = ys
        if plot_metric == "gene_zscore":
            for sid in list(rows.keys()):
                arr = np.asarray(rows[sid], dtype=float)
                mu, sd = float(np.mean(arr)), float(np.std(arr, ddof=0))
                rows[sid] = list((arr - mu) / sd) if sd > 1e-12 else [0.0] * len(arr)

    mat = pd.DataFrame(rows).T
    mat.columns = col_labels
    return mat, col_labels


def plot_cgc_gene_heatmap(
    cgc: CGC,
    matrix: pd.DataFrame,
    out_pdf: Path,
    metric: str,
    gene_deg: pd.DataFrame,
    protein_ids: List[str],
) -> None:
    """Per-CGC gene expression heatmap."""
    if matrix.empty:
        return
    n_genes = matrix.shape[1]
    fig_w = max(6, n_genes * 0.55)
    fig_h = max(3, matrix.shape[0] * 0.45 + 1.5)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    deg_idx = gene_deg.set_index("protein_id") if len(gene_deg) else pd.DataFrame()
    col_colors = []
    for pid in protein_ids:
        if pid in deg_idx.index and bool(deg_idx.loc[pid].get("is_deg", False)):
            lfc = float(deg_idx.loc[pid].get("log2FoldChange", 0) or 0)
            col_colors.append(DEG_UP_COLOR if lfc >= 0 else DEG_DOWN_COLOR)
        else:
            col_colors.append("white")

    sns.heatmap(
        matrix,
        ax=ax,
        cmap="YlOrRd",
        linewidths=0.5,
        linecolor="#EEEEEE",
        cbar_kws={"label": _metric_ylabel(metric)},
    )
    for tick, color in zip(ax.get_xticklabels(), col_colors):
        if color != "white":
            tick.set_color(color)
            tick.set_fontweight("bold")

    ax.set_xlabel("Genes (genomic order)")
    ax.set_ylabel("Sample" if matrix.index.name != "condition" else "Condition")
    ax.set_title(f"{cgc.ID} — gene expression", fontsize=10)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    fig.tight_layout()
    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_pdf, bbox_inches="tight", dpi=150)
    plt.close(fig)
    logger.info(f"Saved {out_pdf}")


def plot_cgc_expression_panel(
    cgc: CGC,
    counts_df: pd.DataFrame,
    gene_deg: pd.DataFrame,
    cgc_deg_row: Optional[pd.Series],
    out_pdf: Path,
    expression_dir: Path,
    dbcan_dir: str,
    metric: str = "log2_tpm",
    deg_marker: str = "both",
    heatmap_rows: str = "sample",
    alpha: float = 0.05,
) -> None:
    """Render dual-panel CGC expression figure + optional gene heatmap."""
    metadata = _load_sample_metadata(expression_dir)

    fig, (ax_top, ax_bot) = plt.subplots(
        2,
        1,
        sharex=True,
        figsize=(14 * 1.2, 8.4),
        gridspec_kw={"height_ratios": [7, 3], "hspace": 0.015},
    )
    fig.patch.set_facecolor("white")

    plot_sem_trap_expression_panel(
        ax_top, cgc, counts_df, metadata, metric, expression_dir, gene_deg
    )
    ax_top.set_xlim(*X_LIM)

    badge = ""
    if cgc_deg_row is not None:
        if bool(cgc_deg_row.get("is_cgc_de", False)):
            padj = cgc_deg_row.get("cgc_padj", "")
            lfc = cgc_deg_row.get("cgc_log2FoldChange", "")
            badge = f"  [DE CGC  padj={padj:.3g}  log2FC={lfc:.2f}]" if pd.notna(padj) else "  [DE CGC]"
        else:
            badge = "  [Not DE]"
    fig.text(
        0.012, 0.985, f"CGCID: {cgc.ID}{badge}",
        ha="left", va="top", fontsize=11, fontweight="bold", color="#333333",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.75, "pad": 2.0},
    )

    _plot_gene_panel_with_deg(ax_bot, cgc, gene_deg, deg_marker, dbcan_dir, alpha=alpha)

    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_pdf, bbox_inches="tight", dpi=150)
    plt.close(fig)
    logger.info(f"Saved {out_pdf}")

    base_name = out_pdf.name.replace(".expression.pdf", "")
    hm_path = out_pdf.with_name(f"{base_name}.gene_heatmap.pdf")
    mat, _ = build_cgc_gene_matrix_with_labels(
        cgc, counts_df, metadata, metric, expression_dir, gene_deg, dbcan_dir, row_by=heatmap_rows
    )
    plot_cgc_gene_heatmap(cgc, mat, hm_path, metric, gene_deg, cgc.get_proteinID())


def _list_cgc_ids(
    dbcan_dir: str, expression_dir: Path, cgcid: Optional[str], only_de: bool, top: Optional[int]
) -> List[str]:
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
    (
        cgcid,
        dbcan_dir,
        expression_dir,
        out_dir,
        metric,
        deg_marker,
        heatmap_rows,
        alpha,
    ) = args_tuple
    pul_ann = Path(dbcan_dir) / CGC_RESULT_FILE
    cgcs = CGC_standard_out_2CGC(CGC_Standard_Out(str(pul_ann)))
    cgc = cgcs.cgcid2CGC().get(cgcid)
    if not cgc:
        return None
    expr = Path(expression_dir)
    counts_df = _load_cgc_gene_counts(expr)
    gene_deg = _load_gene_deg(expr)
    cgc_deg = _load_cgc_deg(expr)
    row = None
    if len(cgc_deg):
        m = cgc_deg[cgc_deg["cgcid"] == cgcid]
        row = m.iloc[0] if len(m) else None
    out_pdf = Path(out_dir) / f"{cgcid.replace('|', '_')}.expression.pdf"
    plot_cgc_expression_panel(
        cgc,
        counts_df,
        gene_deg,
        row,
        out_pdf,
        expr,
        dbcan_dir,
        metric=metric,
        deg_marker=deg_marker,
        heatmap_rows=heatmap_rows,
        alpha=alpha,
    )
    hm = out_pdf.with_name(f"{cgcid.replace('|', '_')}.gene_heatmap.pdf")
    return [str(out_pdf), str(hm)]


def run_cgc_expression_plots(
    dbcan_dir: str,
    expression_dir: str,
    cgcid: Optional[str] = None,
    only_de: bool = False,
    top: Optional[int] = None,
    metric: str = "log2_tpm",
    deg_marker: str = "both",
    output_dir: Optional[str] = None,
    heatmap_rows: str = "sample",
    alpha: float = 0.05,
    max_cgc: int = 500,
    force: bool = False,
    threads: int = 1,
    **kwargs,
) -> List[str]:
    """Plot one or all CGCs; returns list of output PDF paths."""
    if kwargs.get("show_error_bars") or kwargs.get("replicate_groups"):
        logger.warning("show_error_bars/replicate_groups are deprecated; using SeMa-Trap condition means")

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

    outputs: List[str] = []
    task_args = [
        (cid, dbcan_dir, str(expr), str(out), metric, deg_marker, heatmap_rows, alpha)
        for cid in ids
    ]
    if threads <= 1 or len(ids) == 1:
        for t in task_args:
            r = _plot_one_cgc(t)
            if r:
                outputs.extend(r)
    else:
        with ThreadPoolExecutor(max_workers=threads) as pool:
            futs = {pool.submit(_plot_one_cgc, t): t[0] for t in task_args}
            for fut in as_completed(futs):
                r = fut.result()
                if r:
                    outputs.extend(r)
    logger.info(f"Generated {len(outputs)} CGC expression plot files in {out}")
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
    """CGC-level summary heatmap (all CGCs)."""
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
