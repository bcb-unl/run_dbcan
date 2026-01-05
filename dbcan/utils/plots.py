#!/usr/bin/env python3
import os
import sys
import logging
from types import SimpleNamespace

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import seaborn as sns
import rich_click as click

from matplotlib import pyplot
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Polygon

from dbcan.utils.utils import cgc_standard_line
from dbcan.plot.syntenic_plot import (
    syntenic_plot,
    read_blast_result_cgc,
    read_UHGG_CGC_stanrdard_out,
    read_cgcgff,
    Get_parameters_for_plot,
    plot_Polygon_homologous,
    plot_syntenic_block,
)
from dbcan.plot.syntenic_plot import Get_Position as synGet_Position
from dbcan.plot.syntenic_plot import plot_genome_line as synplot_genome_line
from dbcan.constants.plots_constants import (
    CGC_RESULT_FILE,
    CGC_SUB_PREDICTION_FILE,
)

from dbcan.configs.plots_config import PlotsConfig
from dbcan.parameter import logging_options
from dbcan.main import setup_logging

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42

logger = logging.getLogger(__name__)


class CGC_Standard_Out(object):
    def __init__(self, filename):
        hits = open(filename).readlines()[1:]
        self.genes = []
        for line in hits:
            if line.startswith("CGC#"):
                continue
            lines = line.rstrip("\n").split()  #
            gene_obj = cgc_standard_line(lines)
            # newly added
            if len(lines) >= 1:
                setattr(gene_obj, "annotation", lines[-1])
            self.genes.append(gene_obj)

    def __iter__(self):
        return iter(self.genes)

    def CGCID2genes(self):
        cgcdict = {}
        for gene in self:
            cgcdict.setdefault(gene.cgcid, []).append(gene)
        return cgcdict


class CGC(object):
    def __init__(self, genes):
        self.genes = genes
        self.ID = genes[0].cgcid
        self.start = min([gene.gene_start for gene in genes])
        self.end = max([gene.gene_end for gene in genes])
        self.gene_num = len(genes)

    def __iter__(self):
        return iter(self.genes)

    def __repr__(self):
        return "\t".join([self.ID, str(self.start), str(self.end), str(self.gene_num)])

    def __len__(self):
        return len(self.genes)

    def get_positions(self):
        starts = []
        ends = []
        strands = []
        for gene in self:
            starts.append(gene.gene_start)
            ends.append(gene.gene_end)
            strands.append(gene.strand)
        return starts, ends, strands

    def get_proteinID(self):
        return [gene.seqid for gene in self]

    def get_cgc_CAZyme(self):
        return [gene.gene_type for gene in self]

    def get_annotations(self):
        annos = []
        for g in self.genes:
            ann = getattr(g, "annotation", None)
            if not ann:
                # use gene_type as backup
                ann = getattr(g, "gene_type", "")
            annos.append(ann)
        return annos


class CGC_standard_out_2CGC(object):
    def __init__(self, dbcan):
        self.CGCs = []
        cgcdict = dbcan.CGCID2genes()
        for cgc in cgcdict:
            self.CGCs.append(CGC(cgcdict[cgc]))

    def __iter__(self):
        return iter(self.CGCs)

    def cgcid2CGC(self):
        return {cgc.ID: cgc for cgc in self}


def derive_paths(cfg: PlotsConfig):
    return {
        "pul_annotation": os.path.join(cfg.input_dir, CGC_RESULT_FILE),
        "pul_substrate": os.path.join(cfg.input_dir, CGC_SUB_PREDICTION_FILE),
        "blastp": os.path.join(cfg.input_dir, "PUL_blast.out"),
    }


def CGC_plot(cfg: PlotsConfig):
    paths = derive_paths(cfg)
    pul_ann = paths["pul_annotation"]
    if not os.path.exists(pul_ann):
        logger.error(f"CGC annotation file not found: {pul_ann}")
        return
    if not cfg.cgcid:
        logger.error("CGC id is required for CGC_plot (--cgcid).")
        return

    dbCAN_standard_out = CGC_Standard_Out(pul_ann)
    cgcs = CGC_standard_out_2CGC(dbCAN_standard_out)
    cgcid2cgc = cgcs.cgcid2CGC()
    if cfg.cgcid not in cgcid2cgc:
        logger.error(f"CGC id not found in annotation: {cfg.cgcid}")
        return
    cgc = cgcid2cgc[cfg.cgcid]
    starts, ends, strands = cgc.get_positions()
    types = cgc.get_cgc_CAZyme()
    labels = cgc.get_proteinID()
    annotations = cgc.get_annotations() if getattr(cfg, "show_annotation", False) else None  # only when enabled

    out_pdf = f"{cfg.cgcid.replace('|', '_')}.cgc.pdf"
    cgc_fig_plot(starts, ends, strands, types, labels, out_pdf, annotations=annotations)


def read_location_reads_count(filename):
    xs2ys = {}
    with open(filename) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                xs2ys[int(parts[1])] = int(parts[2])
            except Exception:
                continue
    return xs2ys


def CGC_plot_reads_count(cfg: PlotsConfig):
    paths = derive_paths(cfg)
    pul_ann = paths["pul_annotation"]
    if not os.path.exists(pul_ann):
        logger.error(f"CGC annotation file not found: {pul_ann}")
        return
    if not cfg.cgcid:
        logger.error("CGC id is required for CGC_coverage_plot (--cgcid).")
        return
    if not cfg.reads_count or not os.path.exists(cfg.reads_count):
        logger.error("Reads count file is required (--reads-count) and must exist.")
        return

    dbCAN_standard_out = CGC_Standard_Out(pul_ann)
    cgcs = CGC_standard_out_2CGC(dbCAN_standard_out)
    cgcid2cgc = cgcs.cgcid2CGC()
    if cfg.cgcid not in cgcid2cgc:
        logger.error(f"CGC id not found in annotation: {cfg.cgcid}")
        return
    cgc = cgcid2cgc[cfg.cgcid]
    starts, ends, strands = cgc.get_positions()
    types = cgc.get_cgc_CAZyme()
    labels = cgc.get_proteinID()

    out_pdf = f"{cfg.cgcid.replace('|', '_')}.cgc-coverage.pdf"
    cgc_fig_plot_abund(starts, ends, strands, types, labels, cfg, out_pdf)


def Get_Position(starts, ends, strands, labels, yshift=0):
    Width = 1000
    Height = 160
    poly_heigth = 10
    Triangle_length = 4
    plot_start_y = Height / 2 - poly_heigth - yshift

    # copy, do not mutate caller inputs
    s = list(starts)
    e = list(ends)
    st = list(strands)
    shift_pos = min(s) if s else 0
    maxbp = max(e) - min(s) if s and e else 1
    pixeachbp = Width / maxbp if maxbp > 0 else 1.0

    for i in range(len(s)):
        s[i] -= shift_pos
        e[i] -= shift_pos

    lines = []
    polygens = []
    texts = []
    for i in range(len(s)):
        if st[i] == "+":
            positions_str = f"{s[i] * pixeachbp} {plot_start_y} "
            positions_str += f"{e[i] * pixeachbp - Triangle_length} {plot_start_y} "
            positions_str += f"{e[i] * pixeachbp} {plot_start_y + poly_heigth} "
            positions_str += f"{e[i] * pixeachbp - Triangle_length} {plot_start_y + 2*poly_heigth} "
            positions_str += f"{s[i] * pixeachbp} {plot_start_y + 2*poly_heigth}"
        else:
            positions_str = f"{s[i] * pixeachbp} {plot_start_y + poly_heigth} "
            positions_str += f"{s[i] * pixeachbp + Triangle_length} {plot_start_y} "
            positions_str += f"{e[i] * pixeachbp} {plot_start_y} "
            positions_str += f"{e[i] * pixeachbp} {plot_start_y + 2* poly_heigth} "
            positions_str += f"{s[i]* pixeachbp + Triangle_length} {plot_start_y + 2* poly_heigth}"
        polygens.append(positions_str)

        if i < len(s) - 1:
            positions_str = f"{e[i] * pixeachbp} {plot_start_y + poly_heigth} "
            positions_str += f"{s[i+1]*pixeachbp} {plot_start_y + poly_heigth}"
            lines.append(positions_str)

        if i < len(labels):
            texts.append(labels[i].split(".")[0])
        else:
            texts.append(str(i + 1))

    scale_number = 10
    each_scale_bp = maxbp / scale_number
    each_scale_pix = each_scale_bp * pixeachbp

    plot_start_y -= 50
    scale_positions = []
    scale_positions_texts = []
    scale_text = []
    scale_positions.append(f"0 {plot_start_y + 3*poly_heigth} {10*each_scale_pix} {plot_start_y + 3*poly_heigth}")
    plot_start_y -= 1
    for i in range(scale_number + 1):
        positions_str = f"{i*each_scale_pix} {plot_start_y + 3* poly_heigth} "
        positions_str += f"{i*each_scale_pix} {plot_start_y + 3*poly_heigth + 0.6* poly_heigth}"
        scale_positions.append(positions_str)
        positions_str = f"{i*each_scale_pix} {plot_start_y + 3*poly_heigth + 0.6* poly_heigth}"
        scale_positions_texts.append(positions_str)
        scale_text.append(str(int(each_scale_bp * i) + shift_pos))

    return polygens, lines, texts, scale_positions, scale_text


def plot_Polygon(polygens1, types1, ax):
    colors_map = {
        "CAZyme": "#E67E22",
        "TC": "#2ECC71",
        "TF": "#9B59B6",
        "STP": "#F1C40F",
        "Peptidase": "#16A085",
        "Sulfatase": "#34495E",
        "Other": "#95A5A6",
    }
    default_color = "#95A5A6"
    for j in range(len(polygens1)):
        polygen = polygens1[j].split()
        points = []
        color = colors_map.get(types1[j], default_color)
        for i in range(int(len(polygen) / 2)):
            points.append([float(polygen[2 * i]), float(polygen[2 * i + 1])])
        ax.add_patch(Polygon(points, facecolor=color, edgecolor="none", alpha=0.5, lw=0))


def plot_genome_line(lines, ax):
    # draw connectors between adjacent genes as real line segments
    for line in lines:
        x1, y1, x2, y2 = points2(line)
        ax.plot([x1, x2], [y1, y2], color="gray", lw=2, solid_capstyle="round", zorder=5)


def plot_scale_line(lines, label, ax):
    # draw the horizontal scale and tick marks as lines (not polygons)
    for i, line in enumerate(lines):
        x1, y1, x2, y2 = points2(line)
        ax.plot([x1, x2], [y1, y2], color="gray", lw=2, zorder=5)
        if i >= 1:
            # smaller offset so labels stay close to the lines
            ax.text(float(x1), float(y1) - 6, label[i - 1], va="top", ha="center", fontsize=8)


def points2(coord):
    x1, y1, x2, y2 = coord.split()
    return float(x1), float(y1), float(x2), float(y2)


def cgc_fig_plot(starts, ends, strands, types, gene_labels, out_pdf: str, annotations=None):
    genelabelcolor = ["#E67E22", "#2ECC71", "#9B59B6", "#F1C40F", "#16A085", "#34495E", "#95A5A6"]
    geneslabels = ["CAZyme", "TC", "TF", "STP", "Peptidase", "Sulfatase", "Other"]
    genecustom_lines = [Patch(color=c, alpha=0.5) for c in genelabelcolor]

    px = 1 / plt.rcParams["figure.dpi"]
    Width = 1400
    Height = 100
    fig = plt.figure(figsize=(Width * px * 1.2, Height * px * 2))
    ax = fig.add_subplot(111)

    polygens, lines, texts, scale_positions, scale_text = Get_Position(starts, ends, strands, gene_labels)
    plot_Polygon(polygens, types, ax)
    plot_genome_line(lines, ax)
    plot_scale_line(scale_positions, scale_text, ax)

    ylim_top = 150
    if annotations:
        # rotation 60, larger font; dynamically extend ylim based on returned top y
        top_y = add_gene_annotations(
            ax, starts, ends, annotations,
            font_size=8,
            rotation=60,
            y_offset=28
        )
        if top_y:
            ylim_top = max(ylim_top, top_y + 5)

    ax.plot()
    legend = pyplot.legend(genecustom_lines, geneslabels, frameon=False, loc="best", title_fontsize="x-large")
    ax.add_artist(legend)
    plt.ylim(0, ylim_top)
    plt.xlim(-50, 1100)
    # Remove tight_layout to avoid shrinking arrows when text grows
    fig.subplots_adjust(left=0.02, right=0.98, top=0.96, bottom=0.02)
    plt.axis("off")
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.close()
    logger.info(f"Saved figure: {out_pdf}")


def add_gene_annotations(ax, starts, ends, annotations, font_size=10, rotation=60, y_offset=25):
    """
    Add full annotation text above gene arrows (no truncation).
    rotation: text rotation angle (default 60).
    font_size: text size.
    y_offset: vertical offset above arrow top in data coordinates.
    Returns the highest y used by text for dynamic ylim adjustment.
    """
    if not starts or not ends:
        return
    shift_pos = min(starts)
    maxbp = max(ends) - shift_pos
    if maxbp <= 0:
        return
    Width = 1000  # must match Get_Position
    pixeachbp = Width / maxbp
    # Arrow top baseline (from Get_Position geometry):
    # Height/2 - poly_height = 160/2 - 10 = 70; arrow top ~ 70 + 20 = 90
    base_polygon_y_top = 70 + 20
    label_base_y = base_polygon_y_top + y_offset
    max_top_y = label_base_y
    # approximate single-line height in data units (rough factor)
    line_height = font_size * 0.62
    for s, e, ann in zip(starts, ends, annotations):
        try:
            cx = ((s - shift_pos) + (e - shift_pos)) / 2 * pixeachbp
            txt = ax.text(
                cx,
                label_base_y,
                ann,
                ha="center",
                va="bottom",
                fontsize=font_size,
                rotation=rotation,
            )
            top_est = label_base_y + line_height  # single line
            if top_est > max_top_y:
                max_top_y = top_est
        except Exception:
            continue
    return max_top_y


def cgc_fig_plot_abund(starts, ends, strands, types, labels, cfg: PlotsConfig, out_pdf: str):
    genelabelcolor = ["#E67E22", "#2ECC71", "#9B59B6", "#F1C40F", "#16A085", "#34495E", "#95A5A6"]
    geneslabels = ["CAZyme", "TC", "TF", "STP", "Peptidase", "Sulfatase", "Other"]
    genecustom_lines = [Patch(color=c, alpha=0.5) for c in genelabelcolor]

    px = 1 / plt.rcParams["figure.dpi"]
    Width = 1400
    Height = 100
    fig = plt.figure(figsize=(Width * px * 1.2, Height * px * 4))
    ax = fig.add_subplot(212)

    polygens, lines, texts, scale_positions, scale_text = Get_Position(starts, ends, strands, labels)
    plot_Polygon(polygens, types, ax)
    plot_genome_line(lines, ax)
    plot_scale_line(scale_positions, scale_text, ax)
    ax.plot()
    legend = pyplot.legend(genecustom_lines, geneslabels, frameon=False, loc="best", title_fontsize="x-large")
    ax.add_artist(legend)
    plt.ylim(0, 150)
    xlim_x1, xlim_x2 = (-10, 1100)
    plt.xlim(xlim_x1, xlim_x2)
    plt.axis("off")

    xs2ys = read_location_reads_count(cfg.reads_count)
    max_y = max(xs2ys.values()) if xs2ys else 0
    add_readcount_layout(fig, starts, ends, xs2ys, max_y, -3, max_y + 10, xlim_x1, xlim_x2, max(ends) - min(starts))

    plt.savefig(out_pdf, bbox_inches="tight")
    plt.close()
    logger.info(f"Saved figure: {out_pdf}")


def add_readcount_layout(fig, starts, ends, xs2ys, max_y, ylim_y1, ylim_y2, xlim_x1, xlim_x2, syn_maxbp):
    maxbp = max(ends) - min(starts)
    Width = 1000
    pixeachbp = Width / syn_maxbp if syn_maxbp > 0 else 1.0
    ax = fig.add_subplot(211)
    plt.ylim(ylim_y1, ylim_y2)
    plt.xlim(xlim_x1, xlim_x2)
    plt.tight_layout(pad=0.1)
    plt.plot((0, 1000), (0, 0), color="gray", lw=1)
    all_xs = []
    all_ys = []
    start = min(starts)
    for i in range(1, maxbp + 1):
        all_xs.append(pixeachbp * i)
        all_ys.append(xs2ys.get(i + start, 0))
    plt.plot(all_xs, all_ys, "-", alpha=0.5, color="red", lw=1)
    ax.fill_between(all_xs, all_ys, 0, facecolor="red", alpha=0.3, edgecolor="white")
    for pos in ["top", "right", "bottom"]:
        ax.spines[pos].set_visible(False)
    ax.tick_params(bottom=False, top=False, left=True, right=False)
    ax.set_xticks([])


def generate_syntenic_block(cgcpul, cgcpul_blastp, genes1, genes2):
    blocks = []
    for record in cgcpul_blastp.get(cgcpul, []):
        query = record.qseqid
        hit = record.sseqid
        try:
            cgc_proteinid = query.split("|")[2]
        except Exception:
            continue
        # sseqid shape often contains PUL and protein id separated by ':'
        parts = hit.split(":")
        pul_proteinid = parts[2] if len(parts) > 2 else (parts[-1] if parts else "")
        try:
            index1 = genes1.index(cgc_proteinid)
            index2 = genes2.index(pul_proteinid)
            blocks.append(f"{index1}-{index2}-{record.pident}")
        except Exception:
            continue
    return blocks


def add_synteny_scale(ax, starts_cgc, ends_cgc, starts_pul, ends_pul, maxbp,
                      ticks=10, width=1000, y=4, font_size=7):
    """
    Draw a manual horizontal scale (relative length; each track internally shifted).
    Labels are absolute positions using the minimal start of both tracks as origin.
    """
    if maxbp <= 0:
        return
    try:
        min_abs = min(
            [m for m in (min(starts_cgc, default=0), min(starts_pul, default=0))]
        )
    except Exception:
        min_abs = 0
    pix_per_bp = width / maxbp
    # baseline
    ax.plot([0, width], [y, y], color="gray", lw=1)
    step_bp = maxbp / ticks
    for i in range(ticks + 1):
        x = pix_per_bp * step_bp * i
        ax.plot([x, x], [y, y + 4], color="gray", lw=1)
        # avoid overcrowding: only label min, max, and every second tick
        if i in (0, ticks) or i % 2 == 0:
            label_val = int(min_abs + step_bp * i)
            ax.text(x, y - 6, f"{label_val}", ha="center", va="top", fontsize=font_size)
    ax.text(width, y - 16, "bp", ha="right", va="top", fontsize=font_size)


def CGC_syntenic_with_PUL(cfg: PlotsConfig):
    paths = derive_paths(cfg)
    sub_pred = paths["pul_substrate"]
    pul_ann = paths["pul_annotation"]
    if not os.path.exists(sub_pred) or not os.path.exists(pul_ann):
        logger.error("Required files not found (substrate_prediction or cgc_standard_out).")
        return
    if not cfg.cgcid:
        logger.error("CGC id is required (--cgcid).")
        return

    cgcid2pulid = {line.rstrip().split("\t")[0]: line.rstrip().split("\t")[1] for line in open(sub_pred).readlines()[1:]}
    cgc = cfg.cgcid
    pul = cgcid2pulid.get(cgc, "")
    if not pul:
        logger.error(f"Homolog PUL not found for CGC: {cgc}")
        return

    cgcpul_blastp = read_blast_result_cgc(paths["blastp"])
    cgc_proteinid2gene, cgcid2gene, cgcid2geneid = read_UHGG_CGC_stanrdard_out(pul_ann)
    ns = SimpleNamespace(db_dir=cfg.db_dir)
    PULid_proteinid2gene, PULid2gene, PULid2geneid = read_cgcgff(ns)

    cgcpul = cgc + ":" + pul
    bed_cgc = cgcid2gene[cgc]
    bed_pul = PULid2gene[pul]
    starts1, ends1, strands1, types1 = Get_parameters_for_plot(bed_cgc)
    starts2, ends2, strands2, types2 = Get_parameters_for_plot(bed_pul)
    genes1 = cgcid2geneid[cgc]
    genes2 = PULid2geneid[pul]
    blocks = generate_syntenic_block(cgcpul, cgcpul_blastp, genes1, genes2)
    config = {"output_dir": "."}
    syntenic_plot(starts1, starts2, ends1, ends2, strands1, strands2, types1, types2, blocks, cgc, pul, config)
    logger.info(f"Saved figure: {cgc.replace('|','_')}_{pul.replace('|','_')}-syntenic.pdf")


def CGC_syntenic_with_PUL_abund(cfg: PlotsConfig):
    paths = derive_paths(cfg)
    sub_pred = paths["pul_substrate"]
    pul_ann = paths["pul_annotation"]
    if not os.path.exists(sub_pred) or not os.path.exists(pul_ann):
        logger.error("Required files not found (substrate_prediction or cgc_standard_out).")
        return
    if not cfg.cgcid:
        logger.error("CGC id is required (--cgcid).")
        return
    if not cfg.reads_count or not os.path.exists(cfg.reads_count):
        logger.error("Reads count file is required (--reads-count) and must exist.")
        return

    cgcid2pulid = {line.rstrip().split("\t")[0]: line.rstrip().split("\t")[1] for line in open(sub_pred).readlines()[1:]}
    cgc = cfg.cgcid
    pul = cgcid2pulid.get(cgc, "")
    if not pul:
        logger.error(f"Homolog PUL not found for CGC: {cgc}")
        return

    cgcpul_blastp = read_blast_result_cgc(paths["blastp"])
    cgc_proteinid2gene, cgcid2gene, cgcid2geneid = read_UHGG_CGC_stanrdard_out(pul_ann)
    # Load PUL GFFs (wrapper around read_cgcgff which expects a file + dict)
    PULid_proteinid2gene, PULid2gene, PULid2geneid = load_pul_gffs(cfg.db_dir)

    cgcpul = cgc + ":" + pul
    bed_cgc = cgcid2gene[cgc]
    bed_pul = PULid2gene[pul]
    starts1, ends1, strands1, types1 = Get_parameters_for_plot(bed_cgc)
    starts2, ends2, strands2, types2 = Get_parameters_for_plot(bed_pul)
    genes1 = cgcid2geneid[cgc]
    genes2 = PULid2geneid[pul]
    blocks = generate_syntenic_block(cgcpul, cgcpul_blastp, genes1, genes2)

    # draw synteny and coverage
    px = 1 / plt.rcParams["figure.dpi"]
    Width = 1600
    Height = 620 * 2
    fig = plt.figure(figsize=(Width * px, Height * px * 2 / 2.5))
    ax = fig.add_subplot(212)

    maxbp = max([max(ends1) - min(starts1), max(ends2) - min(starts2)])
    # First (CGC) band – reduce yshift to move upward (more compact)
    polygens_cgc, blocks_cgc_coor, lines_cgc_coor, _, _ = synGet_Position(
        starts1, ends1, strands1, maxbp, yshift=30, up=2
    )
    # Second (PUL) band
    polygens_pul, blocks_pul_coor, lines_pul_coor, _, _ = synGet_Position(
        starts2, ends2, strands2, maxbp, yshift=0, up=1
    )
    # draw bands, blocks, lines
    plot_Polygon_homologous(polygens_cgc, polygens_pul, types1, types2, 2, ax)
    plot_syntenic_block(blocks, blocks_cgc_coor, blocks_pul_coor, ax)
    synplot_genome_line(lines_cgc_coor, lines_pul_coor, ax)

    custom_lines = [
        Line2D([0], [0], color="red", lw=4, alpha=0.5),
        Line2D([0], [0], color="blue", lw=4, alpha=0.5),
        Line2D([0], [0], color="green", lw=4, alpha=0.5),
        Line2D([0], [0], color="cyan", lw=4, alpha=0.5),
        Line2D([0], [0], color="gray", lw=4, alpha=0.5),
    ]
    labels = ["80-100", "60-80", "40-60", "20-40", "0-20"]

    genelabelcolor = ["#E67E22", "#2ECC71", "#9B59B6", "#F1C40F", "#16A085", "#34495E", "#95A5A6"]
    geneslabels = ["CAZyme", "TC", "TF", "STP", "Peptidase", "Sulfatase", "Other"]
    genecustom_lines = [Patch(color=c, alpha=0.5) for c in genelabelcolor]

    legend1 = pyplot.legend(custom_lines, labels, frameon=False, loc="lower right", bbox_to_anchor=(1, 0.5), title="Identity")
    ax.add_artist(legend1)
    legend2 = pyplot.legend(genecustom_lines, geneslabels, frameon=False, loc="lower right", bbox_to_anchor=(1, 0.1), title="Gene")
    ax.add_artist(legend2)

    # Move CGC title upward (was y=10). Slightly adjust PUL title for compact layout.
    plt.text(500, 15, cgc, fontsize=20, horizontalalignment="center")
    plt.text(500, 88, pul, fontsize=20, horizontalalignment="center")
    xlim_x1, xlim_x2 = (-10, 1100)
    ylim_y1, ylim_y2 = (0, 100)
    plt.ylim(ylim_y1, ylim_y2)
    plt.xlim(xlim_x1, xlim_x2)
    # Custom scale (syntenic Get_Position does not provide one)
    add_synteny_scale(ax, starts1, ends1, starts2, ends2, maxbp)
    # Hide default axes
    plt.axis("off")
    ax.plot()

    xs2ys = read_location_reads_count(cfg.reads_count)
    max_y = max(xs2ys.values()) if xs2ys else 0
    add_readcount_layout(fig, starts1, ends1, xs2ys, max_y, ylim_y1, max_y, xlim_x1, xlim_x2, maxbp)

    out_pdf = f"{cgc.replace('|','_')}-syntenic-cov.pdf"
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.close()
    logger.info(f"Saved figure: {out_pdf}")


def combined_datafram_based_on_first_col(pd_lists, samples):
    if len(pd_lists) <= 1:
        return pd_lists[0]
    else:
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
    merged_table.fillna(0, inplace=True)
    merged_table["diff_abs"] = np.abs(merged_table[abundance_col + "_x"] - merged_table[abundance_col + "_y"])
    merged_table["diff"] = merged_table[abundance_col + "_x"] - merged_table[abundance_col + "_y"]

    merged_columns = merged_table.columns
    rename_columns = []
    abund_index = 0
    for column in merged_columns:
        if abundance_col in column and column != on_merge_col:
            rename_columns.append(samples[abund_index])
            abund_index += 1
        else:
            rename_columns.append(column)
    merged_table.columns = rename_columns
    merged_table.sort_values("diff_abs", inplace=True, ascending=False)
    return merged_table, on_merge_col


def filter_out_enzyme_number(table):
    bools = []
    for i in table.iloc[:, 0]:
        if i and i[0].isdigit():
            bools.append(True)
        elif i in ["PL0", "GH0", "GT0", "CBM0", "AA0", "CE0"]:
            bools.append(False)
        else:
            bools.append(True)
    table = table[bools]
    return table


def add_column_type(table):
    import re

    cols = []
    for i in table["CAZy"]:
        fam = re.sub(r"[0-9]+", "", i)
        cols.append(fam)
    table["fam"] = cols
    return table


def check_input_files_min_lines(input_files, min_lines=2) -> bool:
    for f in input_files.split(","):
        try:
            with open(f) as fh:
                lines = sum(1 for _ in fh)
            if lines < min_lines:
                logger.warning(f"Input file {f} has less than {min_lines} lines, skip plotting.")
                return False
        except Exception as e:
            logger.warning(f"Error reading file {f}: {e}")
            return False
    return True


def heatmap_plot(cfg: PlotsConfig, input_files: str, samples: str, show_abund: bool, cluster_map: bool, palette: str, col: str, value: str, plot_style: str, top: int, show_fig: bool):
    pds = [filter_out_enzyme_number(pd.read_csv(filename, sep="\t")) for filename in input_files.split(",")]
    samples_list = samples.split(",")
    plt.style.use(plot_style)
    if len(pds) != len(samples_list):
        logger.error("The number of samples is not equal to the number of input files.")
        return
    for i in range(len(pds)):
        pds[i]["sample"] = samples_list[i]
    if len(pds) == 1:
        data = pds[0]
        data = data.rename(columns={data.columns[1]: samples_list[0]})
    else:
        data, _ = combined_datafram_based_on_first_col(pds, samples_list)

    if not col:
        data = data.iloc[0 : int(top), :]
    else:
        if value:
            data = data.loc[data[col].isin(value.split(","))]
        else:
            data = data.iloc[0 : int(top), :]
    data = data.set_index(data.iloc[:, 0])
    data = data[samples_list]
    sns.set_style("whitegrid")
    sns.set_context("paper")

    n_rows, n_cols = data.shape
    cell_width = 0.7
    cell_height = 0.5
    min_width = 6
    min_height = 4
    max_width = 30
    max_height = 30
    fig_width = min(max(n_cols * cell_width, min_width), max_width)
    fig_height = min(max(n_rows * cell_height, min_height), max_height)

    if palette:
        cmap = palette
    else:
        mycolor = ["aliceblue", "skyblue", "deepskyblue", "orange", "tomato", "red"]
        cmap = colors.LinearSegmentedColormap.from_list("my_list", mycolor)
        cmap.set_under("white")

    if cluster_map:
        g = sns.clustermap(
            data,
            cmap=cmap,
            cbar=True,
            vmin=0.1,
            dendrogram_ratio=0.03,
            cbar_pos=(0.1, 1, 0.1, 0.1),
            col_cluster=False,
            cbar_kws={"shrink": 0.3},
            figsize=(fig_width, fig_height),
        )
        plt.setp(g.ax_heatmap.get_xticklabels(), rotation=30, ha="right", fontsize=10)
        plt.setp(g.ax_heatmap.get_yticklabels(), fontsize=10)
        if show_fig:
            plt.show()
        else:
            plt.savefig("heatmap_cluster.pdf", bbox_inches="tight")
    else:
        plt.figure(figsize=(fig_width, fig_height))
        ax = sns.heatmap(
            data,
            cmap=cmap,
            yticklabels=True,
            annot=show_abund,
            fmt=".0f",
            linewidths=0 if not show_abund else 0.5,
            cbar=True,
            vmin=0.1,
            cbar_kws={"shrink": 0.3, "anchor": (0, 0.0)},
        )
        ax.collections[0].colorbar.ax.tick_params(labelsize=8)
        plt.xticks(rotation=30, ha="right", fontsize=10)
        plt.yticks(fontsize=10)
        plt.tight_layout(pad=0.5)
        if show_fig:
            plt.show()
        else:
            plt.savefig("heatmap.pdf", bbox_inches="tight")


def bar_plot(cfg: PlotsConfig, input_files: str, samples: str, plot_style: str, col: str, value: str, top: int, vertical_bar: bool, pdf: str):
    pds = [filter_out_enzyme_number(pd.read_csv(filename, sep="\t")) for filename in input_files.split(",")]
    samples_list = samples.split(",")
    plt.style.use(plot_style)
    if len(pds) != len(samples_list):
        logger.error("The number of samples is not equal to the number of input files.")
        return
    for i in range(len(pds)):
        pds[i]["sample"] = samples_list[i]
    result = combined_datafram_based_on_first_col(pds, samples_list)
    if isinstance(result, tuple):
        data, x = result
    else:
        data = result
        x = data.columns[0]
        if len(samples_list) == 1:
            old_abund_col = data.columns[1]
            data = data.rename(columns={old_abund_col: samples_list[0]})

    if not col:
        data = data.iloc[0 : int(top), :]
    else:
        if value:
            data = data.loc[data[col].isin(value.split(","))]
        else:
            data = data.iloc[0 : int(top), :]

    if vertical_bar:
        ax = data.plot.barh(x=x, y=samples_list)
        plt.ylabel("")
        plt.xlabel("Abundance")
    else:
        ax = data.plot(x=x, y=samples_list, kind="bar")
        plt.xticks(rotation=90)
        plt.xlabel("")
        plt.ylabel("Abundance")

    plt.title(f"The most top {top} different families")
    if not pdf.endswith(".pdf"):
        pdf = pdf + ".pdf"
    plt.savefig(pdf, bbox_inches="tight")
    logger.info(f"Saved plot: {pdf}")


click.rich_click.USE_RICH_MARKUP = True
click.rich_click.USE_MARKDOWN = True
click.rich_click.SHOW_ARGUMENTS = True
click.rich_click.GROUP_ARGUMENTS_OPTIONS = True
click.rich_click.STYLE_ERRORS_SUGGESTION = "magenta italic"
click.rich_click.SHOW_METAVARS_COLUMN = False
click.rich_click.APPEND_METAVARS_HELP = True
click.rich_click.MAX_WIDTH = 100
click.rich_click.COMMAND_GROUPS = {
    "dbcan_plot": [
        {
            "name": "CGC/PUL Visualization",
            "commands": ["CGC_plot", "CGC_coverage_plot", "CGC_synteny_plot", "CGC_synteny_coverage_plot"],
        },
        {
            "name": "Abundance Visualization",
            "commands": ["bar_plot", "heatmap_plot"],
        },
    ]
}

@click.command()
@logging_options
@click.argument(
    "function", 
    type=click.Choice([
        "CGC_plot", "CGC_coverage_plot", "CGC_synteny_plot", "CGC_synteny_coverage_plot",
        "bar_plot", "heatmap_plot"
    ])
)
@click.option("-i", "--input", "input_path", required=True, 
              help="dbCAN CAZyme annotation output folder or input files (bar/heatmap).")
@click.option("--db-dir", "--db_dir", default="db", show_default=True, 
              help="Database directory.")
@click.option("--cgcid", "cgcid", default=None, 
              help="CGC id (contig_ID|cgc_order).")
@click.option("--reads-count", "--reads_count", "--readscount", default=None, 
              help="Read counts file generated by samtools depth.")
@click.option("--samples", "samples", default=None, 
              help='Sample names separated by "," (bar/heatmap).')
@click.option("--top", "top", default=20, show_default=True, type=int, 
              help="Top N entries to plot.")
@click.option("--plot-style", "--plot_style", default="ggplot", show_default=True, 
              type=click.Choice(["ggplot", "seaborn", "seaborn-poster"]), 
              help="Style for barplot and heatmap.")
@click.option("--vertical-bar", "--vertical_bar", is_flag=True, 
              help="Use horizontal bar (True) vs vertical (False).")
@click.option("--show-fig", "--show_fig", is_flag=True, 
              help="Show figure interactively.")
@click.option("--show-abund", "--show_abund", is_flag=True, 
              help="Show abundance values on heatmap cells.")
@click.option("--palette", "palette", default=None, 
              help="Matplotlib colormap/palette name.")
@click.option("--cluster-map", "--cluster_map", is_flag=True, 
              help="Clustered heatmap.")
@click.option("--col", "col", default=None, 
              help="Filter column name.")
@click.option("--value", "value", default=None, 
              help='Filter values separated by ",".')
@click.option("--pdf", "pdf", default="bar_plot.pdf", show_default=True, 
              help="Bar plot output pdf file.")
@click.option("--show-annotation", "--show_annotation", is_flag=True,
              help="Show gene annotation labels in CGC_plot (default: hidden).")
def cli(function, log_level, log_file, verbose, input_path, db_dir, cgcid, reads_count, samples, top, plot_style, vertical_bar, show_fig, show_abund, palette, cluster_map, col, value, pdf, show_annotation):
    """# dbCAN plotting utilities
    
    ## CGC/PUL Visualization Commands
    
    **CGC_plot**: Generate basic CGC gene visualization
    ```
    dbcan_plot CGC_plot -i ./sample.dbCAN --cgcid 'contigX|CGC1'
    ```
    
    **CGC_coverage_plot**: Generate CGC visualization with read coverage
    ```
    dbcan_plot CGC_coverage_plot -i ./sample.dbCAN --cgcid 'contigX|CGC1' --reads-count cgc.depth.txt
    ```
    
    **CGC_synteny_plot**: Generate CGC-PUL synteny comparison
    ```
    dbcan_plot CGC_synteny_plot -i ./sample.dbCAN --cgcid 'contigX|CGC1'
    ```
    
    **CGC_synteny_coverage_plot**: Generate CGC-PUL comparison with coverage
    ```
    dbcan_plot CGC_synteny_coverage_plot -i ./sample.dbCAN --cgcid 'contigX|CGC1' --reads-count cgc.depth.txt
    ```
    
    ## Abundance Visualization Commands
    
    **bar_plot**: Create bar plot of CAZyme abundance
    ```
    dbcan_plot bar_plot -i a.CAZyme_abund,b.CAZyme_abund --samples A,B --vertical-bar
    ```
    
    **heatmap_plot**: Create heatmap of CAZyme abundance
    ```
    dbcan_plot heatmap_plot -i a.CAZyme_abund,b.CAZyme_abund --samples A,B --show-abund
    ```
    """
    setup_logging(log_level, log_file, verbose)
    # Build config
    cfg = PlotsConfig(
        input_dir=input_path if os.path.isdir(input_path) else os.path.dirname(os.path.abspath(input_path)) or ".",
        db_dir=db_dir,
        cgcid=cgcid,
        reads_count=reads_count,
        samples=(samples.split(",") if samples else None),
        top=top,
        plot_style=plot_style,
        vertical_bar=vertical_bar,
        show_fig=show_fig,
        show_abund=show_abund,
        palette=palette,
        cluster_map=cluster_map,
        filter_col=col,
        filter_value=value,
        pdf=pdf,
    )
    # attach optional flag even if PlotsConfig does not define it (fallback)
    setattr(cfg, "show_annotation", show_annotation)

    # Dispatch
    if function == "CGC_plot":
        CGC_plot(cfg)
    elif function == "CGC_coverage_plot":
        CGC_plot_reads_count(cfg)
    elif function == "CGC_synteny_plot":
        CGC_syntenic_with_PUL(cfg)
    elif function == "CGC_synteny_coverage_plot":
        CGC_syntenic_with_PUL_abund(cfg)
    elif function == "bar_plot":
        if not check_input_files_min_lines(input_path, min_lines=2):
            return
        if not samples:
            logger.error("--samples is required for bar_plot")
            return
        bar_plot(cfg, input_files=input_path, samples=samples, plot_style=plot_style, col=col, value=value, top=top, vertical_bar=vertical_bar, pdf=pdf)
    elif function == "heatmap_plot":
        if not check_input_files_min_lines(input_path, min_lines=2):
            return
        if not samples:
            logger.error("--samples is required for heatmap_plot")
            return
        heatmap_plot(cfg, input_files=input_path, samples=samples, show_abund=show_abund, cluster_map=cluster_map, palette=palette, col=col, value=value, plot_style=plot_style, top=top, show_fig=show_fig)


if __name__ == "__main__":
    cli()

# Provide a main() entry-point for console_scripts referencing dbcan.utils.plots:main
def main():
    cli()

# ---- helpers ----
def load_pul_gffs(db_dir: str):
    """
    Collect PUL gene models by scanning db_dir/dbCAN-PUL/*/cgc.gff
    Returns:
      proteinid2gene  (key: PULID:protein_id)
      pulid2genes     (key: PULID -> list[gene])
      pulid2geneids   (key: PULID -> list[protein_id])
    """
    base = os.path.join(db_dir, "dbCAN-PUL")
    proteinid2gene = {}
    if not os.path.isdir(base):
        logger.error(f"PUL directory not found: {base}")
        return {}, {}, {}
    from dbcan.plot.syntenic_plot import CGC_stanrdard  # reuse class
    # replicate logic of syntenic_plot.read_PUL_cgcgff
    for entry in os.scandir(base):
        if entry.is_dir() and entry.name.startswith("PUL") and entry.name.endswith(".out"):
            gff_path = os.path.join(entry.path, "cgc.gff")
            if os.path.exists(gff_path):
                read_cgcgff(gff_path, proteinid2gene)
    pulid2genes = {}
    pulid2geneids = {}
    for k, gene in proteinid2gene.items():
        pulid2genes.setdefault(gene.CGCID, []).append(gene)
        pulid2geneids.setdefault(gene.CGCID, []).append(gene.Protein_ID)
    return proteinid2gene, pulid2genes, pulid2geneids
