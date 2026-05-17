"""Orchestrate expression and DEG analysis pipeline."""

from __future__ import annotations

import logging
import os
from pathlib import Path

from dbcan.configs.expression_config import ExpressionConfig
from dbcan.expression.abundance import (
    build_cgc_gene_counts_and_matrix,
    run_per_sample_abundance,
)
from dbcan.expression.alignment import resolve_bam_paths
from dbcan.expression.cgc_de import summarize_cgc_de
from dbcan.expression.deseq2_runner import run_deseq2, write_sample_metadata
from dbcan.expression.gff_resolve import resolve_expression_gff
from dbcan.expression.quantification import build_count_matrix
from dbcan.expression.samplesheet import load_samplesheet

logger = logging.getLogger(__name__)


def _normalize_paths(config: ExpressionConfig) -> None:
    """Resolve paths so downstream chdir() does not break relative paths."""
    config.input_dir = str(Path(config.input_dir).expanduser().resolve())
    config.output_dir = str(Path(config.output_dir).expanduser().resolve())
    if config.gff:
        config.gff = str(Path(config.gff).expanduser().resolve())
    if config.reference_fasta:
        config.reference_fasta = str(Path(config.reference_fasta).expanduser().resolve())
    if config.samplesheet:
        config.samplesheet = str(Path(config.samplesheet).expanduser().resolve())


def _validate_config(config: ExpressionConfig) -> str:
    _normalize_paths(config)
    if not config.input_dir or not Path(config.input_dir).is_dir():
        raise FileNotFoundError(f"run_dbcan output directory not found: {config.input_dir}")
    if not config.samplesheet:
        raise ValueError("--samplesheet is required")
    overview = Path(config.input_dir) / "overview.tsv"
    if not overview.is_file():
        raise FileNotFoundError(f"overview.tsv not found: {overview}")
    gff_path = resolve_expression_gff(config)
    return gff_path


def run_plots_if_requested(config: ExpressionConfig) -> None:
    """Invoke dbcan_plot helpers after pipeline completes."""
    if not config.run_plots:
        return
    from dbcan.plot.expression_plots import (
        plot_deg_cgc_bar,
        plot_deg_volcano,
        run_cgc_expression_plots,
    )

    expr_dir = Path(config.output_dir)
    dbcan_dir = config.input_dir
    plots_dir = expr_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    gene_deg = expr_dir / "gene_deg.tsv"
    cgc_deg = expr_dir / "cgc_deg.tsv"
    if gene_deg.is_file():
        plot_deg_volcano(gene_deg, plots_dir / "DEG_volcano.pdf")
    if cgc_deg.is_file():
        plot_deg_cgc_bar(cgc_deg, plots_dir / "DEG_CGC_bar.pdf")

    run_cgc_expression_plots(
        dbcan_dir=dbcan_dir,
        expression_dir=str(expr_dir),
        cgcid=None,
        only_de=config.run_plots_only_de,
        top=None,
        metric="counts",
        deg_marker="edge",
        output_dir=str(plots_dir / "CGC_expression"),
        max_cgc=config.max_cgc,
        force=config.force,
        threads=config.threads,
    )
    if config.also_heatmap:
        from dbcan.plot.expression_plots import plot_cgc_expression_heatmap
        plot_cgc_expression_heatmap(expr_dir, plots_dir / "CGC_expression_heatmap.pdf")


def run_expression_pipeline(config: ExpressionConfig) -> None:
    """Run full expression workflow."""
    gff_path = _validate_config(config)
    config.gff = gff_path
    out = Path(config.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    samplesheet = load_samplesheet(config.samplesheet)
    write_sample_metadata(config, samplesheet)

    bam_map = resolve_bam_paths(config, samplesheet)
    count_matrix, depth_paths = build_count_matrix(config, samplesheet, bam_map)

    run_per_sample_abundance(config, samplesheet, depth_paths)
    build_cgc_gene_counts_and_matrix(config, samplesheet, depth_paths)

    if not config.skip_deseq2:
        gene_deg = run_deseq2(config, samplesheet, count_matrix)
        summarize_cgc_de(config, samplesheet, count_matrix, gene_deg)
    else:
        logger.info("Skipping DESeq2 (--skip-deseq2)")

    run_plots_if_requested(config)
    logger.info(f"Expression analysis complete: {config.output_dir}")
