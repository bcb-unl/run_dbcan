"""PyDESeq2 differential expression analysis."""

from __future__ import annotations

import logging
import re
from collections import Counter
from pathlib import Path

import pandas as pd

from dbcan.configs.expression_config import ExpressionConfig
from dbcan.constants.expression_constants import DEG_FILE, GENE_DEG_FILE, GENE_NORM_COUNTS_FILE
from dbcan.expression.samplesheet import Samplesheet

logger = logging.getLogger(__name__)


def _require_pydeseq2():
    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
        return DeseqDataSet, DeseqStats
    except ImportError as e:
        raise ImportError(
            "PyDESeq2 is required for differential expression. "
            "Install with: pip install 'dbcan[expression]'"
        ) from e
    except AttributeError as e:
        if "StringDType" in str(e):
            raise ImportError(
                "PyDESeq2/anndata are incompatible with NumPy 2.4+. "
                "Fix with: pip install 'numpy>=2.0,<2.4' "
                "or pip install -U 'dbcan[expression]'"
            ) from e
        raise


def validate_deseq2_design(samplesheet: Samplesheet, design: str) -> None:
    """Fail fast when the design has no residual df for dispersion (PyDESeq2)."""
    design = design.strip()
    if not design.startswith("~"):
        design = f"~{design}"

    # Only validate simple single-factor designs used by default (~condition).
    m = re.match(r"^~([A-Za-z_][A-Za-z0-9_]*)$", design)
    if not m:
        return
    factor = m.group(1)

    if factor == "condition":
        level_counts = Counter(s.condition for s in samplesheet)
    elif factor == "replicate_group":
        level_counts = Counter(
            (s.replicate_group or s.sample_id) for s in samplesheet
        )
    else:
        return

    n_samples = len(samplesheet)
    n_levels = len(level_counts)
    if n_samples <= n_levels:
        detail = ", ".join(f"{k}({v})" for k, v in sorted(level_counts.items()))
        raise ValueError(
            f"DESeq2 cannot run: {n_samples} sample(s) but {n_levels} distinct "
            f"'{factor}' levels ({detail}). PyDESeq2 needs biological replicates "
            f"(fewer levels than samples, with ≥2 samples per compared group). "
            f"Merge conditions in the samplesheet (e.g. control vs treatment) or "
            f"use --skip-deseq2 for abundance/plots only."
        )

    singletons = sorted(k for k, v in level_counts.items() if v < 2)
    if singletons:
        raise ValueError(
            f"DESeq2 requires at least 2 biological replicates per '{factor}' level. "
            f"Levels with only one sample: {', '.join(singletons)}."
        )


def write_gene_norm_counts(config: ExpressionConfig, dds) -> Path:
    """Write DESeq2 size-factor normalized counts (samples x genes)."""
    out = Path(config.output_dir) / GENE_NORM_COUNTS_FILE
    try:
        layer = dds.layers["normed_counts"]
    except (KeyError, AttributeError):
        try:
            layer = dds.normed_counts
        except AttributeError:
            logger.warning("Could not extract normed counts from DESeq2; skipping gene_norm_counts.tsv")
            return out
    normed = pd.DataFrame(
        layer,
        index=dds.obs_names,
        columns=dds.var_names,
    )
    normed.to_csv(out, sep="\t")
    logger.info(f"Wrote {out}")
    return out


def write_sample_metadata(config: ExpressionConfig, samplesheet: Samplesheet) -> Path:
    """Export sample metadata for reproducibility."""
    rows = []
    for s in samplesheet:
        rows.append({
            "sample_id": s.sample_id,
            "condition": s.condition,
            "replicate_group": s.replicate_group or s.sample_id,
        })
    df = pd.DataFrame(rows).set_index("sample_id")
    out = Path(config.output_dir) / "sample_metadata.tsv"
    df.to_csv(out, sep="\t")
    return out


def run_deseq2(
    config: ExpressionConfig,
    samplesheet: Samplesheet,
    count_matrix: pd.DataFrame,
) -> pd.DataFrame:
    """Run gene-level DESeq2 and write gene_deg.tsv and DEG.tsv."""
    DeseqDataSet, DeseqStats = _require_pydeseq2()
    validate_deseq2_design(samplesheet, config.design)

    metadata = pd.DataFrame(
        [{"condition": s.condition} for s in samplesheet],
        index=[s.sample_id for s in samplesheet],
    )
    counts = count_matrix.loc[metadata.index]
    # Drop genes with zero counts in all samples
    counts = counts.loc[:, (counts.sum(axis=0) > 0)]

    design = config.design.strip()
    if not design.startswith("~"):
        design = f"~{design}"

    logger.info(f"Running PyDESeq2 with design {design}")
    dds = DeseqDataSet(counts=counts, metadata=metadata, design=design)
    dds.deseq2()
    write_gene_norm_counts(config, dds)

    # Contrast: last level vs first (alphabetical) — standard for 2 groups
    conditions = sorted(metadata["condition"].unique())
    if len(conditions) >= 2:
        contrast = [conditions[-1], conditions[0]]
    else:
        contrast = None

    stat_res = DeseqStats(dds, contrast=contrast)
    stat_res.summary()
    results = stat_res.results_df.copy()
    results.index.name = "protein_id"
    results = results.reset_index()

    results["is_deg"] = (
        (results["padj"].fillna(1.0) < config.alpha)
        & (results["log2FoldChange"].abs() >= config.lfc_threshold)
    )

    out_cols = [
        "protein_id", "baseMean", "log2FoldChange", "pvalue", "padj", "is_deg"
    ]
    for c in out_cols:
        if c not in results.columns and c != "protein_id":
            results[c] = float("nan")

    gene_deg = results[out_cols]
    gene_deg_path = Path(config.output_dir) / GENE_DEG_FILE
    gene_deg.to_csv(gene_deg_path, sep="\t", index=False)
    logger.info(f"Wrote {gene_deg_path}")

    deg_only = gene_deg[gene_deg["is_deg"]][["protein_id", "log2FoldChange"]]
    deg_path = Path(config.output_dir) / DEG_FILE
    deg_only.to_csv(deg_path, sep="\t", index=False, header=False)
    logger.info(f"Wrote {deg_path} ({len(deg_only)} significant genes)")

    return gene_deg


def run_cgc_pseudobulk_deseq2(
    config: ExpressionConfig,
    samplesheet: Samplesheet,
    cgc_count_matrix: pd.DataFrame,
) -> pd.DataFrame:
    """Run DESeq2 on CGC pseudobulk summed counts."""
    DeseqDataSet, DeseqStats = _require_pydeseq2()
    validate_deseq2_design(samplesheet, config.design)

    metadata = pd.DataFrame(
        [{"condition": s.condition} for s in samplesheet],
        index=[s.sample_id for s in samplesheet],
    )
    counts = cgc_count_matrix.loc[metadata.index]
    counts = counts.loc[:, (counts.sum(axis=0) > 0)]

    design = config.design.strip()
    if not design.startswith("~"):
        design = f"~{design}"

    dds = DeseqDataSet(counts=counts, metadata=metadata, design=design)
    dds.deseq2()
    conditions = sorted(metadata["condition"].unique())
    contrast = [conditions[-1], conditions[0]] if len(conditions) >= 2 else None
    stat_res = DeseqStats(dds, contrast=contrast)
    stat_res.summary()
    results = stat_res.results_df.copy()
    results.index.name = "cgcid"
    results = results.reset_index()
    results.rename(
        columns={
            "log2FoldChange": "cgc_log2FoldChange",
            "padj": "cgc_padj",
            "pvalue": "cgc_pvalue",
            "baseMean": "cgc_baseMean",
        },
        inplace=True,
    )
    return results
