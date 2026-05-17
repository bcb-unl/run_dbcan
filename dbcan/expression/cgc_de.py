"""CGC-level differential expression summaries."""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

from dbcan.configs.expression_config import ExpressionConfig
from dbcan.constants.expression_constants import CGC_DEG_FILE, CGC_GENE_DEG_MAP_FILE
from dbcan.expression.deseq2_runner import run_cgc_pseudobulk_deseq2
from dbcan.expression.samplesheet import Samplesheet
from dbcan.utils.utils import Read_cgc_standard_out

logger = logging.getLogger(__name__)


def _cgc_de_rule_passes(n_deg: int, n_genes: int, rule: str) -> bool:
    if n_genes == 0:
        return False
    if rule == "any":
        return n_deg >= 1
    if rule == "majority":
        return n_deg > n_genes / 2
    if rule == "all":
        return n_deg == n_genes
    raise ValueError(f"Unknown cgc_de_rule: {rule}")


def build_cgc_pseudobulk_matrix(
    count_matrix: pd.DataFrame,
    cgcid2records: dict,
) -> pd.DataFrame:
    """Sum gene counts per CGC per sample."""
    rows = []
    for sample_id in count_matrix.index:
        row = {}
        for cgcid, records in cgcid2records.items():
            total = 0
            for rec in records:
                if rec.seqid in count_matrix.columns:
                    total += int(count_matrix.loc[sample_id, rec.seqid])
            row[cgcid] = total
        rows.append(row)
    return pd.DataFrame(rows, index=count_matrix.index).fillna(0).astype(int)


def summarize_cgc_de(
    config: ExpressionConfig,
    samplesheet: Samplesheet,
    count_matrix: pd.DataFrame,
    gene_deg: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build cgc_gene_deg_map.tsv and cgc_deg.tsv."""
    cgc_path = Path(config.input_dir) / "cgc_standard_out.tsv"
    _, cgcid2records = Read_cgc_standard_out(str(cgc_path))

    deg_lookup = gene_deg.set_index("protein_id") if "protein_id" in gene_deg.columns else gene_deg

    map_rows = []
    for cgcid, records in cgcid2records.items():
        for rec in records:
            gid = rec.seqid
            if gid in deg_lookup.index:
                g = deg_lookup.loc[gid]
                map_rows.append({
                    "cgcid": cgcid,
                    "protein_id": gid,
                    "gene_type": rec.gene_type,
                    "log2FoldChange": g.get("log2FoldChange", np.nan),
                    "padj": g.get("padj", np.nan),
                    "is_deg": bool(g.get("is_deg", False)),
                })
            else:
                map_rows.append({
                    "cgcid": cgcid,
                    "protein_id": gid,
                    "gene_type": rec.gene_type,
                    "log2FoldChange": np.nan,
                    "padj": np.nan,
                    "is_deg": False,
                })

    map_df = pd.DataFrame(map_rows)
    map_path = Path(config.output_dir) / CGC_GENE_DEG_MAP_FILE
    map_df.to_csv(map_path, sep="\t", index=False)
    logger.info(f"Wrote {map_path}")

    # Pseudobulk DESeq2
    cgc_counts = build_cgc_pseudobulk_matrix(count_matrix, cgcid2records)
    try:
        pseudo = run_cgc_pseudobulk_deseq2(config, samplesheet, cgc_counts)
    except ImportError:
        pseudo = pd.DataFrame({"cgcid": list(cgcid2records.keys())})
        for col in ("cgc_log2FoldChange", "cgc_padj", "cgc_pvalue", "cgc_baseMean"):
            pseudo[col] = np.nan

    cgc_rows = []
    for cgcid in cgcid2records:
        sub = map_df[map_df["cgcid"] == cgcid]
        n_genes = len(sub)
        n_deg = int(sub["is_deg"].sum())
        rule_de = _cgc_de_rule_passes(n_deg, n_genes, config.cgc_de_rule)
        pseudo_row = pseudo[pseudo["cgcid"] == cgcid] if "cgcid" in pseudo.columns else pd.DataFrame()
        cgc_padj = pseudo_row["cgc_padj"].iloc[0] if len(pseudo_row) else np.nan
        cgc_lfc = pseudo_row["cgc_log2FoldChange"].iloc[0] if len(pseudo_row) else np.nan
        pseudo_de = bool(
            pd.notna(cgc_padj)
            and cgc_padj < config.alpha
            and abs(cgc_lfc) >= config.lfc_threshold
        )
        cgc_rows.append({
            "cgcid": cgcid,
            "n_genes": n_genes,
            "n_deg_genes": n_deg,
            "is_cgc_de_rule": rule_de,
            "is_cgc_de_pseudobulk": pseudo_de,
            "is_cgc_de": rule_de or pseudo_de,
            "cgc_log2FoldChange": cgc_lfc,
            "cgc_padj": cgc_padj,
        })

    cgc_deg = pd.DataFrame(cgc_rows)
    if not pseudo.empty and "cgcid" in pseudo.columns:
        extra = [c for c in pseudo.columns if c not in cgc_deg.columns]
        if extra:
            cgc_deg = cgc_deg.merge(pseudo[["cgcid"] + extra], on="cgcid", how="left")

    deg_path = Path(config.output_dir) / CGC_DEG_FILE
    cgc_deg.to_csv(deg_path, sep="\t", index=False)
    logger.info(f"Wrote {deg_path}")
    return map_df, cgc_deg
