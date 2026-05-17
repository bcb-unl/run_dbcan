"""Per-sample CAZyme/CGC abundance estimation."""

from __future__ import annotations

import logging
import os
from pathlib import Path

from dbcan.configs.expression_config import ExpressionConfig
from dbcan.configs.utils_config import AbundanceConfig
from dbcan.constants.expression_constants import (
    CGC_EXPRESSION_MATRIX_FILE,
    CGC_GENE_COUNTS_FILE,
)
from dbcan.expression.samplesheet import Samplesheet
from dbcan.utils.utils import (
    CAZyme_Abundance_estimate,
    AbundParameters,
    ReadBedtools,
    Read_cgc_standard_out,
)

logger = logging.getLogger(__name__)

ABUNDANCE_FUNCTIONS = (
    "fam_abund",
    "fam_substrate_abund",
    "CGC_abund",
    "CGC_substrate_abund",
)


def _run_abundance_function(
    function: str,
    config: ExpressionConfig,
    depth_path: Path,
    sample_dir: Path,
) -> None:
    """Run one dbcan_utils abundance function in sample_dir."""
    cwd = os.getcwd()
    input_dir = str(Path(config.input_dir).resolve())
    depth_abs = str(Path(depth_path).resolve())
    sample_dir = Path(sample_dir).resolve()
    try:
        os.chdir(sample_dir)
        abund_cfg = AbundanceConfig(
            input_dir=input_dir,
            bedtools_depth=depth_abs,
            output_dir=str(sample_dir),
            abundance=config.abundance,
            overlap_base_ratio=config.overlap_base_ratio,
            mapping_quality=config.mapping_quality,
            identity=config.identity,
            threads=config.threads,
            hifi=config.hifi,
        )
        pars = AbundParameters(abund_cfg, function)
        est = CAZyme_Abundance_estimate(pars)
        est.Cal_Seq_Abundance(config.abundance)
        if function == "fam_abund":
            est.Cal_Famliy_Abundance()
            est.output_family_abund()
        elif function == "fam_substrate_abund":
            est.Cal_Substrate_Abundance()
            est.output_substrate_abund()
        elif function == "CGC_abund":
            est.Cal_PUL_Abundance()
            est.output_cgc_abund()
        elif function == "CGC_substrate_abund":
            est.Cal_PUL_Abundance()
            est.Cal_PUL_Substrate_Abundance()
            est.output_cgcsubstrate_abund()
    finally:
        os.chdir(cwd)


def run_per_sample_abundance(
    config: ExpressionConfig,
    samplesheet: Samplesheet,
    depth_paths: dict[str, Path],
) -> None:
    """Compute fam/CGC abundance for each sample."""
    for sample in samplesheet:
        sample_dir = Path(config.output_dir) / sample.sample_id
        sample_dir.mkdir(parents=True, exist_ok=True)
        depth = depth_paths[sample.sample_id]
        for func in ABUNDANCE_FUNCTIONS:
            try:
                _run_abundance_function(func, config, depth, sample_dir)
            except FileNotFoundError as e:
                if func in ("fam_substrate_abund", "CGC_substrate_abund"):
                    logger.warning(f"Skipping {func} for {sample.sample_id}: {e}")
                else:
                    raise


def build_cgc_gene_counts_and_matrix(
    config: ExpressionConfig,
    samplesheet: Samplesheet,
    depth_paths: dict[str, Path],
) -> None:
    """Write long-format cgc_gene_counts.tsv and wide cgc_expression_matrix.tsv."""
    cgc_path = Path(config.input_dir).resolve() / "cgc_standard_out.tsv"
    if not cgc_path.is_file():
        logger.warning("cgc_standard_out.tsv not found; skipping CGC count tables")
        return

    _, cgcid2records = Read_cgc_standard_out(str(cgc_path))
    long_rows = []
    cgc_sample_tpm: dict[str, dict[str, float]] = {}

    for sample in samplesheet:
        depth = depth_paths[sample.sample_id]
        seqid2info, normalized_tpm = ReadBedtools(str(depth))
        norm_reads = max(sum(r.read_count for r in seqid2info.values()), 1.0) / 1e6
        denom_tpm = max(normalized_tpm, 1e-9)

        for cgcid, records in cgcid2records.items():
            tpms = []
            for rec in records:
                info = seqid2info.get(rec.seqid)
                if not info:
                    count, length = 0, max(rec.gene_end - rec.gene_start + 1, 1)
                else:
                    count = int(info.read_count)
                    length = max(int(info.length), 1)
                if config.abundance == "TPM":
                    val = (count / length * 1e6) / denom_tpm
                elif config.abundance == "FPKM":
                    val = count / norm_reads / (length / 1000.0)
                else:
                    val = count / norm_reads
                tpms.append(val)
                long_rows.append({
                    "cgcid": cgcid,
                    "protein_id": rec.seqid,
                    "sample_id": sample.sample_id,
                    "count": count,
                    "tpm": round(val, 4),
                })
            cgc_sample_tpm.setdefault(cgcid, {})[sample.sample_id] = round(
                sum(tpms) / len(tpms) if tpms else 0.0, 4
            )

    import pandas as pd

    long_df = pd.DataFrame(long_rows)
    long_path = Path(config.output_dir) / CGC_GENE_COUNTS_FILE
    long_df.to_csv(long_path, sep="\t", index=False)
    logger.info(f"Wrote {long_path}")

    if cgc_sample_tpm:
        wide = pd.DataFrame(cgc_sample_tpm).T
        wide.index.name = "cgcid"
        wide_path = Path(config.output_dir) / CGC_EXPRESSION_MATRIX_FILE
        wide.to_csv(wide_path, sep="\t")
        logger.info(f"Wrote {wide_path}")
