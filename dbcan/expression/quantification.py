"""Gene-level read counting and count matrix construction."""

from __future__ import annotations

import logging
from pathlib import Path
from types import SimpleNamespace

import pandas as pd

from dbcan.configs.expression_config import ExpressionConfig
from dbcan.constants.base_constants import CGC_RESULT_FILE, OVERVIEW_FILE
from dbcan.constants.expression_constants import COUNT_MATRIX_FILE
from dbcan.expression.samplesheet import Samplesheet
from dbcan.utils.utils import ReadBedtools, cal_coverage

logger = logging.getLogger(__name__)


def _gene_filter_set(config: ExpressionConfig) -> set[str] | None:
    """Genes to include in count matrix; None means all GFF genes."""
    if config.all_genes:
        return None
    genes: set[str] = set()
    dbcan_dir = Path(config.input_dir)
    overview = dbcan_dir / OVERVIEW_FILE
    cgc = dbcan_dir / CGC_RESULT_FILE
    if overview.is_file():
        with overview.open() as fh:
            header = fh.readline().rstrip("\n").split("\t")
            gene_col = next((h for h in header if h.startswith("Gene ID")), header[0])
            idx = header.index(gene_col)
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if parts:
                    genes.add(parts[idx])
    if cgc.is_file():
        with cgc.open() as fh:
            next(fh, None)
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 4:
                    genes.add(parts[3])
    return genes if genes else None


def run_cal_coverage(
    bam_path: Path,
    gff_path: str,
    depth_out: Path,
    config: ExpressionConfig,
) -> Path:
    """Run cal_coverage for one sample."""
    args = SimpleNamespace(
        input=str(bam_path),
        gff=gff_path,
        output=str(depth_out),
        overlap_base_ratio=config.overlap_base_ratio,
        mapping_quality=config.mapping_quality,
        identity=config.identity,
        threads=config.threads,
        hifi=config.hifi,
    )
    depth_out.parent.mkdir(parents=True, exist_ok=True)
    cal_coverage(args)
    return depth_out


def build_count_matrix(
    config: ExpressionConfig,
    samplesheet: Samplesheet,
    bam_map: dict[str, Path],
) -> tuple[pd.DataFrame, dict[str, Path]]:
    """Count reads per gene for each sample and merge into count matrix."""
    gene_filter = _gene_filter_set(config)
    rows = []
    depth_paths: dict[str, Path] = {}

    for sample in samplesheet:
        sample_dir = Path(config.output_dir) / sample.sample_id
        depth_path = sample_dir / f"{sample.sample_id}.depth.txt"
        run_cal_coverage(bam_map[sample.sample_id], config.gff, depth_path, config)
        depth_paths[sample.sample_id] = depth_path

        seqid2info, _ = ReadBedtools(str(depth_path))
        row = {"sample_id": sample.sample_id}
        for gid, rec in seqid2info.items():
            if gene_filter is not None and gid not in gene_filter:
                continue
            row[gid] = int(rec.read_count)
        rows.append(row)

    df = pd.DataFrame(rows).set_index("sample_id").fillna(0).astype(int)
    out_path = Path(config.output_dir) / COUNT_MATRIX_FILE
    df.to_csv(out_path, sep="\t")
    logger.info(f"Wrote count matrix: {out_path} ({df.shape[0]} samples x {df.shape[1]} genes)")
    return df, depth_paths
