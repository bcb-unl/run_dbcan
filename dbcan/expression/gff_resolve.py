"""Resolve GFF path for expression analysis."""

from __future__ import annotations

import logging
from pathlib import Path

from dbcan.configs.expression_config import ExpressionConfig
from dbcan.constants.base_constants import CGC_GFF_FILE

logger = logging.getLogger(__name__)


def resolve_expression_gff(config: ExpressionConfig) -> str:
    """
    Choose GFF for cal_coverage.

    Priority:
    1. Explicit --gff if provided and exists
    2. {input_dir}/cgc.gff (run_dbcan output; protein_id= matches cgc_standard_out)
    3. {input_dir}/*.fix.gff (first match)
    4. {input_dir}/uniInput.gff
    """
    if config.gff:
        p = Path(config.gff)
        if p.is_file():
            logger.info(f"Using GFF: {p}")
            return str(p.resolve())
        raise FileNotFoundError(f"GFF not found: {config.gff}")

    dbcan_dir = Path(config.input_dir)
    cgc_gff = dbcan_dir / CGC_GFF_FILE
    if cgc_gff.is_file():
        logger.info(f"Using run_dbcan CGC GFF: {cgc_gff}")
        return str(cgc_gff.resolve())

    fix_candidates = sorted(dbcan_dir.glob("*.fix.gff"))
    if fix_candidates:
        logger.info(f"Using fix GFF: {fix_candidates[0]}")
        return str(fix_candidates[0].resolve())

    uni = dbcan_dir / "uniInput.gff"
    if uni.is_file():
        logger.info(f"Using uniInput GFF: {uni}")
        return str(uni.resolve())

    raise FileNotFoundError(
        f"No GFF found under {dbcan_dir}. Pass --gff or ensure {CGC_GFF_FILE} exists."
    )
