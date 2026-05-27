"""BWA alignment and BAM indexing for expression analysis."""

from __future__ import annotations

import logging
import subprocess
from pathlib import Path
from typing import Optional

from dbcan.configs.expression_config import ExpressionConfig
from dbcan.expression.samplesheet import SampleRecord, Samplesheet

logger = logging.getLogger(__name__)


def _run_cmd(cmd: list, desc: str) -> None:
    logger.info(f"{desc}: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)


def ensure_bwa_index(reference_fasta: str, threads: int) -> None:
    """Build BWA index if .bwt is missing."""
    ref = Path(reference_fasta)
    if not ref.is_file():
        raise FileNotFoundError(f"Reference FASTA not found: {reference_fasta}")
    bwt = ref.parent / f"{ref.name}.bwt"
    alt_bwt = Path(f"{reference_fasta}.bwt")
    if bwt.exists() or alt_bwt.exists():
        return
    _run_cmd(["bwa", "index", "-p", str(ref), str(ref)], "bwa index")


def align_sample(
    sample: SampleRecord,
    reference_fasta: str,
    output_bam: Path,
    threads: int,
) -> Path:
    """Align FASTQ to reference with bwa mem and sort/index with samtools."""
    ref = Path(reference_fasta)
    sam_path = output_bam.with_suffix(".sam")
    ensure_bwa_index(reference_fasta, threads)

    cmd = ["bwa", "mem", "-t", str(threads), str(ref), sample.r1]
    if sample.r2:
        cmd.append(sample.r2)
    logger.info(f"Aligning sample {sample.sample_id}")
    with sam_path.open("w") as sam_out:
        subprocess.run(cmd, stdout=sam_out, check=True)

    _run_cmd(
        ["samtools", "sort", "-@", str(threads), "-o", str(output_bam), str(sam_path)],
        f"samtools sort {sample.sample_id}",
    )
    _run_cmd(["samtools", "index", str(output_bam)], f"samtools index {sample.sample_id}")
    if sam_path.exists():
        sam_path.unlink()
    return output_bam


def resolve_bam_paths(
    config: ExpressionConfig,
    samplesheet: Samplesheet,
) -> dict[str, Path]:
    """Return sample_id -> BAM path, aligning from FASTQ when needed."""
    out_dir = Path(config.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    bam_map: dict[str, Path] = {}

    for sample in samplesheet:
        target = out_dir / f"{sample.sample_id}.bam"
        if sample.has_bam and (config.skip_alignment or not sample.has_fastq):
            src = Path(sample.bam)
            if not src.is_file():
                raise FileNotFoundError(f"BAM not found for {sample.sample_id}: {src}")
            if src.resolve() != target.resolve():
                if target.exists() or target.is_symlink():
                    target.unlink()
                target.symlink_to(src.resolve())
            bam_map[sample.sample_id] = target
            if not Path(f"{target}.bai").is_file() and not Path(f"{target}.bam.bai").is_file():
                _run_cmd(["samtools", "index", str(target)], f"samtools index {sample.sample_id}")
        elif sample.has_fastq:
            if not config.reference_fasta:
                raise ValueError(f"reference_fasta required to align sample {sample.sample_id}")
            bam_map[sample.sample_id] = align_sample(
                sample, config.reference_fasta, target, config.threads
            )
        else:
            raise ValueError(f"No alignment input for sample {sample.sample_id}")
    return bam_map
