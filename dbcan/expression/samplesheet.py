"""Parse and validate expression analysis samplesheets."""

from __future__ import annotations

import csv
import logging
import re
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

logger = logging.getLogger(__name__)

REQUIRED_COLUMNS = {"sample_id", "condition"}
OPTIONAL_COLUMNS = {"bam", "r1", "r2", "replicate_group"}


@dataclass
class SampleRecord:
    sample_id: str
    condition: str
    bam: Optional[str] = None
    r1: Optional[str] = None
    r2: Optional[str] = None
    replicate_group: Optional[str] = None

    @property
    def has_bam(self) -> bool:
        return bool(self.bam and str(self.bam).strip())

    @property
    def has_fastq(self) -> bool:
        return bool(self.r1 and str(self.r1).strip())


@dataclass
class Samplesheet:
    samples: List[SampleRecord]
    path: Path

    def __iter__(self):
        return iter(self.samples)

    def sample_ids(self) -> List[str]:
        return [s.sample_id for s in self.samples]


def _normalize_header(name: str) -> str:
    """Normalize column names; collapse whitespace so 'condition   ' -> 'condition'."""
    s = name.strip().lower()
    s = re.sub(r"\s+", "_", s)
    s = re.sub(r"_+", "_", s).strip("_")
    return s


# Accepted aliases for required logical columns
_SAMPLE_ID_ALIASES = ("sample_id", "sample", "sampleid", "id")
_CONDITION_ALIASES = ("condition", "group", "treatment", "factor", "phenotype")


def _resolve_column(colmap: dict[str, str], aliases: tuple[str, ...], logical: str) -> str:
    for alias in aliases:
        if alias in colmap:
            return colmap[alias]
    raise KeyError(logical)


def load_samplesheet(path: str) -> Samplesheet:
    """Load TSV samplesheet with required sample_id and condition columns."""
    p = Path(path)
    if not p.is_file():
        raise FileNotFoundError(f"Samplesheet not found: {path}")

    with p.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Samplesheet has no header: {path}")
        colmap = {_normalize_header(h): h for h in reader.fieldnames if h}
        try:
            sample_id_key = _resolve_column(colmap, _SAMPLE_ID_ALIASES, "sample_id")
            condition_key = _resolve_column(colmap, _CONDITION_ALIASES, "condition")
        except KeyError as e:
            hint = ""
            if any("r1" in k or "r2" in k for k in colmap):
                merged = [k for k in colmap if "condition" in k and ("r1" in k or "r2" in k)]
                if merged:
                    hint = (
                        f" Header may use spaces instead of tabs between columns "
                        f"(found merged column {merged!r}). Use TAB between sample_id, condition, r1, r2."
                    )
            raise ValueError(
                f"Samplesheet missing required column '{e.args[0]}' "
                f"(expected one of {_SAMPLE_ID_ALIASES if e.args[0] == 'sample_id' else _CONDITION_ALIASES}): {path}.{hint}"
            ) from e

        samples: List[SampleRecord] = []
        seen_ids = set()
        for row_num, row in enumerate(reader, start=2):
            sid = (row.get(sample_id_key) or "").strip()
            cond = (row.get(condition_key) or "").strip()
            if not sid:
                logger.warning(f"Skipping samplesheet row {row_num}: empty sample_id")
                continue
            if sid in seen_ids:
                raise ValueError(f"Duplicate sample_id '{sid}' at row {row_num}")
            seen_ids.add(sid)

            def _get(col: str) -> Optional[str]:
                key = colmap.get(col)
                if not key:
                    return None
                val = (row.get(key) or "").strip()
                return val or None

            rec = SampleRecord(
                sample_id=sid,
                condition=cond,
                bam=_get("bam"),
                r1=_get("r1"),
                r2=_get("r2"),
                replicate_group=_get("replicate_group"),
            )
            if not rec.has_bam and not rec.has_fastq:
                raise ValueError(
                    f"Sample '{sid}' (row {row_num}) must have bam or r1 fastq path"
                )
            if rec.has_bam and rec.has_fastq:
                logger.warning(
                    f"Sample '{sid}' has both bam and fastq; bam will be preferred when skip_alignment is set"
                )
            samples.append(rec)

    if len(samples) < 2:
        raise ValueError("Samplesheet must contain at least 2 samples for differential analysis")

    conditions = {s.condition for s in samples}
    if len(conditions) < 2:
        logger.warning(
            "Only one condition in samplesheet; DESeq2 contrast may be underpowered"
        )

    logger.info(f"Loaded {len(samples)} samples from {path}")
    return Samplesheet(samples=samples, path=p)
