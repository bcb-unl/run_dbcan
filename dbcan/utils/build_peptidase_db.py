"""Rebuild peptidase_db.dmnd (the MEROPS peptidase-unit DIAMOND database) from MEROPS pepunit.lib.

Maintainer tool; run_dbcan itself only downloads the finished database::

    python -m dbcan.utils.build_peptidase_db --input pepunit.lib --output peptidase_db.dmnd

Steps:

1. Parse each pepunit.lib header. The MEROPS identifier is the bracketed token right before
   ``#<family>#{<type> unit: ...}`` (the *last* one: names and organisms can contain brackets,
   e.g. ``[misleading]`` or ``([Clostridium] irregulare)``).
2. Leave out inhibitor units (``{inhibitor unit: ...}`` or an ``I`` family), entries whose
   identifier does not match the Peptidase signature pattern
   (``process_utils_constants.MEROPS_PEPTIDASE_ID_PATTERN``), entries without a valid
   ``MER`` accession or without sequence, and exact duplicates (same MEROPS identifier and
   same sequence; the first accession is kept). ``--keep-inhibitors`` keeps inhibitor units with a
   valid ``I`` identifier as decoys: run_dbcan drops hits to them, but with --max-target-seqs 1
   they stop a protein whose best match is an inhibitor from taking a weaker peptidase hit.
3. Write ``>MER0000001|S01.152`` records (PeptidaseDiamondProcessor reports the part after
   ``|`` as ``Peptidase ID``) and run ``diamond makedb`` (``--no-makedb`` skips it).
"""
import gzip
import hashlib
import logging
import re
import shutil
import subprocess
import sys
from collections import Counter
from contextlib import ExitStack
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Tuple

import rich_click as click

import dbcan.constants.peptidase_db_constants as P
from dbcan.configs.peptidase_db_config import PeptidaseDbBuildConfig
from dbcan.parameter import build_peptidase_db_options, create_config

logger = logging.getLogger(__name__)

_ACCESSION_RE = re.compile(P.MEROPS_ACCESSION_PATTERN)
_ID_FIELD_RE = re.compile(P.MEROPS_ID_FIELD_PATTERN)
_BRACKET_RE = re.compile(P.MEROPS_BRACKET_PATTERN)
_PEPTIDASE_ID_RE = re.compile(P.MEROPS_PEPTIDASE_ID_PATTERN)
_INHIBITOR_ID_RE = re.compile(P.MEROPS_INHIBITOR_ID_PATTERN)

# Number of left-out examples per reason written to the log.
_LOG_EXAMPLES = 5


class PeptidaseDbBuildError(RuntimeError):
    """The peptidase database could not be built (no input, no entries left, makedb failed)."""


@dataclass
class MeropsEntry:
    """One pepunit.lib header, parsed."""
    accession: str
    merops_id: Optional[str]
    unit_type: Optional[str]          # 'peptidase' / 'inhibitor' / None when the tail is missing
    header: str

    @property
    def is_inhibitor(self) -> bool:
        return (self.unit_type or "").lower() == P.UNIT_TYPE_INHIBITOR or bool(
            self.merops_id and _INHIBITOR_ID_RE.match(self.merops_id)
        )


@dataclass
class BuildSummary:
    total: int = 0
    kept: int = 0
    inhibitors_kept: int = 0
    dropped: Counter = field(default_factory=Counter)
    examples: Dict[str, List[str]] = field(default_factory=dict)
    families: Counter = field(default_factory=Counter)

    def drop(self, reason: str, entry: MeropsEntry) -> None:
        self.dropped[reason] += 1
        ex = self.examples.setdefault(reason, [])
        if len(ex) < _LOG_EXAMPLES:
            ex.append(entry.header)


def parse_header(header: str) -> MeropsEntry:
    """Parse a pepunit.lib header (with or without the leading '>')."""
    header = header.strip()
    if header.startswith(">"):
        header = header[1:]
    accession = header.split(None, 1)[0] if header else ""
    matches = list(_ID_FIELD_RE.finditer(header))
    if matches:
        last = matches[-1]
        merops_id, unit_type = last.group(1).strip(), last.group(3).strip().lower()
    else:
        brackets = _BRACKET_RE.findall(header)
        merops_id = brackets[-1].strip() if brackets else None
        unit_type = None
    return MeropsEntry(accession=accession, merops_id=merops_id, unit_type=unit_type, header=header)


def _open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def iter_fasta(path: Path) -> Iterator[Tuple[str, str]]:
    """Yield (header, sequence) from a FASTA file (plain or .gz)."""
    header: Optional[str] = None
    chunks: List[str] = []
    with _open_text(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header, chunks = line[1:].rstrip("\r\n"), []
            elif header is not None:
                chunks.append("".join(line.split()))
    if header is not None:
        yield header, "".join(chunks)


def _clean_sequence(seq: str) -> str:
    return seq.upper().rstrip("*")


def _write_record(out, record_id: str, seq: str) -> None:
    out.write(f">{record_id}\n")
    for i in range(0, len(seq), P.FASTA_LINE_WIDTH):
        out.write(seq[i:i + P.FASTA_LINE_WIDTH] + "\n")


def write_clean_fasta(input_path: Path, fasta_path: Path, report_path: Optional[Path] = None,
                      keep_inhibitors: bool = False) -> BuildSummary:
    """Filter pepunit.lib into a FASTA with '<MER accession>|<MEROPS id>' records.

    keep_inhibitors: keep inhibitor units with a valid ``I`` identifier as decoys (the
    Peptidase signature step drops hits to them, see process_utils._drop_non_peptidase_units).
    """
    summary = BuildSummary()
    seen_pairs = set()
    seen_accessions = set()
    fasta_path.parent.mkdir(parents=True, exist_ok=True)
    with ExitStack() as stack:
        out = stack.enter_context(open(fasta_path, "w", encoding="utf-8"))
        report = None
        if report_path is not None:
            report_path.parent.mkdir(parents=True, exist_ok=True)
            report = stack.enter_context(open(report_path, "w", encoding="utf-8"))
            report.write("\t".join(P.REPORT_COLUMNS) + "\n")

        def _drop(reason: str, entry: MeropsEntry) -> None:
            summary.drop(reason, entry)
            if report is not None:
                row = (entry.accession, entry.merops_id or "", reason, entry.header.replace("\t", " "))
                report.write("\t".join(row) + "\n")

        for header, raw_seq in iter_fasta(input_path):
            summary.total += 1
            entry = parse_header(header)
            seq = _clean_sequence(raw_seq)
            inhibitor = entry.is_inhibitor
            valid_id_re = _INHIBITOR_ID_RE if inhibitor else _PEPTIDASE_ID_RE
            if not _ACCESSION_RE.match(entry.accession):
                _drop(P.DROP_INVALID_ACCESSION, entry)
            elif inhibitor and not keep_inhibitors:
                _drop(P.DROP_INHIBITOR, entry)
            elif not entry.merops_id or not valid_id_re.match(entry.merops_id):
                _drop(P.DROP_INVALID_ID, entry)
            elif not seq:
                _drop(P.DROP_EMPTY_SEQUENCE, entry)
            else:
                key = hashlib.blake2b(f"{entry.merops_id}\0{seq}".encode(), digest_size=16).digest()
                if key in seen_pairs:
                    _drop(P.DROP_DUPLICATE, entry)
                elif entry.accession in seen_accessions:
                    _drop(P.DROP_DUPLICATE_ACCESSION, entry)
                else:
                    seen_pairs.add(key)
                    seen_accessions.add(entry.accession)
                    _write_record(out, f"{entry.accession}{P.OUTPUT_ID_SEPARATOR}{entry.merops_id}", seq)
                    summary.kept += 1
                    summary.inhibitors_kept += int(inhibitor)
                    summary.families[entry.merops_id.split(".", 1)[0]] += 1
    return summary


def run_makedb(diamond_bin: str, fasta_path: Path, db_path: Path, threads: int) -> None:
    exe = shutil.which(diamond_bin)
    if exe is None:
        raise PeptidaseDbBuildError(
            f"DIAMOND executable not found: {diamond_bin!r} (install DIAMOND or use --no-makedb)"
        )
    db_path.parent.mkdir(parents=True, exist_ok=True)
    # diamond makedb appends '.dmnd' when it is missing; pass the name without it.
    db_arg = str(db_path.with_suffix("")) if db_path.suffix == ".dmnd" else str(db_path)
    cmd = [exe, P.DIAMOND_MAKEDB, "--in", str(fasta_path), "--db", db_arg, "--threads", str(threads)]
    logger.info("Running: %s", " ".join(cmd))
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        raise PeptidaseDbBuildError(f"diamond makedb failed (exit {res.returncode}):\n{res.stderr.strip()}")


def default_fasta_path(output: Path) -> Path:
    return output.with_suffix(P.FASTA_SUFFIX) if output.suffix else Path(f"{output}{P.FASTA_SUFFIX}")


def build_peptidase_db(config: PeptidaseDbBuildConfig) -> BuildSummary:
    """Write the cleaned FASTA and (unless config.makedb is False) the DIAMOND database."""
    input_path = Path(config.input)
    if not input_path.is_file():
        raise PeptidaseDbBuildError(f"MEROPS pepunit.lib not found: {input_path}")
    output = Path(config.output)
    fasta_path = Path(config.fasta) if config.fasta else default_fasta_path(output)
    if fasta_path.resolve() == input_path.resolve():
        raise PeptidaseDbBuildError(f"--fasta must differ from --input: {fasta_path}")
    report_path = Path(config.report) if config.report else None

    summary = write_clean_fasta(input_path, fasta_path, report_path, keep_inhibitors=config.keep_inhibitors)
    _log_summary(summary, fasta_path)
    if summary.kept - summary.inhibitors_kept == 0:
        raise PeptidaseDbBuildError(
            f"No MEROPS peptidase units left in {input_path} ({summary.total} records read); "
            "is this pepunit.lib?"
        )
    if config.makedb:
        run_makedb(config.diamond_bin, fasta_path, output, config.threads)
        if not output.is_file():
            raise PeptidaseDbBuildError(f"diamond makedb finished but {output} was not written")
        logger.info("Wrote %s", output)
    return summary


def _log_summary(summary: BuildSummary, fasta_path: Path) -> None:
    logger.info(
        "pepunit.lib: %d records, %d units kept (%d families; %d inhibitor units kept as decoys) -> %s",
        summary.total, summary.kept, len(summary.families), summary.inhibitors_kept, fasta_path,
    )
    for reason in P.DROP_REASONS:
        n = summary.dropped.get(reason, 0)
        if n:
            logger.info("  left out %-26s %d", reason, n)
            for ex in summary.examples.get(reason, []):
                logger.debug("    e.g. %s", ex)


@click.command(help=(
    "Rebuild peptidase_db.dmnd from MEROPS pepunit.lib: keep peptidase units with a valid MEROPS "
    "identifier, drop inhibitor units, malformed entries and exact duplicates, then run "
    "'diamond makedb'. See the 'Prepare the database' page of the documentation."
))
@build_peptidase_db_options
def cli(**kwargs):
    log_level = kwargs.pop("log_level")
    logging.basicConfig(level=getattr(logging, log_level.upper()),
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", force=True)
    config = create_config(PeptidaseDbBuildConfig, **kwargs)
    try:
        build_peptidase_db(config)
    except PeptidaseDbBuildError as e:
        logger.error("%s", e)
        sys.exit(1)


if __name__ == "__main__":
    cli()
