"""dbcan_utils: read counting per gene (cal_coverage) and CAZyme / substrate / CGC abundance."""
import logging
import os
import re
from collections import Counter
from dataclasses import dataclass, field
from itertools import groupby
from multiprocessing import Pool
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pysam
import rich_click as click
from Bio import SeqIO
from tqdm import tqdm

import dbcan.constants.utils_constants as U
from dbcan.configs.utils_config import AbundanceConfig, CoverageConfig
from dbcan.constants.utils_constants import (
    CGC_RESULT_FILE,
    CGC_SUB_PREDICTION_FILE,
    DBCAN_SUB_RESULT_FILE,
    OVERVIEW_FILE,
)
from dbcan.core import DbcanStepError
from dbcan.IO.gff import _split_gff_attributes
from dbcan.main import setup_logging
from dbcan.parameter import dbcan_utils_options, logging_options

click.rich_click.USE_RICH_MARKUP = True
click.rich_click.USE_MARKDOWN = True
click.rich_click.SHOW_ARGUMENTS = True
click.rich_click.GROUP_ARGUMENTS_OPTIONS = True
click.rich_click.STYLE_ERRORS_SUGGESTION = "magenta italic"
click.rich_click.SHOW_METAVARS_COLUMN = False
click.rich_click.APPEND_METAVARS_HELP = True

logger = logging.getLogger(__name__)


class AbundanceError(DbcanStepError):
    """dbcan_utils cannot produce the requested output (the command exits non-zero)."""


def is_placeholder_id(gene_id: Optional[str]) -> bool:
    """run_dbcan writes 'NA' / 'unknown' as protein ID of non-coding or unresolved features."""
    return gene_id is None or gene_id.strip() in U.PLACEHOLDER_GENE_IDS


# ---------------------------------------------------------------------------
# Gene read-count table (cal_coverage output: id, length, count)
# ---------------------------------------------------------------------------
class bedtools_read_count:
    def __init__(self, lines: List[str]):
        self.seqid = lines[0]
        self.length = int(float(lines[1]))
        self.read_count = float(lines[2])

    def __repr__(self):
        return "\t".join([str(getattr(self, k)) for k in vars(self)])


def ReadBedtools(filename: str) -> Tuple[Dict[str, bedtools_read_count], float, float]:
    """Read the gene count table. Returns (id -> record, sum(count/length), sum(count))."""
    seqid2info: Dict[str, bedtools_read_count] = {}
    normalized_tpm = 0.0
    total = 0.0
    dup_ids: List[str] = []
    with open(filename) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                parts = line.split()
            if len(parts) < 3:
                continue
            try:
                rec = bedtools_read_count(parts[:3])
            except ValueError:
                continue
            if rec.seqid in seqid2info:
                dup_ids.append(rec.seqid)
                old = seqid2info[rec.seqid]
                total -= old.read_count
                if old.length > 0:
                    normalized_tpm -= old.read_count / old.length
            seqid2info[rec.seqid] = rec
            total += rec.read_count
            if rec.length > 0:
                normalized_tpm += rec.read_count / rec.length
    if dup_ids:
        logger.warning(f"{len(dup_ids)} rows of {filename} repeat a gene ID (the last row is kept), examples: "
                       f"{sorted(set(dup_ids))[:5]}")
    return seqid2info, normalized_tpm, total


def total_mapped_reads_count(file_name: str) -> float:
    """Total count assigned to genes in a cal_coverage table (the RPM/FPKM library size)."""
    if not os.path.exists(file_name):
        logger.warning(f"Read count file does not exist: {file_name}")
        return 0.0
    return ReadBedtools(file_name)[2]


# ---------------------------------------------------------------------------
# Annotation parsers
# ---------------------------------------------------------------------------
def Is_EC(ec: str) -> bool:
    return bool(ec) and ec[0].isdigit()


def parse_ec_field(ec_field: str) -> List[str]:
    """All EC numbers of an overview.tsv 'EC#' cell, in order of appearance, without duplicates.

    The cell is '<domain 1>|<domain 2>|...'; each domain is '-' or EC tokens joined by ';'
    and each token is 'a.b.c.d' optionally followed by ':<number of supporting proteins>'
    (e.g. '3.2.1.1:3;3.2.1.20:13|3.2.1.1;3.2.1.133').
    """
    out: List[str] = []
    for domain in (ec_field or "").split("|"):
        for token in domain.split(";"):
            ec = token.split(":")[0].strip()
            if ec and ec != "-" and Is_EC(ec) and ec not in out:
                out.append(ec)
    return out


# ';' separates the entries of different ECs, ',' several substrates of one entry. A comma
# between digits belongs to a chemical name (alpha-1,3-glucan) and is not a separator.
_SUBSTRATE_SPLIT = re.compile(r";|,(?!\d)")


def split_substrates(field_value: str) -> List[str]:
    """Unique substrates of a dbCAN-sub / CGC substrate cell ('-' and empty entries dropped)."""
    out: List[str] = []
    for piece in _SUBSTRATE_SPLIT.split(field_value or ""):
        piece = piece.strip()
        if piece and piece != "-" and piece not in out:
            out.append(piece)
    return out


def Clean_Hmmer_sub(pred: str) -> List[str]:
    if not pred or pred == "-":
        return []
    preds = pred.split("+")
    clean = [p.split("(")[0] for p in preds if p and p != "-"]
    return sorted(set(clean))


def Clean_diamond(pred: str) -> List[str]:
    if not pred or pred == "-":
        return []
    preds = pred.split("+")
    clean = [p.split("(")[0] for p in preds if p and not Is_EC(p)]
    return sorted(set(clean))


class OverView:
    def __init__(self, seqid: str, ec_field: str, hmmer_field: str, sub_field: str, dia_field: str):
        self.seqid = seqid
        self.ECs = parse_ec_field(ec_field)
        self.hmmer = Clean_Hmmer_sub(hmmer_field)
        self.dbcan_sub = Clean_Hmmer_sub(sub_field)
        self.diamond = Clean_diamond(dia_field)
        self.justify_final_pred()

    def __repr__(self):
        if not self.preds:
            return ""
        outline = []
        for value in vars(self):
            attr = getattr(self, value)
            outline.append("+".join(attr) if isinstance(attr, list) else str(attr))
        return "\t".join(outline)

    def justify_final_pred(self):
        union_fam = sorted(set(self.hmmer) | set(self.dbcan_sub)) + self.ECs
        self.preds = union_fam if union_fam else []


def ReadOverView(filename: str) -> Dict[str, OverView]:
    """Parse overview.tsv by header name."""
    with open(filename) as f:
        header = f.readline().rstrip("\n").split("\t")
        hmap = {h: i for i, h in enumerate(header)}
        gene_col = next((h for h in header if h.startswith("Gene ID")), None)
        ec_col = next((h for h in header if "EC#" in h), None)
        hmm_col = next((h for h in header if "dbCAN_hmm" in h), None)
        sub_col = next((h for h in header if "dbCAN_sub" in h or "dbcan_sub" in h), None)
        dia_col = next((h for h in header if "DIAMOND" in h), None)
        if not all([gene_col, ec_col, hmm_col, sub_col, dia_col]):
            logger.warning(f"overview.tsv missing required columns, header={header}")
        out: Dict[str, OverView] = {}
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            try:
                seqid = parts[hmap.get(gene_col, 0)]
                ec = parts[hmap[ec_col]] if ec_col else "-"
                hmm = parts[hmap[hmm_col]] if hmm_col else "-"
                sub = parts[hmap[sub_col]] if sub_col else "-"
                dia = parts[hmap[dia_col]] if dia_col else "-"
                out[seqid] = OverView(seqid, ec, hmm, sub, dia)
            except (IndexError, KeyError):
                continue
        return out


class DbcanSub_line:
    def __init__(self, lines: List[str]):
        self.dbcan_subfam = lines[0]
        self.subfam_comp = lines[1]
        self.subfam_EC = lines[2]
        self.substrate = lines[3]
        self.hmmlen = lines[4]
        self.seqid = lines[5]
        self.protlen = lines[6]
        self.evalue = lines[7]

    def __repr__(self):
        return "\t".join([str(getattr(self, k)) for k in vars(self)])


class DbcanSubProtein:
    """All dbCAN-sub hits of one protein (one row per domain in dbCANsub_hmm_results.tsv)."""

    def __init__(self, seqid: str):
        self.seqid = seqid
        self.hits: List[DbcanSub_line] = []
        self.substrates: List[str] = []

    def add(self, rec: DbcanSub_line) -> None:
        self.hits.append(rec)
        for sub in split_substrates(rec.substrate):
            if sub not in self.substrates:
                self.substrates.append(sub)


def Read_dbcansub_out(filename: str) -> Dict[str, DbcanSubProtein]:
    out: Dict[str, DbcanSubProtein] = {}
    with open(filename) as f:
        next(f, None)
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            rec = DbcanSub_line(parts[:8])
            out.setdefault(rec.seqid, DbcanSubProtein(rec.seqid)).add(rec)
    return out


class cgc_standard_line:
    def __init__(self, lines: List[str]):
        self.cgcid = lines[2] + "|" + lines[0]
        self.cgc_order = lines[0]
        self.gene_type = lines[1]
        self.contig_id = lines[2]
        self.seqid = lines[3]
        self.gene_start = int(lines[4])
        self.gene_end = int(lines[5])
        self.strand = lines[6]
        self.protfam = lines[7].split("|")[1] if lines[7] != "null" and "|" in lines[7] else (lines[7] if lines[7] != "null" else "null")

    def __repr__(self):
        return "\t".join([str(getattr(self, k)) for k in vars(self)])


def Read_cgc_standard_out(filename: str) -> Tuple[Dict[str, cgc_standard_line], Dict[str, List[cgc_standard_line]]]:
    """Return (protein ID -> record, CGC ID -> records in file order).

    Rows with a placeholder protein ID ('NA' / 'unknown': tRNA and other non-coding features)
    have no protein and no count in a cal_coverage table, so they are kept out of both maps;
    otherwise all of them would collapse into a single 'NA' entry.
    """
    seqid2records: Dict[str, cgc_standard_line] = {}
    cgcid2records: Dict[str, List[cgc_standard_line]] = {}
    skipped = 0
    with open(filename) as f:
        next(f, None)
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            rec = cgc_standard_line(parts)
            if is_placeholder_id(rec.seqid):
                skipped += 1
                continue
            seqid2records[rec.seqid] = rec
            cgcid2records.setdefault(rec.cgcid, []).append(rec)
    if skipped:
        logger.info(f"{skipped} non-coding CGC member(s) without a protein ID (e.g. tRNA) are not part of CGC abundance")
    return seqid2records, cgcid2records


class cgc_substrate:
    def __init__(self, lines: List[str]):
        self.cgcid = lines[0]
        self.homo_pul = lines[1] if len(lines) > 1 else ""
        self.homo_sub = lines[2] if len(lines) > 2 else ""
        self.bitscore = lines[3] if len(lines) > 3 else ""
        self.signature_pairs = lines[4] if len(lines) > 4 else ""
        self.major_voting_sub = lines[5] if len(lines) > 5 else ""
        self.major_voting_score = lines[6] if len(lines) > 6 else ""

    def __repr__(self):
        return "\t".join([str(getattr(self, k)) for k in vars(self)])


def Read_cgc_substrate(filename: str) -> Dict[str, cgc_substrate]:
    out: Dict[str, cgc_substrate] = {}
    with open(filename) as f:
        next(f, None)
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            rec = cgc_substrate(parts)
            out[rec.cgcid] = rec
    return out


# ---------------------------------------------------------------------------
# Abundance
# ---------------------------------------------------------------------------
class AbundParameters:
    """Derive and validate input/output paths for one *_abund function."""

    def __init__(self, config: AbundanceConfig, function: str):
        self.config = config
        self.function = function
        self.input_dir = config.input_dir
        self.bedtools = config.bedtools_depth
        self.output_dir = config.output_dir or "."

        self.CAZyme_annotation = os.path.join(self.input_dir, OVERVIEW_FILE)
        self.dbCANsub_substrate_annotation = os.path.join(self.input_dir, DBCAN_SUB_RESULT_FILE)
        self.PUL_substrate_annotation = os.path.join(self.input_dir, CGC_SUB_PREDICTION_FILE)
        self.PUL_annotation = os.path.join(self.input_dir, CGC_RESULT_FILE)

        self.output = self.out_path({
            "fam_abund": U.FAM_ABUND_FILE,
            "fam_substrate_abund": U.FAM_SUBSTRATE_ABUND_FILE,
            "CGC_abund": U.CGC_ABUND_FILE,
        }.get(function, U.CGC_SUBSTRATE_HOMOLOGY_FILE))

        self.parameters_check()

    def out_path(self, name: str) -> str:
        return os.path.join(self.output_dir, name)

    def parameters_check(self):
        def _require(path: str, desc: str):
            if not path or not os.path.exists(path):
                raise FileNotFoundError(f"{desc} does not exist: {path or '(not given)'}")

        if self.function == "fam_abund":
            logger.info("Estimating abundance of CAZyme")
            _require(self.CAZyme_annotation, "CAZyme overview")
        if self.function == "fam_substrate_abund":
            logger.info("Estimating abundance of Substrate according to dbCAN-sub")
            _require(self.dbCANsub_substrate_annotation, "dbCAN-sub result")
        if self.function == "CGC_abund":
            logger.info("Estimating abundance of CGC")
            _require(self.PUL_annotation, "CGC annotation")
        if self.function == "CGC_substrate_abund":
            logger.info("Estimating abundance of CGC substrate")
            _require(self.PUL_substrate_annotation, "CGC substrate prediction")
            _require(self.PUL_annotation, "CGC annotation")
        _require(self.bedtools, "gene read count file (-bt, from dbcan_utils cal_coverage)")
        os.makedirs(self.output_dir, exist_ok=True)


def get_length_readcount(seqid2dbcan_annotation: Dict[str, object],
                         seqid2readcount: Dict[str, bedtools_read_count],
                         annotation_desc: str = "annotation",
                         depth_file: str = "") -> int:
    """Attach length/read_count to annotation records. Returns the number of IDs found.

    Raises AbundanceError when none of the annotated proteins is in the count table: that is
    always an ID mismatch between the GFF given to cal_coverage and the run_dbcan output, and
    continuing would only write a table of zeros.
    """
    missing = []
    for seqid, ann in seqid2dbcan_annotation.items():
        rec = seqid2readcount.get(seqid)
        if rec is None:
            missing.append(seqid)
            ann.length = 0
            ann.read_count = 0.0
            continue
        ann.length = rec.length
        ann.read_count = rec.read_count
    found = len(seqid2dbcan_annotation) - len(missing)
    if seqid2dbcan_annotation and found == 0:
        raise AbundanceError(
            f"None of the {len(seqid2dbcan_annotation)} protein IDs of the {annotation_desc} is in the read count "
            f"table {depth_file} ({len(seqid2readcount)} genes). Example annotation IDs: {list(seqid2dbcan_annotation)[:3]}; "
            f"example count-table IDs: {list(seqid2readcount)[:3]}. Run cal_coverage with a GFF whose gene IDs are "
            f"the protein IDs used by run_dbcan (the cgc.gff in the run_dbcan output directory always is), or "
            f"choose the GFF attribute with --id-attr."
        )
    if missing:
        logger.warning(f"{len(missing)} of {len(seqid2dbcan_annotation)} annotated proteins are not in the read "
                       f"count table and get abundance 0, examples: {missing[:5]}")
    return found


class CAZyme_Abundance_estimate:
    def __init__(self, parameters: AbundParameters):
        self.pars = parameters
        seqid2readcount, normalized_tpm, total = ReadBedtools(parameters.bedtools)
        if not seqid2readcount:
            raise AbundanceError(f"Read count table {parameters.bedtools} has no gene rows (expected: id, length, count)")
        if total <= 0:
            raise AbundanceError(f"Read count table {parameters.bedtools} has no counted reads (all {len(seqid2readcount)} "
                                 f"genes are 0); check the BAM, the GFF contig names and the cal_coverage filters")
        self.fq_reads_count = total
        self.normalized_tpm = normalized_tpm
        logger.info(f"Total count assigned to genes: {total}")

        if parameters.function == "fam_abund":
            seqid2dbcan_annotation = ReadOverView(parameters.CAZyme_annotation)
            desc = OVERVIEW_FILE
        elif parameters.function == "fam_substrate_abund":
            seqid2dbcan_annotation = Read_dbcansub_out(parameters.dbCANsub_substrate_annotation)
            desc = DBCAN_SUB_RESULT_FILE
        elif parameters.function in ("CGC_abund", "CGC_substrate_abund"):
            seqid2dbcan_annotation, self.cgcid2cgc_standard = Read_cgc_standard_out(parameters.PUL_annotation)
            desc = CGC_RESULT_FILE
            if parameters.function == "CGC_substrate_abund":
                self.cgcid2cgc_substrate = Read_cgc_substrate(parameters.PUL_substrate_annotation)
        else:
            raise ValueError(f"Unknown abundance function: {parameters.function}")

        get_length_readcount(seqid2dbcan_annotation, seqid2readcount, desc, parameters.bedtools)
        self.seqid2dbcan_annotation = seqid2dbcan_annotation

    def Cal_Seq_Abundance(self, method: str = "FPKM"):
        method = method.upper()
        if method in ("FPKM", "RPKM"):
            norm_reads = self.fq_reads_count / 1e6
            for ann in self.seqid2dbcan_annotation.values():
                length_kb = ann.length / 1000.0
                ann.abund = ann.read_count / norm_reads / length_kb if length_kb > 0 else 0.0
        elif method == "RPM":
            norm_reads = self.fq_reads_count / 1e6
            for ann in self.seqid2dbcan_annotation.values():
                ann.abund = ann.read_count / norm_reads
        elif method == "TPM":
            denom = max(self.normalized_tpm, 1e-12)
            for ann in self.seqid2dbcan_annotation.values():
                ann.abund = (ann.read_count / ann.length * 1e6) / denom if ann.length > 0 else 0.0
        else:
            raise ValueError(f"Unknown abundance method: {method}")

    def Cal_Famliy_Abundance(self):
        family2seqid: Dict[str, List[str]] = {}
        for seqid, ann in self.seqid2dbcan_annotation.items():
            for fam in getattr(ann, "preds", []):
                family2seqid.setdefault(fam, []).append(seqid)
        self.family2seqid = family2seqid
        self.family2abund = {
            fam: float(sum(self.seqid2dbcan_annotation[sid].abund for sid in seqs))
            for fam, seqs in family2seqid.items()
        }

    def Cal_Substrate_Abundance(self):
        substrate2seqid: Dict[str, List[str]] = {}
        for seqid, ann in self.seqid2dbcan_annotation.items():
            for sub in getattr(ann, "substrates", []):
                substrate2seqid.setdefault(sub, []).append(seqid)
        self.substrate2seqid = substrate2seqid
        self.substrate2abund = {
            sub: float(sum(self.seqid2dbcan_annotation[sid].abund for sid in seqs))
            for sub, seqs in substrate2seqid.items()
        }

    def Cal_PUL_Abundance(self):
        """CGC abundance = mean abundance of its protein-coding members (file order kept)."""
        cgcid2seqid: Dict[str, List[str]] = {}
        cgcid2seqabund: Dict[str, List[float]] = {}
        for cgcid, records in self.cgcid2cgc_standard.items():
            ids = [rec.seqid for rec in records]
            cgcid2seqid[cgcid] = ids
            cgcid2seqabund[cgcid] = [float(self.seqid2dbcan_annotation[i].abund) for i in ids]
        self.cgcid2seqid = cgcid2seqid
        self.cgcid2seqabund = cgcid2seqabund
        self.cgcid2abund = {k: (float(np.mean(v)) if v else 0.0) for k, v in cgcid2seqabund.items()}

    def Cal_PUL_Substrate_Abundance(self):
        cgcsubstrate2cgcid_homo: Dict[str, List[str]] = {}
        cgcsubstrate2cgcid_major: Dict[str, List[str]] = {}
        for cgcid, rec in self.cgcid2cgc_substrate.items():
            homo = (rec.homo_sub or "").strip()
            if homo and homo not in ("X", "-"):
                cgcsubstrate2cgcid_homo.setdefault(homo, []).append(cgcid)
            for s in split_substrates(rec.major_voting_sub):
                cgcsubstrate2cgcid_major.setdefault(s, []).append(cgcid)

        def _abunds(sub2cgcs: Dict[str, List[str]]) -> Dict[str, List[float]]:
            out: Dict[str, List[float]] = {}
            for sub, cgcs in sub2cgcs.items():
                for c in cgcs:
                    if c in self.cgcid2abund:
                        out.setdefault(sub, []).append(self.cgcid2abund[c])
                    else:
                        logger.warning(f"CGC {c} of {CGC_SUB_PREDICTION_FILE} is not in {CGC_RESULT_FILE}; skipped")
            return out

        # Only CGCs with an abundance are listed, so cgcs and cgcs_abunds stay aligned.
        self.cgcsubstrate2abunds_homo = _abunds(cgcsubstrate2cgcid_homo)
        self.cgcsubstrate2cgcid_homo = {s: [c for c in cs if c in self.cgcid2abund] for s, cs in cgcsubstrate2cgcid_homo.items()}
        self.cgcsubstrate2abunds_major_votting = _abunds(cgcsubstrate2cgcid_major)
        self.cgcsubstrate2cgcid_major_votting = {s: [c for c in cs if c in self.cgcid2abund] for s, cs in cgcsubstrate2cgcid_major.items()}

    @staticmethod
    def _order(values: List[float]) -> List[int]:
        # stable sort, highest first
        return sorted(range(len(values)), key=lambda i: -values[i])

    def _write_cgc_substrate(self, path: str, sub2abunds: Dict[str, List[float]], sub2cgcs: Dict[str, List[str]]):
        subs = list(sub2abunds)
        sums = [float(np.sum(sub2abunds[s])) for s in subs]
        with open(path, "w") as f:
            f.write("Substrate\tAbundance(sum of CGC)\tcgcs\tcgcs_abunds\n")
            for i in self._order(sums):
                s = subs[i]
                f.write(f"{s}\t{round(sums[i], 3)}\t{';'.join(sub2cgcs.get(s, []))}\t"
                        f"{';'.join(str(round(a, 3)) for a in sub2abunds[s])}\n")

    def output_cgcsubstrate_abund(self):
        homo = self.pars.out_path(U.CGC_SUBSTRATE_HOMOLOGY_FILE)
        major = self.pars.out_path(U.CGC_SUBSTRATE_MAJORITY_FILE)
        logger.info(f"Writing CGC substrate abundance -> {homo}, {major}")
        self._write_cgc_substrate(homo, self.cgcsubstrate2abunds_homo, self.cgcsubstrate2cgcid_homo)
        self._write_cgc_substrate(major, self.cgcsubstrate2abunds_major_votting, self.cgcsubstrate2cgcid_major_votting)

    def output_family_abund(self):
        fams = list(self.family2abund)
        abunds = [self.family2abund[f] for f in fams]
        fam_path = self.pars.output
        subfam_path = self.pars.out_path(U.SUBFAM_ABUND_FILE)
        ec_path = self.pars.out_path(U.EC_ABUND_FILE)
        logger.info(f"Writing family/subfamily/EC abundance -> {fam_path}, {subfam_path}, {ec_path}")
        with open(fam_path, "w") as fam_file, open(subfam_path, "w") as subfam_file, open(ec_path, "w") as ec_file:
            fam_file.write("Family\tAbundance\tSeqNum\n")
            subfam_file.write("Subfamily\tAbundance\tSeqNum\n")
            ec_file.write("EC\tAbundance\tSeqNum\n")
            for i in self._order(abunds):
                famid = fams[i]
                line = f"{famid}\t{round(abunds[i], 3)}\t{len(self.family2seqid.get(famid, []))}\n"
                if Is_EC(famid):
                    ec_file.write(line)
                elif "_e" in famid:
                    subfam_file.write(line)
                else:
                    fam_file.write(line)

    def output_substrate_abund(self):
        subs = list(self.substrate2abund)
        abunds = [self.substrate2abund[s] for s in subs]
        logger.info(f"Writing substrate abundance -> {self.pars.output}")
        with open(self.pars.output, "w") as f:
            f.write("Substrate\tAbundance\tGeneID\n")
            for i in self._order(abunds):
                f.write(f"{subs[i]}\t{round(abunds[i], 3)}\t{';'.join(self.substrate2seqid.get(subs[i], []))}\n")

    def output_cgc_abund(self):
        cgcids = list(self.cgcid2abund)
        abunds = [self.cgcid2abund[c] for c in cgcids]
        logger.info(f"Writing CGC abundance -> {self.pars.output}")
        with open(self.pars.output, "w") as f:
            f.write("#CGCID\tAbundance(mean)\tSeqid\tSeq_abund\tFams\n")
            for i in self._order(abunds):
                cgcid = cgcids[i]
                records = self.cgcid2cgc_standard.get(cgcid, [])
                fams = ";".join(rec.protfam if rec.gene_type == "CAZyme" else rec.gene_type for rec in records)
                f.write(f"{cgcid}\t{round(abunds[i], 3)}\t{';'.join(self.cgcid2seqid.get(cgcid, []))}\t"
                        f"{';'.join(str(round(a, 3)) for a in self.cgcid2seqabund.get(cgcid, []))}\t{fams}\n")


def run_abundance(cfg: AbundanceConfig, function: str) -> None:
    pars = AbundParameters(cfg, function)
    est = CAZyme_Abundance_estimate(pars)
    est.Cal_Seq_Abundance(cfg.abundance)
    if function == "fam_abund":
        est.Cal_Famliy_Abundance()
        est.output_family_abund()
    elif function == "fam_substrate_abund":
        est.Cal_Substrate_Abundance()
        est.output_substrate_abund()
    elif function == "CGC_abund":
        est.Cal_PUL_Abundance()
        est.output_cgc_abund()
    else:
        est.Cal_PUL_Abundance()
        est.Cal_PUL_Substrate_Abundance()
        est.output_cgcsubstrate_abund()


# ---------------------------------------------------------------------------
# cal_coverage: GFF -> genes
# ---------------------------------------------------------------------------
@dataclass
class GeneRegion:
    """One gene: every GFF row with the same ID (e.g. the CDS exons of one protein)."""
    seqid: str
    strand: str
    # 0-based half-open intervals, merged, per contig
    intervals: List[Tuple[str, int, int]] = field(default_factory=list)
    read_count: int = 0

    @property
    def length(self) -> int:
        return sum(e - s for _, s, e in self.intervals)

    def add(self, contig: str, start1: int, end1: int) -> None:
        """Add a 1-based inclusive GFF interval."""
        self.intervals.append((contig, start1 - 1, end1))

    def merge(self) -> None:
        merged: List[Tuple[str, int, int]] = []
        for contig, s, e in sorted(self.intervals):
            if merged and merged[-1][0] == contig and s <= merged[-1][2]:
                merged[-1] = (contig, merged[-1][1], max(merged[-1][2], e))
            else:
                merged.append((contig, s, e))
        self.intervals = merged


def _gff_rows(gff: str) -> Iterable[Tuple[str, str, int, int, str, Dict[str, str]]]:
    with open(gff, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("##FASTA") or line.startswith(">"):
                break
            if not line.strip() or line.startswith("#"):
                continue
            p = line.rstrip("\n\r").split("\t")
            if len(p) < 9:
                continue
            if p[2].strip().lower() in U.GFF_SKIP_FEATURE_TYPES:
                continue
            try:
                start, end = int(p[3]), int(p[4])
            except ValueError:
                continue
            yield p[0], p[2], min(start, end), max(start, end), p[6], _split_gff_attributes(p[8])


def read_gene_regions(gff: str, id_attr: str = U.GFF_ID_ATTR_AUTO,
                      feature_type: Optional[str] = None) -> List[GeneRegion]:
    """Genes to count from a GFF, keyed by the protein ID that run_dbcan uses.

    * ``id_attr='auto'``: if any row carries ``protein_id`` (NCBI GFF CDS rows, run_dbcan's
      cgc.gff) or ``proteinId`` (JGI), only those rows are used and keyed by it; otherwise each
      row is keyed by ``ID``, then ``Name``, then ``locus_tag`` (Prodigal/pyrodigal GFF).
    * When the chosen rows include CDS rows only CDS rows are used (gene/mRNA/exon rows of the
      same locus would otherwise be counted again under another ID).
    * Rows with the same ID are one gene; the gene length is the union of its intervals.
    * Placeholder IDs ('NA', 'unknown') of non-coding features are skipped.
    """
    types = {t.strip() for t in feature_type.split(",") if t.strip()} if feature_type else None
    rows = [r for r in _gff_rows(gff) if types is None or r[1] in types]
    if not rows:
        raise AbundanceError(f"No usable feature rows in GFF {gff}" + (f" with type(s) {sorted(types)}" if types else ""))

    if id_attr and id_attr != U.GFF_ID_ATTR_AUTO:
        keys: Tuple[str, ...] = (id_attr,)
    else:
        keys = next(((k,) for k in U.GFF_ID_ATTR_PRIMARY
                     if any(not is_placeholder_id(r[5].get(k)) for r in rows)), U.GFF_ID_ATTR_FALLBACK)

    def _key(attrs: Dict[str, str]) -> Optional[str]:
        for k in keys:
            v = attrs.get(k)
            if not is_placeholder_id(v):
                return v.strip()
        return None

    keyed = [(r, _key(r[5])) for r in rows]
    keyed = [(r, k) for r, k in keyed if k is not None]
    if types is None and any(r[1] == U.GFF_CDS_FEATURE for r, _ in keyed):
        keyed = [(r, k) for r, k in keyed if r[1] == U.GFF_CDS_FEATURE]
    if not keyed:
        raise AbundanceError(
            f"No feature in GFF {gff} has a usable gene ID (tried attribute(s) {', '.join(keys)}); "
            f"set the attribute that holds the protein ID with --id-attr")

    genes: Dict[str, GeneRegion] = {}
    for (contig, _ftype, start, end, strand, _attrs), gid in keyed:
        g = genes.get(gid)
        if g is None:
            g = genes[gid] = GeneRegion(seqid=gid, strand=strand if strand in ("+", "-") else ".")
        g.add(contig, start, end)
    for g in genes.values():
        g.merge()
    logger.info(f"{len(genes)} genes from {len(keyed)} GFF rows (gene ID attribute: {'/'.join(keys)})")
    return list(genes.values())


# ---------------------------------------------------------------------------
# cal_coverage: alignment filters
# ---------------------------------------------------------------------------
_CIGAR_MATCH_OPS = (0, 7, 8)  # M, =, X
_CIGAR_INS, _CIGAR_DEL, _CIGAR_DIFF = 1, 2, 8


def cal_identity_based_on_MDtag(MDtag: List[Tuple[str, str]]) -> float:
    """Identity from an MD tag alone (matches / (matches + mismatches + deleted bases)); -1 if absent."""
    md_val = next((val for tag, val in MDtag if tag == "MD"), "")
    if not md_val:
        return -1.0
    match_base = mismatch_base = 0
    for token in ("".join(g) for _, g in groupby(md_val, key=lambda x: x.isdigit())):
        if token.isdigit():
            match_base += int(token)
        elif token.startswith("^"):
            mismatch_base += max(len(token) - 1, 0)
        else:
            mismatch_base += len(token)
    denom = match_base + mismatch_base
    return float(match_base) / denom if denom > 0 else -1.0


def _md_edits(md: str) -> int:
    """Mismatched + deleted bases encoded in an MD tag."""
    edits = 0
    for token in ("".join(g) for _, g in groupby(md, key=lambda x: x.isdigit())):
        if not token.isdigit():
            edits += len(token) - 1 if token.startswith("^") else len(token)
    return edits


def alignment_identity(read: pysam.AlignedSegment) -> Optional[float]:
    """Identity over the alignment columns (M/=/X + I + D): 1 - edits / columns.

    Edits come from the NM tag, else from MD (+ inserted bases from the CIGAR), else from
    '='/'X' CIGAR operations. Returns None when the alignment carries none of these.
    """
    cig = read.cigartuples
    if not cig:
        return None
    ops = Counter()
    for op, n in cig:
        ops[op] += n
    matched = sum(ops[o] for o in _CIGAR_MATCH_OPS)
    cols = matched + ops[_CIGAR_INS] + ops[_CIGAR_DEL]
    if cols <= 0:
        return None
    if read.has_tag("NM"):
        edits = int(read.get_tag("NM"))
    elif read.has_tag("MD"):
        edits = _md_edits(str(read.get_tag("MD"))) + ops[_CIGAR_INS]
    elif ops[0] == 0 and (ops[7] or ops[8]):
        edits = ops[_CIGAR_DIFF] + ops[_CIGAR_INS] + ops[_CIGAR_DEL]
    else:
        return None
    return max(0.0, 1.0 - edits / cols)


def _read_on_gene_strand(read: pysam.AlignedSegment, gene_strand: str, strandedness: str) -> bool:
    if strandedness == U.STRANDEDNESS_NONE or gene_strand not in ("+", "-"):
        return True
    minus = read.is_reverse
    if read.is_paired and read.is_read2:
        minus = not minus
    if strandedness == U.STRANDEDNESS_REVERSE:
        minus = not minus
    return minus == (gene_strand == "-")


def read_passes(cfg: CoverageConfig, read: pysam.AlignedSegment, gene: GeneRegion, stats: Counter) -> bool:
    """Whether one alignment counts for a gene. Rejection reasons are tallied in ``stats``."""
    if read.is_unmapped:
        return False
    if read.is_secondary or read.is_supplementary:
        stats["secondary_or_supplementary"] += 1
        return False
    if read.is_qcfail:
        stats["qc_fail"] += 1
        return False
    if read.is_duplicate and not cfg.include_duplicates:
        stats["duplicate"] += 1
        return False
    if read.mapping_quality < cfg.mapping_quality:
        stats["low_mapq"] += 1
        return False
    if not _read_on_gene_strand(read, gene.strand, cfg.strandedness):
        stats["wrong_strand"] += 1
        return False
    overlap = 0
    for contig, s, e in gene.intervals:
        if contig == read.reference_name:
            overlap += read.get_overlap(s, e) or 0
    if overlap <= 0:
        return False
    # read length including soft/hard clips (as before); long reads: relative to the shorter of read and gene
    read_len = read.infer_read_length() or read.query_length or read.query_alignment_length or 1
    denom = min(gene.length, read_len) if cfg.hifi else read_len
    if overlap / max(denom, 1) < cfg.overlap_base_ratio:
        stats["low_overlap"] += 1
        return False
    if cfg.identity > 0:
        ident = alignment_identity(read)
        if ident is None:
            stats["no_identity_tag"] += 1  # identity filter not applicable; read is kept
        elif ident < cfg.identity:
            stats["low_identity"] += 1
            return False
    return True


def count_gene(cfg: CoverageConfig, bam: pysam.AlignmentFile, gene: GeneRegion, refs: set, stats: Counter) -> int:
    """Reads (or fragments: both mates of a pair count once) assigned to one gene."""
    fragments = cfg.count_mode == U.COUNT_MODE_FRAGMENTS
    seen = set()
    for contig, s, e in gene.intervals:
        if contig not in refs:
            continue
        for read in bam.fetch(contig, s, e):
            if fragments and read.is_paired:
                key = read.query_name
            else:
                key = (read.query_name, read.is_read1, read.is_read2)
            if key in seen:
                continue
            if read_passes(cfg, read, gene, stats):
                seen.add(key)
    return len(seen)


def accomplete_function(cfg: CoverageConfig, genes: List[GeneRegion], m: int = 0,
                        samfile_name: Optional[str] = None) -> Tuple[List[GeneRegion], Counter]:
    stats: Counter = Counter()
    with pysam.AlignmentFile(samfile_name or cfg.input, "rb") as bam:
        refs = set(bam.references)
        for gene in tqdm(genes, desc="Counting reads per gene", disable=m > 0 or len(genes) < 1000):
            gene.read_count = count_gene(cfg, bam, gene, refs, stats)
    return genes, stats


def slice_list(items: List, n: int) -> List[List]:
    """Split a list into n near-equal chunks."""
    n = max(1, int(n))
    k, r = divmod(len(items), n)
    out, start = [], 0
    for i in range(n):
        end = start + k + (1 if i < r else 0)
        out.append(items[start:end])
        start = end
    return out


def cal_coverage(cfg: CoverageConfig) -> List[GeneRegion]:
    """Count reads (default: fragments) per gene; write 'id<TAB>length<TAB>count' to cfg.output."""
    if cfg.count_mode not in U.COUNT_MODES:
        raise ValueError(f"count_mode must be one of {U.COUNT_MODES}")
    if cfg.strandedness not in U.STRANDEDNESS_CHOICES:
        raise ValueError(f"strandedness must be one of {U.STRANDEDNESS_CHOICES}")
    genes = read_gene_regions(cfg.gff, cfg.id_attr, cfg.feature_type)

    with pysam.AlignmentFile(cfg.input, "rb") as bam:
        if not bam.has_index():
            raise AbundanceError(f"BAM {cfg.input} has no index; run 'samtools sort' and 'samtools index' first")
        refs = set(bam.references)
    absent = [g.seqid for g in genes if not any(c in refs for c, _, _ in g.intervals)]
    if len(absent) == len(genes):
        contigs = sorted({c for g in genes for c, _, _ in g.intervals})
        raise AbundanceError(
            f"None of the GFF contigs is a reference of BAM {cfg.input} (GFF: {contigs[:3]}, BAM: {sorted(refs)[:3]}); "
            f"the reads must be mapped to the same contigs the genes were predicted on")
    if absent:
        logger.warning(f"{len(absent)} genes lie on contigs absent from the BAM header and get count 0, examples: {absent[:5]}")

    stats: Counter = Counter()
    threads = max(1, int(cfg.threads or 1))
    if threads >= 2 and len(genes) > 1:
        chunks = [c for c in slice_list(genes, threads) if c]
        logger.info(f"Counting reads with {len(chunks)} worker processes")
        with Pool(len(chunks)) as pool:
            jobs = [pool.apply_async(accomplete_function, (cfg, chunk, i + 1)) for i, chunk in enumerate(chunks)]
            genes = []
            for job in jobs:
                chunk_genes, chunk_stats = job.get()
                genes.extend(chunk_genes)
                stats.update(chunk_stats)
    else:
        genes, stats = accomplete_function(cfg, genes)

    out_dir = os.path.dirname(os.path.abspath(cfg.output))
    os.makedirs(out_dir, exist_ok=True)
    with open(cfg.output, "w") as fh:
        for g in genes:
            fh.write(f"{g.seqid}\t{g.length}\t{g.read_count}\n")

    total = sum(g.read_count for g in genes)
    unit = "fragments" if cfg.count_mode == U.COUNT_MODE_FRAGMENTS else "reads"
    logger.info(f"Counted {total} {unit} on {sum(1 for g in genes if g.read_count)} of {len(genes)} genes -> {cfg.output}")
    if stats:
        logger.info("Alignments rejected (per gene tested): " + ", ".join(f"{k}={v}" for k, v in sorted(stats.items())))
    if stats.get("no_identity_tag"):
        logger.warning(f"{stats['no_identity_tag']} alignments carry no NM/MD tag and no =/X CIGAR, so the --identity "
                       f"filter was not applied to them (they were counted). Add tags with 'samtools calmd' to filter them.")
    if total == 0:
        logger.warning("No read passed the filters; every gene count is 0")
    return genes


# ---------------------------------------------------------------------------
# gff_fix
# ---------------------------------------------------------------------------
def get_attributes_value(attribute: str, ID: str = "ID=") -> str:
    for attr in attribute.split(";"):
        attr = attr.strip()
        if attr.startswith(ID):
            return attr[len(ID):]
    return ""


def pep_fasta_analysis(filename: str) -> Dict[str, str]:
    ID2geneID: Dict[str, str] = {}
    for seq in SeqIO.parse(filename, 'fasta'):
        description = seq.description.split(";")[0]
        geneID = description.split()[0]
        tmpID = description.split()[-1].split("=")[-1]
        ID2geneID[tmpID] = geneID
    return ID2geneID


def gff_refine(input_faa: str, gff: str, output: Optional[str] = None) -> str:
    ID2geneID = pep_fasta_analysis(input_faa)
    if not output:
        output = os.path.splitext(gff)[0] + ".fix.gff"
    missing = 0
    with open(output, 'w') as out, open(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            lines = line.rstrip("\n").split("\t")
            tmpID = get_attributes_value(lines[-1], "ID=")
            geneID = ID2geneID.get(tmpID, "")
            if not geneID:
                missing += 1
                logger.debug(f"Cannot find real gene id of {tmpID} in {input_faa}")
            out.write(line.replace("ID=" + tmpID, "ID=" + (geneID or tmpID)))
    if missing:
        logger.warning(f"{missing} GFF rows kept their ID (not found in {input_faa})")
    logger.info(f"Writing fixed gff to file {output}")
    return output


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
@click.command()
@logging_options
@click.argument(
    "function",
    type=click.Choice(list(U.ABUND_FUNCTIONS) + ["cal_coverage", "gff_fix"])
)
@dbcan_utils_options
def cli(function, log_level, log_file, verbose, input_path, bedtools_file, gff_file, output_path, abundance,
        overlap_base_ratio, mapping_quality, identity, threads, hifi, count_mode, strandedness,
        include_duplicates, id_attr, feature_type, r1, r2, db_dir):
    """dbCAN abundance utilities

    # Functions:

    **cal_coverage**: Count reads per gene from a sorted, indexed BAM and a GFF

    Example: `dbcan_utils cal_coverage -g sample.dbCAN/cgc.gff -i sample.bam -o sample.depth.txt -t 6`

    **fam_abund**: Compute CAZyme family / subfamily / EC abundance

    Example: `dbcan_utils fam_abund -bt sample.depth.txt -i sample.dbCAN -a TPM -o sample_abund`

    **fam_substrate_abund**: Compute substrate abundance (dbCAN-sub)

    Example: `dbcan_utils fam_substrate_abund -bt sample.depth.txt -i sample.dbCAN -a TPM -o sample_abund`

    **CGC_abund**: Compute CGC/PUL abundance

    Example: `dbcan_utils CGC_abund -bt sample.depth.txt -i sample.dbCAN -a TPM -o sample_abund`

    **CGC_substrate_abund**: Compute CGC substrate abundance

    Example: `dbcan_utils CGC_substrate_abund -bt sample.depth.txt -i sample.dbCAN -a TPM -o sample_abund`

    **gff_fix**: Fix Prodigal GFF IDs (ID -> geneID)

    Example: `dbcan_utils gff_fix -i UT30.3.faa -g UT30.3.gff`
    """
    setup_logging(log_level, log_file, verbose)
    for name, value in (("-1/--r1", r1), ("-2/--r2", r2), ("--db-dir", db_dir)):
        if value:
            logger.warning(f"{name} is not used by dbcan_utils and is ignored (deprecated option)")

    try:
        if function in U.ABUND_FUNCTIONS:
            if not bedtools_file:
                raise click.UsageError(f"{function} requires -bt/--bedtools (the cal_coverage output)")
            input_dir = input_path if os.path.isdir(input_path) else (os.path.dirname(os.path.abspath(input_path)) or ".")
            cfg = AbundanceConfig(
                input_dir=input_dir,
                bedtools_depth=bedtools_file,
                output_dir=output_path or ".",
                abundance=abundance,
            )
            run_abundance(cfg, function)

        elif function == "cal_coverage":
            if not gff_file or not os.path.exists(gff_file):
                raise click.UsageError("cal_coverage requires an existing -g/--gff")
            if not os.path.exists(input_path):
                raise click.UsageError(f"BAM file not found: {input_path}")
            if not output_path:
                raise click.UsageError("cal_coverage requires -o/--output (the read count table to write)")
            cal_coverage(CoverageConfig(
                input=input_path,
                gff=gff_file,
                output=output_path,
                overlap_base_ratio=overlap_base_ratio,
                mapping_quality=mapping_quality,
                identity=identity,
                threads=threads,
                hifi=hifi,
                count_mode=count_mode,
                strandedness=strandedness,
                include_duplicates=include_duplicates,
                id_attr=id_attr,
                feature_type=feature_type,
            ))

        elif function == "gff_fix":
            if not gff_file or not os.path.exists(gff_file) or not os.path.exists(input_path):
                raise click.UsageError("gff_fix requires -g/--gff and -i protein FASTA file")
            gff_refine(input_path, gff_file, output_path)
    except click.UsageError:
        raise
    except (AbundanceError, FileNotFoundError, ValueError, OSError) as e:
        logger.error(f"dbcan_utils {function} failed: {e}")
        raise click.ClickException(str(e))


if __name__ == "__main__":
    cli()
