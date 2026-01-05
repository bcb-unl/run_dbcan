from Bio import SeqIO
from itertools import (takewhile, repeat, groupby)
from multiprocessing import Pool
import rich_click as click

# 配置rich_click以优化显示
click.rich_click.USE_RICH_MARKUP = True
click.rich_click.USE_MARKDOWN = True
click.rich_click.SHOW_ARGUMENTS = True
click.rich_click.GROUP_ARGUMENTS_OPTIONS = True
click.rich_click.STYLE_ERRORS_SUGGESTION = "magenta italic"
click.rich_click.SHOW_METAVARS_COLUMN = False
click.rich_click.APPEND_METAVARS_HELP = True

import gzip
import logging
import os
from types import SimpleNamespace
from typing import Dict, List, Tuple

import numpy as np
import pysam
from tqdm import tqdm

from dbcan.constants.utils_constants import (
    CGC_RESULT_FILE,
    CGC_SUB_PREDICTION_FILE,
    DBCAN_SUB_RESULT_FILE,
    OVERVIEW_FILE,
)
from dbcan.configs.utils_config import AbundanceConfig
from dbcan.parameter import logging_options
from dbcan.main import setup_logging

logger = logging.getLogger(__name__)


def fq_file_line_count(file_name: str) -> int:
    """Estimate number of reads in FASTQ (line_count/4). Support .gz and plain text."""
    def _count_lines_text(path: str) -> int:
        bufsize = 1024 * 1024
        total = 0
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            for buf in takewhile(lambda x: x, (f.read(bufsize) for _ in repeat(None))):
                total += buf.count("\n")
        return total

    def _count_lines_gz(path: str) -> int:
        bufsize = 1024 * 1024
        total = 0
        with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as f:
            for buf in takewhile(lambda x: x, (f.read(bufsize) for _ in repeat(None))):
                total += buf.count("\n")
        return total

    try:
        if file_name.endswith(".gz"):
            lines = _count_lines_gz(file_name)
        else:
            lines = _count_lines_text(file_name)
        return lines // 4
    except Exception as e:
        logger.error(f"FASTQ line counting failed: {file_name}, {e}", exc_info=True)
        return 0


def total_mapped_reads_count(file_name: str) -> float:
    """Sum the third column in bedtools/cal_coverage-like output as total mapped reads."""
    total = 0.0
    if not os.path.exists(file_name):
        logger.warning(f"Read count file does not exist: {file_name}")
        return total
    try:
        with open(file_name) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                try:
                    total += float(parts[-1])
                except Exception:
                    continue
    except Exception as e:
        logger.error(f"Failed to read total mapped reads: {file_name}, {e}", exc_info=True)
    return total


class AbundParameters:
    """Derive and validate parameters (driven by AbundanceConfig)."""
    def __init__(self, config: AbundanceConfig, function: str):
        self.config = config
        self.function = function

        self.input_dir = config.input_dir
        self.bedtools = config.bedtools_depth

        # Use filenames from constants
        self.CAZyme_annotation = os.path.join(self.input_dir, OVERVIEW_FILE)
        self.dbCANsub_substrate_annotation = os.path.join(self.input_dir, DBCAN_SUB_RESULT_FILE)
        self.PUL_substrate_annotation = os.path.join(self.input_dir, CGC_SUB_PREDICTION_FILE)
        self.PUL_annotation = os.path.join(self.input_dir, CGC_RESULT_FILE)

        self.output = {
            "fam_abund": "fam_abund.out",
            "fam_substrate_abund": "fam_substrate_abund.out",
            "CGC_abund": "CGC_abund.out",
            "CGC_substrate_abund": "CGC_substrate"
        }.get(function, "output.out")

        self.parameters_check()

    def parameters_check(self):
        def _require(path: str, desc: str):
            if not os.path.exists(path):
                raise FileNotFoundError(f"{desc} does not exist: {path}")

        if self.function.startswith("fam_abund"):
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
        _require(self.bedtools, "gene read count file")


class bedtools_read_count:
    def __init__(self, lines: List[str]):
        self.seqid = lines[0]
        self.length = int(float(lines[1]))
        self.read_count = float(lines[2])

    def __repr__(self):
        return "\t".join([str(getattr(self, k)) for k in vars(self)])


def ReadBedtools(filename: str) -> Tuple[Dict[str, bedtools_read_count], float]:
    seqid2info: Dict[str, bedtools_read_count] = {}
    normalized_tpm = 0.0
    with open(filename) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            rec = bedtools_read_count(parts[:3])
            seqid2info[rec.seqid] = rec
            if rec.length > 0:
                normalized_tpm += rec.read_count / rec.length
    return seqid2info, normalized_tpm


def Is_EC(ec: str) -> bool:
    return bool(ec) and ec[0].isdigit()


def Clean_Hmmer_sub(pred: str) -> List[str]:
    if not pred or pred == "-":
        return []
    preds = pred.split("+")
    clean = [p.split("(")[0] for p in preds if p]
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
        self.ECs = list(set([i.split(":")[0] for i in ec_field.split("|") if i and i.split(":")[0] != '-']))
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
        union_fam = list(set(self.hmmer) | set(self.dbcan_sub)) + self.ECs
        self.preds = union_fam if union_fam else []


def ReadOverView(filename: str) -> Dict[str, OverView]:
    """Parse overview.tsv by header; compatible with the added Substrate column."""
    with open(filename) as f:
        header = f.readline().rstrip("\n").split("\t")
        hmap = {h: i for i, h in enumerate(header)}
        # Try generic header names
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
                ec = parts[hmap.get(ec_col, 1)] if ec_col else "-"
                hmm = parts[hmap.get(hmm_col, 2)] if hmm_col else "-"
                sub = parts[hmap.get(sub_col, 3)] if sub_col else "-"
                dia = parts[hmap.get(dia_col, 4)] if dia_col else "-"
                out[seqid] = OverView(seqid, ec, hmm, sub, dia)
            except Exception:
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


def Read_dbcansub_out(filename: str) -> Dict[str, DbcanSub_line]:
    out: Dict[str, DbcanSub_line] = {}
    with open(filename) as f:
        next(f, None)
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            rec = DbcanSub_line(parts[:8])
            out[rec.seqid] = rec
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
    seqid2records: Dict[str, cgc_standard_line] = {}
    cgcid2records: Dict[str, List[cgc_standard_line]] = {}
    with open(filename) as f:
        next(f, None)
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            rec = cgc_standard_line(parts)
            seqid2records[rec.seqid] = rec
            cgcid2records.setdefault(rec.cgcid, []).append(rec)
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


def get_length_readcount(seqid2dbcan_annotation: Dict[str, object], seqid2readcount: Dict[str, bedtools_read_count]) -> None:
    missing = []
    for seqid in seqid2dbcan_annotation:
        read_count = seqid2readcount.get(seqid, None)
        if not read_count:
            missing.append(seqid)
            continue
        seqid_annotation = seqid2dbcan_annotation[seqid]
        seqid_annotation.length = read_count.length
        seqid_annotation.read_count = read_count.read_count
    if missing:
        logger.warning(f"{len(missing)} sequences missing read_count, examples: {missing[:5]}")


class CAZyme_Abundance_estimate:
    def __init__(self, parameters: AbundParameters):
        self.pars = parameters
        # Total reads (or total counts)
        self.fq_reads_count = total_mapped_reads_count(self.pars.bedtools)
        logger.info(f"Total reads count: {self.fq_reads_count}")

        seqid2readcount, normalized_tpm = ReadBedtools(parameters.bedtools)
        self.normalized_tpm = normalized_tpm

        # Load annotations
        if parameters.function == "fam_abund":
            seqid2dbcan_annotation = ReadOverView(parameters.CAZyme_annotation)
        elif parameters.function == "fam_substrate_abund":
            seqid2dbcan_annotation = Read_dbcansub_out(parameters.dbCANsub_substrate_annotation)
        elif parameters.function == "CGC_abund":
            seqid2dbcan_annotation, cgcid2cgc_standard = Read_cgc_standard_out(parameters.PUL_annotation)
            self.cgcid2cgc_standard = cgcid2cgc_standard
        elif parameters.function == "CGC_substrate_abund":
            seqid2dbcan_annotation, cgcid2cgc_standard = Read_cgc_standard_out(parameters.PUL_annotation)
            cgcid2cgc_substrate = Read_cgc_substrate(parameters.PUL_substrate_annotation)
            self.cgcid2cgc_standard = cgcid2cgc_standard
            self.cgcid2cgc_substrate = cgcid2cgc_substrate
        else:
            seqid2dbcan_annotation = {}

        get_length_readcount(seqid2dbcan_annotation, seqid2readcount)
        self.seqid2dbcan_annotation = seqid2dbcan_annotation

    def Cal_Seq_Abundance(self, method: str = "FPKM"):
        if method == "FPKM":
            norm_reads = max(self.fq_reads_count, 1.0) / 1e6
            for seqid, ann in self.seqid2dbcan_annotation.items():
                length_kb = max(getattr(ann, "length", 0) / 1000.0, 1e-9)
                ann.abund = getattr(ann, "read_count", 0.0) / norm_reads / length_kb
        elif method == "RPM":
            norm_reads = max(self.fq_reads_count, 1.0) / 1e6
            for seqid, ann in self.seqid2dbcan_annotation.items():
                ann.abund = getattr(ann, "read_count", 0.0) / norm_reads
        elif method == "TPM":
            denom = max(self.normalized_tpm, 1e-9)
            for seqid, ann in self.seqid2dbcan_annotation.items():
                ann.abund = (getattr(ann, "read_count", 0.0) / max(getattr(ann, "length", 0), 1e-9) * 1e6) / denom

    def Cal_Famliy_Abundance(self):
        family2seqid: Dict[str, List[str]] = {}
        for seqid, ann in self.seqid2dbcan_annotation.items():
            preds = getattr(ann, "preds", [])
            if preds:
                for fam in preds:
                    family2seqid.setdefault(fam, []).append(seqid)
        self.family2seqid = family2seqid
        family2abund = {fam: 0.0 for fam in family2seqid}
        for fam, seqs in family2seqid.items():
            for sid in seqs:
                family2abund[fam] += getattr(self.seqid2dbcan_annotation[sid], "abund", 0.0)
        self.family2abund = family2abund

    def Cal_Substrate_Abundance(self):
        substrate2seqid: Dict[str, List[str]] = {}
        for seqid, ann in self.seqid2dbcan_annotation.items():
            subs_raw = getattr(ann, "substrate", "-")
            subs = [s.strip() for s in subs_raw.replace("and", ",").split(",") if s and s.strip() and s.strip() != "-"]
            for sub in set(subs):
                substrate2seqid.setdefault(sub, []).append(seqid)

        substrate2abund = {sub: 0.0 for sub in substrate2seqid}
        for sub, seqs in substrate2seqid.items():
            for sid in seqs:
                substrate2abund[sub] += getattr(self.seqid2dbcan_annotation[sid], "abund", 0.0)
        self.substrate2abund = substrate2abund
        self.substrate2seqid = substrate2seqid

    def Cal_PUL_Abundance(self):
        cgcid2seqid: Dict[str, List[str]] = {}
        cgcid2seqabund: Dict[str, List[float]] = {}
        for seqid, ann in self.seqid2dbcan_annotation.items():
            cgcid = getattr(ann, "cgcid", None)
            if not cgcid:
                continue
            cgcid2seqid.setdefault(cgcid, []).append(seqid)
            cgcid2seqabund.setdefault(cgcid, []).append(getattr(ann, "abund", 0.0))

        cgcid2abund = {k: (np.mean(v) if v else 0.0) for k, v in cgcid2seqabund.items()}
        self.cgcid2abund = cgcid2abund
        self.cgcid2seqid = cgcid2seqid
        self.cgcid2seqabund = cgcid2seqabund

    def Cal_PUL_Substrate_Abundance(self):
        cgcsubstrate2cgcid_homo: Dict[str, List[str]] = {}
        cgcsubstrate2cgcid_major: Dict[str, List[str]] = {}
        for cgcid, rec in self.cgcid2cgc_substrate.items():
            if getattr(rec, "homo_sub", "") and getattr(rec, "homo_sub") != "X":
                cgcsubstrate2cgcid_homo.setdefault(rec.homo_sub, []).append(cgcid)
            if getattr(rec, "major_voting_sub", ""):
                for s in rec.major_voting_sub.split(","):
                    s = s.strip()
                    if s:
                        cgcsubstrate2cgcid_major.setdefault(s, []).append(cgcid)

        # Aggregate CGC abundance by substrate (use mean abundance per CGC)
        cgcsubstrate2abunds_homo: Dict[str, List[float]] = {}
        for sub, cgcs in cgcsubstrate2cgcid_homo.items():
            for c in cgcs:
                vals = self.cgcid2seqabund.get(c, [])
                if vals:
                    cgcsubstrate2abunds_homo.setdefault(sub, []).append(float(np.mean(vals)))
        self.cgcsubstrate2cgcid_homo = cgcsubstrate2cgcid_homo
        self.cgcsubstrate2abunds_homo = cgcsubstrate2abunds_homo

        cgcsubstrate2abunds_major: Dict[str, List[float]] = {}
        for sub, cgcs in cgcsubstrate2cgcid_major.items():
            for c in cgcs:
                vals = self.cgcid2seqabund.get(c, [])
                if vals:
                    cgcsubstrate2abunds_major.setdefault(sub, []).append(float(np.mean(vals)))
        self.cgcsubstrate2cgcid_major_votting = cgcsubstrate2cgcid_major
        self.cgcsubstrate2abunds_major_votting = cgcsubstrate2abunds_major

    def output_cgcsubstrate_abund(self):
        # Homology-based
        subs, sums, cgcs, cgcs_abund = [], [], [], []
        for sub, vals in self.cgcsubstrate2abunds_homo.items():
            subs.append(sub)
            sums.append(float(np.sum(vals)))
            cgcs.append(self.cgcsubstrate2cgcid_homo.get(sub, []))
            cgcs_abund.append(vals)
        order = np.argsort(sums)[::-1]
        logger.info("Writing CGC substrate abundance (homology) -> CGC_substrate_PUL_homology.out")
        with open("CGC_substrate_PUL_homology.out", "w") as f:
            f.write("Substrate\tAbundance(sum of CGC)\tcgcs\tcgcs_abunds\n")
            for i in order:
                f.write(f"{subs[i]}\t{round(sums[i],3)}\t{';'.join(cgcs[i])}\t{';'.join([str(round(a,3)) for a in cgcs_abund[i]])}\n")

        # Majority voting
        subs, sums, cgcs, cgcs_abund = [], [], [], []
        for sub, vals in self.cgcsubstrate2abunds_major_votting.items():
            subs.append(sub)
            sums.append(float(np.sum(vals)))
            cgcs.append(self.cgcsubstrate2cgcid_major_votting.get(sub, []))
            cgcs_abund.append(vals)
        order = np.argsort(sums)[::-1]
        logger.info("Writing CGC substrate abundance (majority) -> CGC_substrate_majority_voting.out")
        with open("CGC_substrate_majority_voting.out", "w") as f:
            f.write("Substrate\tAbundance(sum of CGC)\tcgcs\tcgcs_abunds\n")
            for i in order:
                f.write(f"{subs[i]}\t{round(sums[i],3)}\t{';'.join(cgcs[i])}\t{';'.join([str(round(a,3)) for a in cgcs_abund[i]])}\n")

    def output_family_abund(self, method: str = "family"):
        fams, abunds, seqs = [], [], []
        for fam, val in self.family2abund.items():
            fams.append(fam)
            abunds.append(val)
            seqs.append(self.family2seqid.get(fam, []))
        order = np.argsort(abunds)[::-1]
        logger.info(f"Writing family/subfamily/EC abundance -> fam_abund.out / subfam_abund.out / EC_abund.out")
        with open(self.pars.output, "w") as fam_file, open("subfam_abund.out", "w") as subfam_file, open("EC_abund.out", "w") as ec_file:
            fam_file.write("Family\tAbundance\tSeqNum\n")
            subfam_file.write("Subfamily\tAbundance\tSeqNum\n")
            ec_file.write("EC\tAbundance\tSeqNum\n")
            for i in order:
                famid = fams[i]
                if Is_EC(famid):
                    ec_file.write(f"{famid}\t{round(abunds[i],3)}\t{len(seqs[i])}\n")
                elif "_e" in famid:
                    subfam_file.write(f"{famid}\t{round(abunds[i],3)}\t{len(seqs[i])}\n")
                else:
                    fam_file.write(f"{famid}\t{round(abunds[i],3)}\t{len(seqs[i])}\n")

    def output_substrate_abund(self):
        subs, abunds, genes = [], [], []
        for sub, val in self.substrate2abund.items():
            subs.append(sub)
            abunds.append(val)
            genes.append(self.substrate2seqid.get(sub, []))
        order = np.argsort(abunds)[::-1]
        logger.info(f"Writing substrate abundance -> {self.pars.output}")
        with open(self.pars.output, "w") as f:
            f.write("Substrate\tAbundance\tGeneID\n")
            for i in order:
                if subs[i]:
                    f.write(f"{subs[i]}\t{round(abunds[i],3)}\t{';'.join(genes[i])}\n")

    def output_cgc_abund(self):
        cgcids, abunds, seqids, seq_abunds, cgc_records = [], [], [], [], []
        for cgcid, val in self.cgcid2abund.items():
            cgcids.append(cgcid)
            abunds.append(val)
            seqids.append(self.cgcid2seqid.get(cgcid, []))
            seq_abunds.append(self.cgcid2seqabund.get(cgcid, []))
            cgc_records.append(self.cgcid2cgc_standard.get(cgcid, []))
        order = np.argsort(abunds)[::-1]
        logger.info(f"Writing CGC abundance -> {self.pars.output}")
        with open(self.pars.output, "w") as f:
            f.write("#CGCID\tAbundance(mean)\tSeqid\tSeq_abund\tFams\n")
            for i in order:
                fams = ";".join(rec.protfam if getattr(rec, "gene_type", "") == "CAZyme" else getattr(rec, "gene_type", "") for rec in cgc_records[i])
                f.write(f"{cgcids[i]}\t{round(abunds[i],3)}\t{';'.join(seqids[i])}\t{';'.join([str(round(a,3)) for a in seq_abunds[i]])}\t{fams}\n")


def get_attributes_value(attribute: str, ID: str = "ID=") -> str:
    for attr in attribute.split(";"):
        if attr.startswith(ID):
            return attr[len(ID):]
    return ""


class GFF_record:
    def __init__(self, lines: List[str]):
        self.contig_id = lines[0]
        self.source = lines[1]
        self.type = lines[2]
        self.start = int(lines[3])
        self.end = int(lines[4])
        self.score = lines[5]
        self.strand = lines[6]
        self.phase = lines[7]
        self.attribute = lines[8]
        self.seqid = get_attributes_value(self.attribute)
        self.length = self.end - self.start + 1
        self.read_count = 0

    def __repr__(self):
        return "\t".join([str(getattr(self, k)) for k in vars(self)])


def cal_identity_based_on_MDtag(MDtag: List[Tuple[str, str]]) -> float:
    """Approximate alignment identity based on MD tag (ignoring indels).
    Fix deletion length: consecutive letters after '^' represent deletion length."""
    md_val = ""
    for tag, val in MDtag:
        if tag == "MD":
            md_val = val
            break
    if not md_val:
        return -1.0

    groups = [''.join(list(g)) for _, g in groupby(md_val, key=lambda x: x.isdigit())]
    match_base = 0
    mismatch_base = 0
    i = 0
    while i < len(groups):
        token = groups[i]
        if token.isdigit():
            match_base += int(token)
            i += 1
        else:
            # Handle deletion: starts with '^' followed by bases
            if token.startswith("^"):
                mismatch_base += max(len(token) - 1, 0)
            else:
                mismatch_base += len(token)
            i += 1
    denom = match_base + mismatch_base
    if denom <= 0:
        return -1.0
    return float(match_base) / denom


def justify_reads_alignment_properties(args, read, gene: GFF_record) -> bool:
    """Alignment filtering: overlap ratio / mapping quality / sequence identity."""
    overlap_base_numer = read.get_overlap(gene.start, gene.end)
    if not overlap_base_numer:
        return False

    # Overlap ratio for HiFi vs short reads
    if args.hifi:
        query_length = read.query_length or read.infer_read_length() or 0
        gene_len = min(gene.end - gene.start + 1, query_length) if query_length else (gene.end - gene.start + 1)
        denom = max(gene_len, 1)
        overlap_base_ratio = overlap_base_numer / denom
    else:
        denom = max(read.query_length or 0, 1)
        overlap_base_ratio = overlap_base_numer / denom

    if overlap_base_ratio < args.overlap_base_ratio:
        return False

    if read.mapping_quality < args.mapping_quality:
        return False

    sequence_identity = cal_identity_based_on_MDtag(read.get_tags())
    if sequence_identity < 0:
        logger.debug(f"MD tag missing: {read.query_name}")
        return False
    if sequence_identity < args.identity:
        return False
    return True


def cal_coverage(args):
    """Count reads per gene (bedtools-like)."""
    genes = [GFF_record(line.rstrip("\n").split("\t")) for line in open(args.gff) if line.strip() and not line.startswith("#")]
    with open(args.output, 'w') as coverage_file:
        if args.threads >= 2:
            results = multi_masks(args, genes, args.input)
            for async_res in results:
                for gene in async_res.get():
                    coverage_file.write(f"{gene.seqid}\t{gene.length}\t{gene.read_count}\n")
        else:
            samfile = pysam.AlignmentFile(args.input, "rb")
            for i in tqdm(range(len(genes)), desc="Processing gene"):
                gene = genes[i]
                reads = samfile.fetch(gene.contig_id, gene.start, gene.end)
                aligned_reads_num = [1 if justify_reads_alignment_properties(args, read, gene) else 0 for read in reads]
                gene.read_count = int(np.sum(aligned_reads_num))
            for gene in genes:
                coverage_file.write(f"{gene.seqid}\t{gene.length}\t{gene.read_count}\n")
    logger.info(f"Writing read count to file {args.output}.")


def accomplete_function(args, genes: List[GFF_record], m: int, samfile_name: str):
    samfile = pysam.AlignmentFile(samfile_name, "rb")
    # Disable tqdm in subprocess
    for gene in genes:
        reads = samfile.fetch(gene.contig_id, gene.start, gene.end)
        aligned_reads_num = [1 if justify_reads_alignment_properties(args, read, gene) else 0 for read in reads]
        gene.read_count = int(np.sum(aligned_reads_num))
    return genes


def slice_list(items: List, n: int) -> List[List]:
    """Split a list into n near-equal chunks."""
    n = max(1, int(n))
    L = len(items)
    k, r = divmod(L, n)
    out = []
    start = 0
    for i in range(n):
        end = start + k + (1 if i < r else 0)
        out.append(items[start:end])
        start = end
    return out


def multi_masks(paras, genes: List[GFF_record], samfile_name: str):
    chunks = slice_list(genes, paras.threads)
    logger.info(f"Parent process {os.getpid()}, spawning {paras.threads} workers")
    p = Pool(paras.threads)
    jobs = []
    for i, chunk in enumerate(chunks):
        jobs.append(p.apply_async(func=accomplete_function, args=(paras, chunk, i, samfile_name)))
    p.close()
    p.join()
    return jobs


def pep_fasta_analysis(filename: str) -> Dict[str, str]:
    ID2geneID: Dict[str, str] = {}
    for seq in SeqIO.parse(filename, 'fasta'):
        description = seq.description.split(";")[0]
        geneID = description.split()[0]
        tmpID = description.split()[-1].split("=")[-1]
        ID2geneID[tmpID] = geneID
    return ID2geneID


def gff_refine(args):
    ID2geneID = pep_fasta_analysis(args.input)
    gff_filename, _ = os.path.splitext(args.gff)
    with open(gff_filename + ".fix.gff", 'w') as gff:
        for line in open(args.gff):
            if line.startswith("#"):
                continue
            lines = line.rstrip("\n").split("\t")
            tmpID = get_attributes_value(lines[-1], "ID=")
            geneID = ID2geneID.get(tmpID, "")
            if not geneID:
                logger.warning(f"Cannot find real gene id of {tmpID} in file: {args.gff}")
            tmp_line = line.replace("ID=" + tmpID, "ID=" + (geneID or tmpID))
            gff.write(tmp_line)
    logger.info(f"Writing fixed gff to file {gff_filename}.fix.gff")



@click.command()
@logging_options
@click.argument(
    "function", 
    type=click.Choice([
        "fam_abund", "fam_substrate_abund", 
        "CGC_abund", "CGC_substrate_abund", 
        "cal_coverage", "gff_fix"
    ])
)
@click.option("-i", "--input", "input_path", required=True, 
              type=click.Path(exists=False), 
              help="run_dbcan output dir or file (depends on function).")
@click.option("-bt", "--bedtools", "bedtools_file", 
              type=click.Path(exists=False), 
              help="Gene read count file generated by cal_coverage.")
@click.option("-g", "--gff", "gff_file", 
              type=click.Path(exists=False), 
              help="Gene annotation GFF (for cal_coverage).")
@click.option("-1", "--r1", "r1", 
              type=click.Path(exists=False), 
              help="R1 reads (optional, supports .gz).")
@click.option("-2", "--r2", "r2", default=None, 
              type=click.Path(exists=False), 
              help="R2 reads (optional, supports .gz).")
@click.option("-o", "--output", "output_path", default="output", 
              show_default=True, 
              help="Output directory/file.")
@click.option("-a", "--abundance", "abundance", 
              default="RPM", 
              type=click.Choice(["FPKM", "RPM", "TPM"]), 
              show_default=True, 
              help="Normalization method.")
@click.option("--db-dir", "db_dir", default="db", 
              show_default=True, 
              help="dbCAN database directory.")
@click.option("--overlap-base-ratio", "--overlap_base_ratio", 
              default=0.2, type=float, 
              show_default=True, 
              help="Minimum read-gene overlap ratio.")
@click.option("--mapping-quality", "--mapping_quality", 
              default=30, type=int, 
              show_default=True, 
              help="Minimum mapping quality.")
@click.option("-c", "--identity", "identity", 
              default=0.98, type=float, 
              show_default=True, 
              help="Minimum alignment identity.")
@click.option("-t", "--threads", "threads", 
              default=1, type=int, 
              show_default=True, 
              help="Number of worker processes.")
@click.option("--hifi", "hifi", is_flag=True, 
              help="Input reads are HiFi.")
def cli(function, log_level, log_file, verbose, input_path, bedtools_file, gff_file, r1, r2, output_path, abundance, db_dir, overlap_base_ratio, mapping_quality, identity, threads, hifi):
    """dbCAN abundance utilities
    
    # Functions:
    
    **fam_abund**: Compute CAZyme abundance
    
    Example: `dbcan_utils fam_abund -bt samfiles/fefifo_8002_7.depth.txt -i fefifo_8002_7.dbCAN`
    
    **fam_substrate_abund**: Compute substrate abundance (dbCAN-sub)
    
    Example: `dbcan_utils fam_substrate_abund -bt samfiles/fefifo_8002_7.depth.txt -i fefifo_8002_7.dbCAN`
    
    **CGC_abund**: Compute CGC/PUL abundance
    
    Example: `dbcan_utils CGC_abund -bt samfiles/fefifo_8002_7.depth.txt -i fefifo_8002_7.dbCAN`
    
    **CGC_substrate_abund**: Compute CGC substrate abundance
    
    Example: `dbcan_utils CGC_substrate_abund -bt samfiles/fefifo_8002_7.depth.txt -i fefifo_8002_7.dbCAN`
    
    **cal_coverage**: Count reads per gene from a BAM and GFF
    
    Example: `dbcan_utils cal_coverage -g Wet2014.gff -i Wet2014.bam -o Wet2014.depth.txt -t 6`
    
    **gff_fix**: Fix Prodigal GFF IDs (ID -> geneID)
    
    Example: `dbcan_utils gff_fix -i UT30.3.faa -g UT30.3.gff`
    """
    setup_logging(log_level, log_file, verbose)
    # Build AbundanceConfig
    cfg = AbundanceConfig(
        input_dir=input_path if os.path.isdir(input_path) else os.path.dirname(os.path.abspath(input_path)) or ".",
        bedtools_depth=bedtools_file if bedtools_file else "",
        output_dir=output_path,
        gff=gff_file,
        R1=r1,
        R2=r2,
        db_dir=db_dir,
        abundance=abundance,
        overlap_base_ratio=overlap_base_ratio,
        mapping_quality=mapping_quality,
        identity=identity,
        threads=threads,
        hifi=hifi,
    )

    if function in {"fam_abund", "fam_substrate_abund", "CGC_abund", "CGC_substrate_abund"}:
        try:
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
        except Exception as e:
            logger.error(f"Abundance workflow failed: {e}", exc_info=True)

    elif function == "cal_coverage":
        if not gff_file or not os.path.exists(gff_file):
            logger.error("cal_coverage requires a valid --gff")
            return
        if not output_path:
            logger.error("cal_coverage requires --output")
            return
        args_ns = SimpleNamespace(
            input=input_path,
            gff=gff_file,
            output=output_path,
            overlap_base_ratio=overlap_base_ratio,
            mapping_quality=mapping_quality,
            identity=identity,
            threads=threads,
            hifi=hifi,
        )
        cal_coverage(args_ns)

    elif function == "gff_fix":
        if not gff_file or not os.path.exists(gff_file) or not input_path or not os.path.exists(input_path):
            logger.error("gff_fix requires --gff and -i protein fasta file")
            return
        args_ns = SimpleNamespace(input=input_path, gff=gff_file)
        gff_refine(args_ns)


if __name__ == "__main__":
    cli()
