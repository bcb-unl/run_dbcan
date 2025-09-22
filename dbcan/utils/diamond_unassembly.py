import os
import sys
from pathlib import Path
import gzip
import rich_click as click
from rich import print as rprint

# rich-click styling (can be adjusted)
click.rich_click.USE_RICH_MARKUP = True
click.rich_click.MAX_WIDTH = 100
click.rich_click.STYLE_HELPTEXT = "green"
click.rich_click.STYLE_METAVAR = "bold cyan"
click.rich_click.STYLE_OPTION = "bold yellow"
click.rich_click.STYLE_SWITCH = "bold yellow"
click.rich_click.STYLE_COMMAND = "bold magenta"
click.rich_click.STYLE_USAGE = "bold white"

def HLError(mess):
    return f"\033[1;31;40m{mess}:\033[0m"

'''
PAF record sample
1 string  Query sequence name
2 int     Query sequence length
3 int     Query start (0-based; BED-like; closed)
4 int     Query end (0-based; BED-like; open)
5 char    Relative strand: "+" or "-"
6 string  Target sequence name
7 int     Target sequence length
8 int     Target start on original strand (0-based)
9 int     Target end on original strand (0-based)
10 int    Number of residue matches
11 int    Alignment block length
12 int    Mapping quality (0-255; 255 for missing)
13 attr
'''

'''
DIAMOND (custom blast tabular)
1.  qseqid
2.  sseqid
3.  pident
4.  length
5.  mismatch
6.  gapopen
7.  qstart
8.  qend
9.  sstart
10. send
11. evalue
12. bitscore
13. qlen
14. slen
'''

def CAZy_filter(cazy):
    return set([aa for aa in cazy])

class PafRecord(object):
    def __init__(self, lines):
        self.Qsn = lines[0]
        self.Qsl = lines[12]
        self.Qs  = int(lines[6]) - 1
        self.Qe  = lines[7]
        self.Strand = lines[4]
        self.Tsn = lines[1]
        self.Tsl = lines[13]
        self.Ts  = int(lines[8]) - 1
        self.Te  = lines[9]
        self.Nrm = lines[11]
        self.Abl = lines[3]
        self.Mq  = lines[10]
        self.SeqID = self.Tsn.split('|')[0]
        self.CAZys = CAZy_filter(self.Tsn.strip("|").split("|")[1:])
        self.UniReadId = lines[0].split(".")[0]

    def __str__(self):
        return "\t".join([str(getattr(self, value)) for value in vars(self) if value != "CAZys"])

class Paf(object):
    def __init__(self, filename):
        self.records = [PafRecord(line.split()) for line in open(filename)]
    def __iter__(self):
        return iter(self.records)
    def GetReadId(self):
        return [record.Qsn for record in self]
    def GetSeqId(self):
        return [record.SeqID for record in self]
    def GetSeqLen(self):
        return {record.SeqID:record.Tsl for record in self}
    def CAZy2SeqID(self, CazySeqId):
        for record in self:
            for cazy in record.CAZys:
                CazySeqId.setdefault(cazy, []).append(record.SeqID)
    def SeqID2ReadID(self, aa):
        for record in self:
            aa.setdefault(record.SeqID, []).append(record.Qsn)
    def ReadID2Record(self):
        return {record.Qsn:record for record in self}
    def Output(self):
        [print(record) for record in self]
    def Assign_CAZy_megahit(self):
        for cazy in self:
            cazy.CAZys = CAZy_filter(cazy.Qsn.strip("|").split("|")[1:])
    def Assign_subfam(self, CAZyID2subfam):
        for hit in self:
            hit.subfams = CAZyID2subfam.get(hit.Tsn, "")
    def Get_subfam2SeqID(self, subfam2SeqID):
        for record in self:
            for cazy in record.subfams:
                subfam2SeqID.setdefault(cazy, []).append(record.SeqID)

def CAZyReadCount(cazyid, cazy2seqid, readtable):
    tmp_sum = 0
    for seqid in cazy2seqid[cazyid]:
        tmp_sum += readtable[seqid]
    return tmp_sum

def FPKMToCsv(args, tool, cazyfpkm, readtable, cazy2seqid):
    outfilename = args.output
    with open(outfilename, 'w') as f:
        f.write("Family\tAbundance\tSeqNum\tReadCount\n")
        for cazyid in cazyfpkm:
            seqnum = len(cazy2seqid[cazyid])
            readcount = CAZyReadCount(cazyid, cazy2seqid, readtable)
            fpkm = cazyfpkm[cazyid]
            if not cazyid[0].isdigit():
                f.write(f"{cazyid}\t{fpkm}\t{seqnum}\t{readcount}\n")

def check_read_type(filename):
    if filename.endswith("fq") or filename.endswith("fq.gz"):
        return "fq"
    elif filename.endswith("fa") or filename.endswith("fa.gz"):
        return "fa"
    else:
        sys.stderr.write(HLError("Error") + " File type not supported, please provide .fa(fa.gz) or (fq)fq.gz reads file.\n")
        exit(1)

def _count_fastq(path: Path) -> int:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, 'rt') as fh:
        # Count lines / 4
        lines = sum(1 for _ in fh)
    return lines // 4

def _count_fasta(path: Path) -> int:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, 'rt') as fh:
        return sum(1 for line in fh if line.startswith('>'))

def get_count_reads(file):
    p = Path(file)
    if p.suffix in [".gz"]:
        # Look at previous suffix
        stem_suffix = "".join(p.name.split(".")[-2:])
        # Not strictly needed; rely on pattern
    if file.endswith("fq.gz"):
        return float(_count_fastq(p))
    elif file.endswith(".fq"):
        return float(_count_fastq(p))
    elif file.endswith("fa.gz"):
        return float(_count_fasta(p))
    elif file.endswith(".fa"):
        return float(_count_fasta(p))
    return 0.0

def diamond_unassemble_data(args):
    check_read_type(args.raw_reads)
    paf1 = Paf(args.paf1)
    if args.paf2:
        paf2 = Paf(args.paf2)
    totalreadnumber = get_count_reads(args.raw_reads)
    if args.paf2:
        totalreadnumber = float(totalreadnumber) * 2
    cazyfpkm, readtable, cazy2seqid = Cal_FPKM(paf1, paf2, totalreadnumber, args.normalized)
    FPKMToCsv(args, "Diamond", cazyfpkm, readtable, cazy2seqid)

def diamond_filter(args):
    print_seqids = {}
    for line in open(args.paf1):
        lines = line.split()
        if lines[0] not in print_seqids:
            print(line.rstrip("\n"))
            print_seqids[lines[0]] = 1

def getSeqlen(paf1, paf2):
    x = paf1.GetSeqLen()
    y = paf2.GetSeqLen()
    return merge_two_dicts(x, y)

def getCazySeqId(paf1, paf2):
    cazy2seqid = {}
    paf1.CAZy2SeqID(cazy2seqid)
    paf2.CAZy2SeqID(cazy2seqid)
    for cazy in cazy2seqid:
        cazy2seqid[cazy] = set(cazy2seqid[cazy])
    return cazy2seqid

def get_subfam2seqid(paf1, paf2):
    subfam2seqid = {}
    paf1.Get_subfam2SeqID(subfam2seqid)
    paf2.Get_subfam2SeqID(subfam2seqid)
    for subfam in subfam2seqid:
        subfam2seqid[subfam] = set(subfam2seqid[subfam])
    return subfam2seqid

def getSeqReadID(paf1, paf2):
    seqid2readid = {}
    paf1.SeqID2ReadID(seqid2readid)
    paf2.SeqID2ReadID(seqid2readid)
    return seqid2readid

def SeqReadCount(seqid2readid):
    return {seqid: len(seqid2readid[seqid]) for seqid in seqid2readid}

def merge_two_dicts(x, y):
    z = x.copy()
    z.update(y)
    return z

def SequenceFPKM(readtable, seq2len, totalreadnumber):
    seqfpkm = {}
    for seqid in readtable:
        tmp_total_read = float(totalreadnumber) / pow(10, 6)
        tmp_trans_len  = float(seq2len[seqid]) / 1000
        read_count = float(readtable[seqid])
        tmp_fpkm = read_count / tmp_total_read / tmp_trans_len
        seqfpkm[seqid] = tmp_fpkm
    return seqfpkm

def SequenceTPM(readtable, seq2len, totalreadnumber):
    seqtpm = {}
    normalized_tpm = 0.0
    for seqid in readtable:
        read_count = float(readtable[seqid])
        seqlen = float(seq2len[seqid])
        normalized_tpm += read_count / seqlen
    for seqid in readtable:
        read_count = float(readtable[seqid])
        seqlen = float(seq2len[seqid])
        normalized_reads_counts = read_count / seqlen * pow(10, 6)
        tmp_seqtpm = normalized_reads_counts / normalized_tpm if normalized_tpm else 0.0
        seqtpm[seqid] = tmp_seqtpm
    return seqtpm

def SequenceRPM(readtable, seq2len, totalreadnumber):
    seqrpm = {}
    for seqid in readtable:
        read_count = float(readtable[seqid])
        rpm = read_count * pow(10, 6) / totalreadnumber if totalreadnumber else 0.0
        seqrpm[seqid] = rpm
    return seqrpm

def CAZyFPKM(seqfpkm, cazy2seqid):
    cazyfpkm = {}
    for cazy in cazy2seqid:
        tmp_fpkm = 0.0
        for seqid in cazy2seqid[cazy]:
            tmp_fpkm += float(seqfpkm[seqid])
        cazyfpkm[cazy] = tmp_fpkm
    return cazyfpkm

def Cal_FPKM(paf1, paf2, totalreadnumber, normalized):
    seq2len = getSeqlen(paf1, paf2)
    cazy2seqid = getCazySeqId(paf1, paf2)
    seqid2readid = getSeqReadID(paf1, paf2)
    readtable = SeqReadCount(seqid2readid)
    if normalized == "FPKM":
        seqfpkm = SequenceFPKM(readtable, seq2len, totalreadnumber)
    elif normalized == "RPM":
        seqfpkm = SequenceRPM(readtable, seq2len, totalreadnumber)
    else:
        seqfpkm = SequenceTPM(readtable, seq2len, totalreadnumber)
    cazyfpkm = CAZyFPKM(seqfpkm, cazy2seqid)
    return cazyfpkm, readtable, cazy2seqid

def read_EC2substrate_table(args):
    famEC2substrate = {}
    db_path = Path(args.db)
    map_table = db_path / "fam-substrate-mapping.tsv"
    map_table_lines = map_table.read_text().splitlines()
    for line in map_table_lines[1:]:
        lines = line.rstrip("\n").split("\t")
        substrates = [sub_tmp.strip(" ") for sub_tmp in lines[0].strip().replace("and", "").split(',')]
        famEC2substrate.setdefault(lines[-1], []).extend(substrates)
    for fam in famEC2substrate:
        famEC2substrate[fam] = list(set(famEC2substrate[fam]))
    return famEC2substrate

def read_CAZyID2subfam_table(args):
    CAZyID2subfam = {}
    db_path = Path(args.db)
    map_table = db_path / "CAZyID_subfam_mapping.tsv"
    map_table_lines = map_table.read_text().splitlines()
    for line in map_table_lines:
        lines = line.rstrip("\n").split("\t")
        CAZyID2subfam.setdefault(lines[-1], []).append(lines[0])
    return CAZyID2subfam

def read_subfam2ECosub_table(args):
    subfam2EC = {}
    subfam2subtrate = {}
    db_path = Path(args.db)
    map_table = db_path / "subfam_EC_mapping.tsv"
    map_table_lines = map_table.read_text().splitlines()
    for line in map_table_lines:
        lines = line.rstrip("\n").split("\t")
        if lines[-1] != "-":
            substrates = [sub.strip() for sub in lines[-1].strip().replace("and", "").split(",")]
            subfam2subtrate.setdefault(lines[0], []).extend(substrates)
        if lines[1] != "-":
            subfam2EC.setdefault(lines[0], []).append(lines[1])
    for subfam in subfam2EC:
        subfam2EC[subfam] = list(set(subfam2EC[subfam]))
    for subfam in subfam2subtrate:
        subfam2subtrate[subfam] = list(set(subfam2subtrate[subfam]))
    return subfam2EC, subfam2subtrate

def diamond_EC_abund(args):
    if not args.db.endswith("/"):
        args.db += "/"
    subfam2EC, subfam2subtrate = read_subfam2ECosub_table(args)
    EC2Abund = {}
    EC2subfam = {}
    for line in open(args.input):
        subfam, FPKM, ReadCount, SeqNum = line.rstrip("\n").split("\t")
        if subfam in subfam2EC:
            ECs = subfam2EC[subfam]
            for EC in ECs:
                subfams = EC2subfam.get(EC, [])
                if subfam not in subfams:
                    EC2subfam.setdefault(EC, []).append(subfam)
                    EC2Abund.setdefault(EC, []).append(float(FPKM))
    outfilename = args.output
    with open(outfilename, 'w') as f:
        f.write("EC\tAbundance\tsubfam\n")
        for sub in EC2Abund:
            f.write(sub + "\t" + str(sum(EC2Abund[sub])) + "\t" + ";".join(EC2subfam[sub]) + "\n")

def CAZyme_substrate(args):
    if not args.db.endswith("/"):
        args.db += "/"
    EC2substrate = read_EC2substrate_table(args)
    subfam2EC, subfam2subtrate = read_subfam2ECosub_table(args)
    Sub2Abund = {}
    Sub2subfam = {}
    for line in open(args.input):
        subfam, FPKM, SeqNum, ReadCount = line.rstrip("\n").split("\t")
        if subfam in subfam2EC:
            ECs = subfam2EC[subfam]
            if ECs:
                for EC in ECs:
                    substrates = EC2substrate.get(EC, "")
                    if substrates:
                        for sub in substrates:
                            subfams = Sub2subfam.get(sub, [])
                            if subfam not in subfams:
                                Sub2Abund.setdefault(sub, []).append(float(FPKM))
                                Sub2subfam.setdefault(sub, []).append(subfam)
        substrates = subfam2subtrate.get(subfam, "")
        if substrates:
            for sub in substrates:
                subfams = Sub2subfam.get(sub, [])
                if subfam not in subfams:
                    Sub2Abund.setdefault(sub, []).append(float(FPKM))
                    Sub2subfam.setdefault(sub, []).append(subfam)
    outfilename = args.output
    with open(outfilename, 'w') as f:
        f.write("Substrate\tAbundance\tsubfam\n")
        for sub in Sub2Abund:
            f.write(sub + "\t" + str(sum(Sub2Abund[sub])) + "\t" + ";".join(Sub2subfam[sub]) + "\n")

def Cal_subfam_FPKM(paf1, paf2, totalreadnumber, normalized):
    seq2len = getSeqlen(paf1, paf2)
    subfam2seqid = get_subfam2seqid(paf1, paf2)
    seqid2readid = getSeqReadID(paf1, paf2)
    readtable = SeqReadCount(seqid2readid)
    if normalized == "FPKM":
        seqfpkm = SequenceFPKM(readtable, seq2len, totalreadnumber)
    elif normalized == "RPM":
        seqfpkm = SequenceRPM(readtable, seq2len, totalreadnumber)
    else:
        seqfpkm = SequenceTPM(readtable, seq2len, totalreadnumber)
    cazyfpkm = CAZyFPKM(seqfpkm, subfam2seqid)
    return cazyfpkm, readtable, subfam2seqid

def diamond_subfam_abund(args):
    if not args.db.endswith("/"):
        args.db += "/"
    check_read_type(args.raw_reads)
    CAZyID2subfam = read_CAZyID2subfam_table(args)
    paf1 = Paf(args.paf1)
    if args.paf2:
        paf2 = Paf(args.paf2)
    paf1.Assign_subfam(CAZyID2subfam)
    paf2.Assign_subfam(CAZyID2subfam)
    totalreadnumber = get_count_reads(args.raw_reads)
    if args.paf2:
        totalreadnumber = float(totalreadnumber) * 2
    subfamfpkm, readtable, subfam2seqid = Cal_subfam_FPKM(paf1, paf2, totalreadnumber, args.normalized)
    FPKMToCsv(args, "Diamond", subfamfpkm, readtable, subfam2seqid)

# -------------------- CLICK / RICH-CLICK CLI --------------------

GROUP_HELP = """
[bold]Assembly-free CAZyme abundance utilities (DIAMOND)[/bold]

[bold]Examples[/bold]:

  1) Family abundance (paired):
     dbcan_asmfree diamond_fam_abund \\
       -paf1 Dry2014_1.blastx -paf2 Dry2014_2.blastx \\
       --raw_reads Dry2014_1_val_1.fq.gz -n FPKM -o Dry2014_fam_abund

  2) Family abundance (paired):
     dbcan_asmfree diamond_fam_abund \\
       -paf1 Wet2014_1.blastx -paf2 Wet2014_2.blastx \\
       --raw_reads Wet2014_1_val_1.fq.gz -n FPKM -o Wet2014_fam_abund

  3) Sub-family abundance:
     dbcan_asmfree diamond_subfam_abund \\
       -paf1 Dry2014_1.blastx -paf2 Dry2014_2.blastx \\
       --raw_reads Dry2014_1_val_1.fq.gz -n FPKM -o Dry2014_subfam_abund

  4) EC abundance:
     dbcan_asmfree diamond_EC_abund \\
       -i Dry2014_subfam_abund -o Dry2014_EC_abund

  5) Substrate abundance:
     dbcan_asmfree diamond_substrate_abund \\
       -i Dry2014_subfam_abund -o Dry2014_substrate_abund
"""

@click.group(help=GROUP_HELP)
def cli():
    pass

def common_norm_option(f):
    return click.option(
        "-n", "--normalized",
        type=click.Choice(["FPKM", "TPM", "RPM"], case_sensitive=False),
        default="TPM",
        show_default=True,
        help="Normalization method"
    )(f)

@cli.command("diamond_fam_abund", help="Compute CAZy family abundance (FPKM/TPM/RPM).")
@click.option("-paf1", required=True, type=click.Path(exists=True), help="R1 DIAMOND blastx output.")
@click.option("-paf2", type=click.Path(exists=True), default="", help="R2 DIAMOND blastx output (optional).")
@click.option("--raw_reads", required=True, type=click.Path(exists=True), help="Raw reads file (fq/fa[.gz]).")
@click.option("-d", "--db", default="./db", show_default=True, type=click.Path(), help="Database directory.")
@click.option("-o", "--output", default="asmfree_fam_abund", show_default=True, help="Output file.")
@common_norm_option
def cmd_fam_abund(paf1, paf2, raw_reads, db, output, normalized):
    class Args: pass
    args = Args()
    args.paf1 = paf1
    args.paf2 = paf2 if paf2 else None
    args.raw_reads = raw_reads
    args.db = db
    args.output = output
    args.normalized = normalized.upper()
    rprint("[bold green]Running family abundance...[/bold green]")
    diamond_unassemble_data(args)
    rprint(f"[bold cyan]Done -> {output}[/bold cyan]")

@cli.command("diamond_subfam_abund", help="Compute CAZy sub-family abundance.")
@click.option("-paf1", required=True, type=click.Path(exists=True))
@click.option("-paf2", type=click.Path(exists=True), default="", help="R2 DIAMOND blastx output (optional).")
@click.option("--raw_reads", required=True, type=click.Path(exists=True))
@click.option("-d", "--db", default="./db", show_default=True, type=click.Path())
@click.option("-o", "--output", default="asmfree_subfam_abund", show_default=True)
@common_norm_option
def cmd_subfam_abund(paf1, paf2, raw_reads, db, output, normalized):
    class Args: pass
    args = Args()
    args.paf1 = paf1
    args.paf2 = paf2 if paf2 else None
    args.raw_reads = raw_reads
    args.db = db
    args.output = output
    args.normalized = normalized.upper()
    rprint("[bold green]Running sub-family abundance...[/bold green]")
    diamond_subfam_abund(args)
    rprint(f"[bold cyan]Done -> {output}[/bold cyan]")

@cli.command("diamond_EC_abund", help="Summarize EC abundance from sub-family abundance file.")
@click.option("-i", "--input", required=True, type=click.Path(exists=True), help="Sub-family abundance input file.")
@click.option("-d", "--db", default="./db", show_default=True, type=click.Path())
@click.option("-o", "--output", default="EC_abund.tsv", show_default=True)
def cmd_ec_abund(input, db, output):
    class Args: pass
    args = Args()
    args.input = input
    args.db = db
    args.output = output
    rprint("[bold green]Running EC abundance...[/bold green]")
    diamond_EC_abund(args)
    rprint(f"[bold cyan]Done -> {output}[/bold cyan]")

@cli.command("diamond_substrate_abund", help="Infer substrate abundance from sub-family abundance file.")
@click.option("-i", "--input", required=True, type=click.Path(exists=True), help="Sub-family abundance input file.")
@click.option("-d", "--db", default="./db", show_default=True, type=click.Path())
@click.option("-o", "--output", default="substrate_abund.tsv", show_default=True)
def cmd_substrate_abund(input, db, output):
    class Args: pass
    args = Args()
    args.input = input
    args.db = db
    args.output = output
    rprint("[bold green]Running substrate abundance...[/bold green]")
    CAZyme_substrate(args)
    rprint(f"[bold cyan]Done -> {output}[/bold cyan]")

if __name__ == "__main__":
    cli()
