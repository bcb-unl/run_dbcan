import os
import sys
from pathlib import Path
import gzip
import logging
import gc
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
import rich_click as click
from rich import print as rprint
from dbcan.parameter import logging_options
from dbcan.main import setup_logging

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

def extract_CAZy_from_Tsn(tsn):
    """
    Extract CAZy annotations from target sequence name (Tsn).
    
    Normal case: CAD7272370.1|GH33 -> CAZy is ['GH33'] (content after |)
    Special case: CE16|390012|Krezon1_GeneCatalog_proteins_20220127.aa.fasta 
                  -> When ID contains '.fasta' (case-insensitive), if the first part after | is a pure number,
                     then CAZy is the part before | (i.e., 'CE16')
                  -> Otherwise, extract all parts before the first pure number
    
    Args:
        tsn: Target sequence name from DIAMOND result (sseqid)
    
    Returns:
        List of CAZy annotations (stripped of whitespace)
    """
    # Check if '.fasta' is in the ID (case-insensitive, more precise than just 'fasta')
    tsn_lower = tsn.lower()
    if '.fasta' in tsn_lower:
        # Split by |
        parts = tsn.split('|')
        if len(parts) < 2:
            return []
        
        # Check if the first part after | is a pure number
        first_part_after_pipe = parts[1].strip()
        if first_part_after_pipe.isdigit():
            # Special case: first part after | is a pure number
            # Return the part before | as CAZy annotation
            return [parts[0].strip()]
        
        # Otherwise, collect all parts before the first pure number
        cazy_parts = []
        for part in parts[1:]:
            stripped_part = part.strip()
            if stripped_part.isdigit():
                # Found first pure number, stop here
                break
            else:
                # Collect non-numeric parts (stripped)
                cazy_parts.append(stripped_part)
        
        # Return collected parts, or empty list if none found
        return cazy_parts
    else:
        # Normal case: return everything after the first |, with each part stripped
        parts = tsn.strip("|").split("|")[1:]
        return [part.strip() for part in parts]

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
        # Use special extraction function for CAZy annotations
        cazy_list = extract_CAZy_from_Tsn(self.Tsn)
        self.CAZys = CAZy_filter(cazy_list) if cazy_list else set()
        self.UniReadId = lines[0].split(".")[0]

    def __str__(self):
        return "\t".join([str(getattr(self, value)) for value in vars(self) if value != "CAZys"])

class Paf(object):
    """PAF file processor using streaming iterator pattern to avoid loading all records into memory"""
    def __init__(self, filename):
        self.filename = filename
        # No longer load all records at once, use streaming processing instead
    
    def __iter__(self):
        """Generator mode, read and parse line by line to avoid memory overflow"""
        with open(self.filename) as f:
            for line in f:
                line = line.strip()
                if line:  # Skip empty lines
                    try:
                        yield PafRecord(line.split())
                    except (IndexError, ValueError) as e:
                        # Skip malformed lines
                        logging.warning(f"Skipping malformed line in {self.filename}: {e}")
                        continue
    
    def GetReadId(self):
        """Get all Read ID list (maintain compatibility, but use generator)"""
        return [record.Qsn for record in self]
    
    def GetSeqId(self):
        """Get all Seq ID list (maintain compatibility, but use generator)"""
        return [record.SeqID for record in self]
    
    def GetSeqLen(self):
        """Get sequence ID to length mapping (single traversal)"""
        seq2len = {}
        for record in self:
            # Keep only the maximum length for each sequence ID (if duplicates exist)
            seqid = record.SeqID
            if seqid not in seq2len or int(record.Tsl) > int(seq2len[seqid]):
                seq2len[seqid] = record.Tsl
        return seq2len
    
    def CAZy2SeqID(self, CazySeqId):
        """Map CAZy ID to sequence ID.
        
        Accepts any dict-like object and uses setdefault for compatibility.
        """
        for record in self:
            for cazy in record.CAZys:
                CazySeqId.setdefault(cazy, []).append(record.SeqID)
    
    def SeqID2ReadID(self, aa):
        """Map sequence ID to Read ID.
        
        Accepts any dict-like object and uses setdefault for compatibility.
        """
        for record in self:
            aa.setdefault(record.SeqID, []).append(record.Qsn)
    
    def ReadID2Record(self):
        """Get Read ID to record mapping (Note: this will load all records into memory)"""
        # This method requires all records, but kept for compatibility
        # For large files, it is recommended to avoid using this method
        result = {}
        for record in self:
            result[record.Qsn] = record
        return result
    
    def Output(self):
        """Output all records"""
        for record in self:
            print(record)
    
    def Assign_CAZy_megahit(self):
        """Assign CAZy for megahit format (streaming processing)"""
        # Note: This method modifies records, but due to streaming processing, modifications are not persistent
        # If persistence is needed, redesign is required
        for record in self:
            record.CAZys = CAZy_filter(record.Qsn.strip("|").split("|")[1:])
    
    def Assign_subfam(self, CAZyID2subfam):
        """Assign subfamily mapping.
        
        .. note:: **API Change (Breaking)**: This method no longer modifies individual 
           record objects. Due to streaming processing, records are generated on-the-fly 
           and cannot be modified persistently. Instead, this method stores the mapping 
           as an instance variable, which is used by Get_subfam2SeqID() during traversal.
        
        Previous behavior: Modified each record's subfams attribute directly.
        Current behavior: Stores mapping for later use during iteration.
        
        Args:
            CAZyID2subfam: Dictionary mapping CAZy IDs to subfamily lists
        """
        self.CAZyID2subfam = CAZyID2subfam
    
    def Get_subfam2SeqID(self, subfam2SeqID):
        """Get subfamily to sequence ID mapping.
        
        Dynamically looks up subfams during traversal using the mapping stored by 
        Assign_subfam(). This is necessary because records are generated on-the-fly 
        in streaming mode and cannot have persistent attributes.
        
        Args:
            subfam2SeqID: Dictionary to populate (accepts any dict-like object)
        """
        CAZyID2subfam = getattr(self, 'CAZyID2subfam', {})
        for record in self:
            # Dynamically look up subfam from stored mapping
            subfams = CAZyID2subfam.get(record.Tsn, [])
            for subfam in subfams:
                subfam2SeqID.setdefault(subfam, []).append(record.SeqID)

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
    """Process unassembled DIAMOND data (supports multithreading)"""
    check_read_type(args.raw_reads)
    threads = getattr(args, 'threads', 1)
    
    paf1 = Paf(args.paf1)
    paf2 = Paf(args.paf2) if args.paf2 else None
    totalreadnumber = get_count_reads(args.raw_reads)
    if args.paf2:
        totalreadnumber = float(totalreadnumber) * 2
    
    # Process PAF files
    logging.info("Processing PAF files...")
    cazyfpkm, readtable, cazy2seqid = Cal_FPKM(paf1, paf2, totalreadnumber, args.normalized, threads)
    logging.info("PAF file processing completed")
    
    FPKMToCsv(args, "Diamond", cazyfpkm, readtable, cazy2seqid)
    
    # Free memory
    del paf1, paf2, cazyfpkm, readtable, cazy2seqid
    gc.collect()

def diamond_filter(args):
    """Filter duplicate sequence IDs (use set to optimize memory)"""
    print_seqids = set()  # Use set instead of dict to save memory
    for line in open(args.paf1):
        lines = line.split()
        if lines[0] not in print_seqids:
            print(line.rstrip("\n"))
            print_seqids.add(lines[0])

def getSeqlen(paf1, paf2):
    """Get sequence length mapping (merge two PAF files)"""
    x = paf1.GetSeqLen()
    y = paf2.GetSeqLen() if paf2 else {}
    return merge_two_dicts(x, y)

def getCazySeqId(paf1, paf2):
    """Get CAZy to sequence ID mapping (optimized with defaultdict and set)"""
    cazy2seqid = defaultdict(set)  # Use set directly to avoid subsequent conversion
    for record in paf1:
        for cazy in record.CAZys:
            cazy2seqid[cazy].add(record.SeqID)
    if paf2:
        for record in paf2:
            for cazy in record.CAZys:
                cazy2seqid[cazy].add(record.SeqID)
    return dict(cazy2seqid)  # Convert to regular dict to maintain compatibility

def get_subfam2seqid(paf1, paf2):
    """Get subfamily to sequence ID mapping (optimized with defaultdict and set)"""
    subfam2seqid = defaultdict(set)  # Use set directly to avoid subsequent conversion
    paf1.Get_subfam2SeqID(subfam2seqid)
    if paf2:
        paf2.Get_subfam2SeqID(subfam2seqid)
    return dict(subfam2seqid)  # Convert to regular dict to maintain compatibility

def getSeqReadID(paf1, paf2):
    """Get sequence ID to Read ID mapping (optimized with defaultdict)"""
    seqid2readid = defaultdict(list)  # Use defaultdict to avoid setdefault calls
    for record in paf1:
        seqid2readid[record.SeqID].append(record.Qsn)
    if paf2:
        for record in paf2:
            seqid2readid[record.SeqID].append(record.Qsn)
    return dict(seqid2readid)  # Convert to regular dict to maintain compatibility

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

def _process_paf_chunk(filename, start_line, chunk_size, CAZyID2subfam=None):
    """Process a chunk of PAF file (for multithreading)"""
    seq2len = {}
    cazy2seqid = defaultdict(set)
    seqid2readid = defaultdict(list)
    subfam2seqid = defaultdict(set) if CAZyID2subfam else None
    
    try:
        with open(filename) as f:
            # Skip preceding lines
            for _ in range(start_line):
                try:
                    next(f)
                except StopIteration:
                    break
            
            # Process current chunk
            for i, line in enumerate(f):
                if i >= chunk_size:
                    break
                line = line.strip()
                if not line:
                    continue
                try:
                    parts = line.split()
                    if len(parts) < 14:  # Ensure sufficient columns
                        continue
                    record = PafRecord(parts)
                    # Collect sequence length
                    seqid = record.SeqID
                    tsl_int = int(record.Tsl)
                    if seqid not in seq2len or tsl_int > int(seq2len[seqid]):
                        seq2len[seqid] = record.Tsl
                    
                    # Collect CAZy to sequence ID mapping
                    for cazy in record.CAZys:
                        cazy2seqid[cazy].add(seqid)
                    
                    # Collect sequence ID to Read ID mapping
                    seqid2readid[seqid].append(record.Qsn)
                    
                    # If subfam mapping exists, collect subfam information
                    if CAZyID2subfam and subfam2seqid is not None:
                        subfams = CAZyID2subfam.get(record.Tsn, [])
                        for subfam in subfams:
                            subfam2seqid[subfam].add(seqid)
                except (IndexError, ValueError) as e:
                    logging.debug(f"Skipping malformed line in chunk: {e}")
                    continue
    except Exception as e:
        logging.error(f"Error processing chunk starting at line {start_line}: {e}")
        raise
    
    result = {
        'seq2len': dict(seq2len),
        'cazy2seqid': dict(cazy2seqid),
        'seqid2readid': dict(seqid2readid)
    }
    if subfam2seqid is not None:
        result['subfam2seqid'] = dict(subfam2seqid)
    
    return result

def _count_lines(filename):
    """Quickly count file lines"""
    with open(filename) as f:
        return sum(1 for _ in f)

def Cal_FPKM(paf1, paf2, totalreadnumber, normalized, threads=1):
    """Calculate FPKM/TPM/RPM (supports multithreading)"""
    # If thread count is 1 or file is small, use single-threaded processing
    use_threading = threads > 1
    
    if use_threading:
        # Multithreaded processing
        chunk_size = 10000  # Process 10000 lines per chunk
        
        def process_file(filename, CAZyID2subfam=None):
            total_lines = _count_lines(filename)
            num_chunks = max(1, (total_lines + chunk_size - 1) // chunk_size)
            
            # Use thread pool for processing
            results = []
            with ThreadPoolExecutor(max_workers=threads) as executor:
                futures = []
                for i in range(num_chunks):
                    start_line = i * chunk_size
                    future = executor.submit(_process_paf_chunk, filename, start_line, chunk_size, CAZyID2subfam)
                    futures.append(future)
                
                for future in as_completed(futures):
                    results.append(future.result())
            
            # Merge results
            seq2len = {}
            cazy2seqid = defaultdict(set)
            seqid2readid = defaultdict(list)
            subfam2seqid = defaultdict(set)
            
            for result in results:
                # Merge seq2len (keep maximum length)
                for seqid, length in result['seq2len'].items():
                    if seqid not in seq2len or int(length) > int(seq2len[seqid]):
                        seq2len[seqid] = length
                
                # Merge cazy2seqid
                for cazy, seqids in result['cazy2seqid'].items():
                    cazy2seqid[cazy].update(seqids)
                
                # Merge seqid2readid
                for seqid, readids in result['seqid2readid'].items():
                    seqid2readid[seqid].extend(readids)
                
                # Merge subfam2seqid
                if 'subfam2seqid' in result and result['subfam2seqid']:
                    for subfam, seqids in result['subfam2seqid'].items():
                        subfam2seqid[subfam].update(seqids)
            
            subfam_result = dict(subfam2seqid) if subfam2seqid else None
            return dict(seq2len), dict(cazy2seqid), dict(seqid2readid), subfam_result
        
        seq2len1, cazy2seqid1, seqid2readid1, _ = process_file(paf1.filename)
        if paf2:
            seq2len2, cazy2seqid2, seqid2readid2, _ = process_file(paf2.filename)
            # Merge results from two files
            seq2len = merge_two_dicts(seq2len1, seq2len2)
            # For seq2len, keep maximum length
            for seqid, length in seq2len2.items():
                if seqid not in seq2len1 or int(length) > int(seq2len1[seqid]):
                    seq2len[seqid] = length
            
            cazy2seqid = defaultdict(set)
            # Merge cazy2seqid
            for cazy in set(cazy2seqid1.keys()) | set(cazy2seqid2.keys()):
                seqids1 = cazy2seqid1.get(cazy, set())
                seqids2 = cazy2seqid2.get(cazy, set())
                if isinstance(seqids1, list):
                    seqids1 = set(seqids1)
                if isinstance(seqids2, list):
                    seqids2 = set(seqids2)
                cazy2seqid[cazy] = seqids1 | seqids2
            
            seqid2readid = defaultdict(list)
            # Merge seqid2readid
            for seqid in set(seqid2readid1.keys()) | set(seqid2readid2.keys()):
                seqid2readid[seqid] = seqid2readid1.get(seqid, []) + seqid2readid2.get(seqid, [])
            
            cazy2seqid = dict(cazy2seqid)
            seqid2readid = dict(seqid2readid)
        else:
            seq2len = seq2len1
            cazy2seqid = cazy2seqid1
            seqid2readid = seqid2readid1
    else:
        # Single-threaded processing (original logic)
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
    
    # Free memory
    del seq2len, seqid2readid
    gc.collect()
    
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

def Cal_subfam_FPKM(paf1, paf2, totalreadnumber, normalized, threads=1):
    """Calculate subfamily FPKM/TPM/RPM (supports multithreading)"""
    # Get CAZyID2subfam mapping
    CAZyID2subfam = getattr(paf1, 'CAZyID2subfam', None)
    if not CAZyID2subfam and paf2:
        CAZyID2subfam = getattr(paf2, 'CAZyID2subfam', None)
    
    # If thread count is 1 or file is small, use single-threaded processing
    use_threading = threads > 1 and CAZyID2subfam
    
    if use_threading:
        # Multithreaded processing
        chunk_size = 10000
        
        def process_file(filename, CAZyID2subfam):
            total_lines = _count_lines(filename)
            num_chunks = max(1, (total_lines + chunk_size - 1) // chunk_size)
            
            results = []
            with ThreadPoolExecutor(max_workers=threads) as executor:
                futures = []
                for i in range(num_chunks):
                    start_line = i * chunk_size
                    future = executor.submit(_process_paf_chunk, filename, start_line, chunk_size, CAZyID2subfam)
                    futures.append(future)
                
                for future in as_completed(futures):
                    results.append(future.result())
            
            # Merge results
            seq2len = {}
            subfam2seqid = defaultdict(set)
            seqid2readid = defaultdict(list)
            
            for result in results:
                for seqid, length in result['seq2len'].items():
                    if seqid not in seq2len or int(length) > int(seq2len[seqid]):
                        seq2len[seqid] = length
                
                for subfam, seqids in result['subfam2seqid'].items():
                    subfam2seqid[subfam].update(seqids)
                
                for seqid, readids in result['seqid2readid'].items():
                    seqid2readid[seqid].extend(readids)
            
            return dict(seq2len), dict(subfam2seqid), dict(seqid2readid)
        
        seq2len1, subfam2seqid1, seqid2readid1 = process_file(paf1.filename, CAZyID2subfam)
        if paf2:
            seq2len2, subfam2seqid2, seqid2readid2 = process_file(paf2.filename, CAZyID2subfam)
            seq2len = merge_two_dicts(seq2len1, seq2len2)
            subfam2seqid = defaultdict(set)
            for subfam in set(subfam2seqid1.keys()) | set(subfam2seqid2.keys()):
                subfam2seqid[subfam] = subfam2seqid1.get(subfam, set()) | subfam2seqid2.get(subfam, set())
            seqid2readid = defaultdict(list)
            for seqid in set(seqid2readid1.keys()) | set(seqid2readid2.keys()):
                seqid2readid[seqid] = seqid2readid1.get(seqid, []) + seqid2readid2.get(seqid, [])
            subfam2seqid = dict(subfam2seqid)
            seqid2readid = dict(seqid2readid)
        else:
            seq2len = seq2len1
            subfam2seqid = subfam2seqid1
            seqid2readid = seqid2readid1
    else:
        # Single-threaded processing (original logic)
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
    
    # Free memory
    del seq2len, seqid2readid
    gc.collect()
    
    return cazyfpkm, readtable, subfam2seqid

def diamond_subfam_abund(args):
    """Calculate subfamily abundance (supports multithreading)"""
    if not args.db.endswith("/"):
        args.db += "/"
    check_read_type(args.raw_reads)
    threads = getattr(args, 'threads', 1)
    
    CAZyID2subfam = read_CAZyID2subfam_table(args)
    paf1 = Paf(args.paf1)
    paf2 = Paf(args.paf2) if args.paf2 else None
    paf1.Assign_subfam(CAZyID2subfam)
    if paf2:
        paf2.Assign_subfam(CAZyID2subfam)
    
    totalreadnumber = get_count_reads(args.raw_reads)
    if args.paf2:
        totalreadnumber = float(totalreadnumber) * 2
    
    # Process subfamily abundance
    logging.info("Processing subfamily abundance...")
    subfamfpkm, readtable, subfam2seqid = Cal_subfam_FPKM(paf1, paf2, totalreadnumber, args.normalized, threads)
    logging.info("Subfamily abundance processing completed")
    
    FPKMToCsv(args, "Diamond", subfamfpkm, readtable, subfam2seqid)
    
    # Free memory
    del paf1, paf2, CAZyID2subfam, subfamfpkm, readtable, subfam2seqid
    gc.collect()

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
@logging_options
@click.pass_context
def cli(ctx, log_level, log_file, verbose):
    setup_logging(log_level, log_file, verbose)
    ctx.ensure_object(dict)
    ctx.obj['log_level'] = log_level
    ctx.obj['log_file'] = log_file
    ctx.obj['verbose'] = verbose

def common_norm_option(f):
    return click.option(
        "-n", "--normalized",
        type=click.Choice(["FPKM", "TPM", "RPM"], case_sensitive=False),
        default="TPM",
        show_default=True,
        help="Normalization method"
    )(f)

def common_threads_option(f):
    return click.option(
        "-t", "--threads",
        type=int,
        default=min(os.cpu_count() or 1, 4),
        show_default=True,
        help="Number of threads for parallel processing (default: min(cpu_count, 4) for I/O-bound tasks)"
    )(f)

@cli.command("diamond_fam_abund", help="Compute CAZy family abundance (FPKM/TPM/RPM).")
@click.option("-paf1", required=True, type=click.Path(exists=True), help="R1 DIAMOND blastx output.")
@click.option("-paf2", type=click.Path(exists=True), default="", help="R2 DIAMOND blastx output (optional).")
@click.option("--raw_reads", required=True, type=click.Path(exists=True), help="Raw reads file (fq/fa[.gz]).")
@click.option("-d", "--db", default="./db", show_default=True, type=click.Path(), help="Database directory.")
@click.option("-o", "--output", default="asmfree_fam_abund", show_default=True, help="Output file.")
@common_norm_option
@common_threads_option
def cmd_fam_abund(paf1, paf2, raw_reads, db, output, normalized, threads):
    class Args: pass
    args = Args()
    args.paf1 = paf1
    args.paf2 = paf2 if paf2 else None
    args.raw_reads = raw_reads
    args.db = db
    args.output = output
    args.normalized = normalized.upper()
    args.threads = threads
    rprint(f"[bold green]Running family abundance (threads={threads})...[/bold green]")
    diamond_unassemble_data(args)
    rprint(f"[bold cyan]Done -> {output}[/bold cyan]")

@cli.command("diamond_subfam_abund", help="Compute CAZy sub-family abundance.")
@click.option("-paf1", required=True, type=click.Path(exists=True))
@click.option("-paf2", type=click.Path(exists=True), default="", help="R2 DIAMOND blastx output (optional).")
@click.option("--raw_reads", required=True, type=click.Path(exists=True))
@click.option("-d", "--db", default="./db", show_default=True, type=click.Path())
@click.option("-o", "--output", default="asmfree_subfam_abund", show_default=True)
@common_norm_option
@common_threads_option
def cmd_subfam_abund(paf1, paf2, raw_reads, db, output, normalized, threads):
    class Args: pass
    args = Args()
    args.paf1 = paf1
    args.paf2 = paf2 if paf2 else None
    args.raw_reads = raw_reads
    args.db = db
    args.output = output
    args.normalized = normalized.upper()
    args.threads = threads
    rprint(f"[bold green]Running sub-family abundance (threads={threads})...[/bold green]")
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
