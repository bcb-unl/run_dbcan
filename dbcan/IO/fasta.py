from pathlib import Path
import os
import logging
from multiprocessing.pool import ThreadPool
import pyrodigal
import contextlib
import gzip
from Bio import SeqIO
from Bio.Seq import Seq
import re
import shutil
import tempfile
from io import StringIO

from dbcan.configs.base_config import GeneralConfig
from dbcan.constants.base_constants import (
    INPUT_PROTEIN_NAME, INPUT_PRODIGAL_GFF_NAME,
    EUK_WARN_MIN_GENES, EUK_WARN_MEDIAN_AA, EUK_WARN_SHORT_AA, EUK_WARN_SHORT_FRACTION,
    PROTEIN_STOP_CHARS, PROTEIN_UNKNOWN_AA,
)

logger = logging.getLogger(__name__)

class BaseProcessor:
    """Base processor for sequence data."""

    def __init__(self, config: GeneralConfig):
        self.config = config
        self._setup_processor()

    # Paths as properties (single source: config)
    @property
    def input_raw_data_path(self) -> Path:
        return Path(self.config.input_raw_data)

    @property
    def output_dir_path(self) -> Path:
        return Path(self.config.output_dir)

    # Prefer using Path objects; call str(path) only at I/O boundaries if needed.
    def _setup_processor(self):
        # Ensure output directory exists
        self.output_dir_path.mkdir(parents=True, exist_ok=True)
        # Hook for subclasses
        self._additional_setup()
        # Validate inputs
        self._validate_attributes()

    def _additional_setup(self):
        """Subclass hook."""
        pass

    def _validate_attributes(self):
        if not self.input_raw_data_path.exists():
            raise FileNotFoundError(f"Input file not found: {str(self.input_raw_data_path)}")
        if self.input_raw_data_path.stat().st_size == 0:
            raise ValueError(f"Input file is empty: {str(self.input_raw_data_path)}")

    def parse(self, path) -> tuple:
        """Generic FASTA parser that handles gzipped files."""
        p = Path(path)

        def zopen(f: Path, mode="r"):
            return gzip.open(f, mode) if f.suffix == ".gz" else open(f, mode)

        try:
            with contextlib.ExitStack() as ctx:
                file = ctx.enter_context(zopen(p, "rt"))
                id_, desc, seq = None, None, []
                for line in file:
                    if line.startswith(">"):
                        if id_ is not None:
                            yield (id_, "".join(seq), desc)
                        fields = line[1:].strip().split(maxsplit=1)
                        id_ = fields[0] if fields else ""
                        desc = fields[1] if len(fields) > 1 else ""
                        seq = []
                    else:
                        seq.append(line.strip())
                if id_ is not None:
                    yield (id_, "".join(seq), desc)
        except Exception as e:
            logger.exception(f"Error parsing file {p}: {str(e)}")
            raise

    def process_input(self):
        raise NotImplementedError("Subclasses must implement process_input method")

    def _verify_output(self, output_path, expected_content: bool = True) -> bool:
        p = Path(output_path)
        if not p.exists():
            logger.error(f"Output file was not created: {p}")
            return False
        if expected_content and p.stat().st_size == 0:
            logger.error(f"Output file is empty: {p}")
            return False
        return True


class ProkProcessor(BaseProcessor):
    """Processor for prokaryotic genomic data."""

    @property
    def output_faa_path(self) -> Path:
        return self.output_dir_path / INPUT_PROTEIN_NAME

    @property
    def output_gff_path(self) -> Path:
        return self.output_dir_path / INPUT_PRODIGAL_GFF_NAME

    @property
    def threads(self) -> int:
        threads = getattr(self.config, 'threads', None)
        if not isinstance(threads, int) or threads is None or threads <= 0:
            threads = os.cpu_count() or 1
        return max(1, threads)

    def process_input(self):
        return self.process_fna(False)

    def process_fna(self, is_meta: bool):
        processor_type = "metagenomic" if is_meta else "prokaryotic"
        logger.info(f'Processing {processor_type} genome with Pyrodigal')

        try:
            gene_finder = pyrodigal.GeneFinder(meta=is_meta)

            # For non-meta mode, we need to collect sequences for training
            # For meta mode, we can process streamingly
            if not is_meta:
                # First pass: collect sequences for training
                logger.info(f"Reading sequences for training from {self.input_raw_data_path}")
                training_sequences = []
                sequence_count = 0
                for rec_id, seq, _ in self.parse(self.input_raw_data_path):
                    training_sequences.append(seq.encode('utf-8'))
                    sequence_count += 1
                
                if not training_sequences:
                    logger.error(f"No valid sequences found in {self.input_raw_data_path}")
                    return None, None
                
                logger.info(f"Found {sequence_count} sequences for training")
                logger.info("Training Pyrodigal on input sequences")
                gene_finder.train(*training_sequences)
                # Release training data from memory
                del training_sequences

            logger.info(f"Finding genes using {self.threads} threads")
            logger.info(f"Writing protein translations to {self.output_faa_path}")
            logger.info(f"Writing gene annotations to {self.output_gff_path}")

            # Second pass: process sequences streamingly
            genes_found = 0
            complete_lengths = []  # protein lengths (aa) of complete genes, for the eukaryote check
            sequence_count = 0
            
            with self.output_faa_path.open('w') as prot_file, self.output_gff_path.open('w') as gff_file:
                # Process sequences in batches to balance memory and performance
                batch_size = max(100, min(1000, self.threads * 50))  # Adaptive batch size
                batch = []
                batch_ids = []
                
                def process_batch(batch_data, batch_ids_data):
                    """Process a batch of sequences."""
                    if self.threads == 1:
                        return [gene_finder.find_genes(seq_bytes) for seq_bytes in batch_data]
                    else:
                        with ThreadPool(self.threads) as pool:
                            return pool.map(gene_finder.find_genes, batch_data)
                
                for rec_id, seq, _ in self.parse(self.input_raw_data_path):
                    batch.append(seq.encode('utf-8'))
                    batch_ids.append(rec_id)
                    sequence_count += 1
                    
                    if len(batch) >= batch_size:
                        # Process batch
                        results = process_batch(batch, batch_ids)
                        
                        # Write results
                        for (ori_seq_id, seq_bytes), genes in zip(zip(batch_ids, batch), results):
                            genes.write_gff(gff_file, sequence_id=ori_seq_id)
                            temp_output = StringIO()
                            genes.write_translations(temp_output, sequence_id=ori_seq_id)
                            temp_output.seek(0)
                            
                            for line in temp_output:
                                if line.startswith(">"):
                                    if " # " in line:
                                        clean_id = line[1:].strip().split(" # ")[0]
                                        prot_file.write(f">{clean_id}\n")
                                    else:
                                        clean_id = line[1:].strip().split()[0]
                                        prot_file.write(f">{clean_id}\n")
                                else:
                                    prot_file.write(line)
                            genes_found += len(genes)
                            _collect_complete_lengths(genes, complete_lengths)
                        
                        # Clear batch
                        batch = []
                        batch_ids = []
                
                # Process remaining sequences
                if batch:
                    results = process_batch(batch, batch_ids)
                    for (ori_seq_id, seq_bytes), genes in zip(zip(batch_ids, batch), results):
                        genes.write_gff(gff_file, sequence_id=ori_seq_id)
                        temp_output = StringIO()
                        genes.write_translations(temp_output, sequence_id=ori_seq_id)
                        temp_output.seek(0)
                        
                        for line in temp_output:
                            if line.startswith(">"):
                                if " # " in line:
                                    clean_id = line[1:].strip().split(" # ")[0]
                                    prot_file.write(f">{clean_id}\n")
                                else:
                                    clean_id = line[1:].strip().split()[0]
                                    prot_file.write(f">{clean_id}\n")
                            else:
                                prot_file.write(line)
                        genes_found += len(genes)
                        _collect_complete_lengths(genes, complete_lengths)

            logger.info(f"Processed {sequence_count} sequences, found {genes_found} genes")
            warn_if_eukaryotic_proteome(complete_lengths)

            if not self._verify_output(self.output_faa_path):
                return None, None
            if not self._verify_output(self.output_gff_path):
                return None, None

            return str(self.output_faa_path), str(self.output_gff_path)

        except Exception as e:
            logger.exception(f"Error processing {processor_type} genome: {str(e)}")
            return None, None


def _collect_complete_lengths(genes, out: list) -> None:
    """Append the protein length (aa, without stop) of every gene with both ends on the contig."""
    for g in genes:
        if getattr(g, "partial_begin", False) or getattr(g, "partial_end", False):
            continue
        out.append((g.end - g.begin + 1) // 3 - 1)


def warn_if_eukaryotic_proteome(lengths) -> bool:
    """Log a WARNING when pyrodigal proteins look like a eukaryotic genome split at introns.

    pyrodigal cannot model introns: on a fungal genome it predicts many short single-exon ORFs
    (Xylona heveae: median 105 aa, 46% < 100 aa; E. coli: median 282 aa, 10% < 100 aa).
    Returns True when the warning was logged. The run is not stopped.
    """
    n = len(lengths)
    if n < EUK_WARN_MIN_GENES:
        return False
    ordered = sorted(lengths)
    median = ordered[n // 2] if n % 2 else (ordered[n // 2 - 1] + ordered[n // 2]) / 2
    short = sum(1 for x in ordered if x < EUK_WARN_SHORT_AA) / n
    if median < EUK_WARN_MEDIAN_AA and short > EUK_WARN_SHORT_FRACTION:
        logger.warning(
            f"pyrodigal proteins look like a eukaryotic genome split at introns: median {median:.0f} aa, "
            f"{100 * short:.0f}% shorter than {EUK_WARN_SHORT_AA} aa ({n} complete genes; bacteria are typically "
            f"~250-300 aa and ~10%). --mode prok/meta (pyrodigal) is for bacteria, archaea and metagenomes only. "
            "For a fungal or other eukaryotic genome, annotate it with funannotate/funannotate2 or BRAKER "
            "(or use the NCBI/JGI annotation) and run --mode protein on its proteins, with --input_gff and "
            "--gff_type funannotate/BRAKER/GFF3/GTF/NCBI_euk/JGI for CGCs."
        )
        return True
    return False


def clean_protein_sequence(seq: str):
    """Make a protein sequence acceptable to DIAMOND: strip trailing stop symbols ('*' or '.',
    as written by gffread -y) and replace internal ones with 'X'.

    Returns (sequence, trailing_stripped, n_internal_replaced).
    """
    stripped = seq.rstrip(PROTEIN_STOP_CHARS)
    trailing = len(stripped) != len(seq)
    internal = sum(stripped.count(c) for c in PROTEIN_STOP_CHARS)
    if internal:
        for c in PROTEIN_STOP_CHARS:
            stripped = stripped.replace(c, PROTEIN_UNKNOWN_AA)
    return stripped, trailing, internal


class MetaProcessor(ProkProcessor):
    """Processor for metagenomic data."""

    def process_input(self):
        return self.process_fna(True)


class ProteinProcessor(BaseProcessor):
    """Processor for protein sequence data."""

    @property
    def output_faa_path(self) -> Path:
        return self.output_dir_path / INPUT_PROTEIN_NAME

    def process_input(self):
        logger.info('Processing protein sequences')

        try:
            input_path = self.input_raw_data_path
            opener = gzip.open if input_path.suffix == '.gz' else open

            record_count = 0
            n_trailing = 0
            n_internal = 0
            logger.info(f"Reading protein sequences from {input_path}")

            with opener(input_path, "rt") as input_handle, self.output_faa_path.open("w") as output_handle:
                for record in SeqIO.parse(input_handle, "fasta"):
                    raw_id = record.id.split()[0]
                    if getattr(self.config, "simplify_ids", True):
                        clean_id = clean_sequence_id(raw_id)
                    else:
                        clean_id = raw_id
                    record.id = clean_id
                    record.description = ''
                    seq, trailing, internal = clean_protein_sequence(str(record.seq))
                    if trailing:
                        n_trailing += 1
                    if internal:
                        n_internal += 1
                    if trailing or internal:
                        record.seq = Seq(seq)
                    SeqIO.write(record, output_handle, "fasta")
                    record_count += 1

            if record_count == 0:
                logger.warning(f"No valid protein sequences found in {input_path}")
                return None

            if n_trailing or n_internal:
                logger.info(
                    f"Protein FASTA cleaning: removed a trailing stop ('*' or '.') from {n_trailing} sequences; "
                    f"replaced internal stops ('*' or '.') with 'X' in {n_internal} sequences "
                    "(DIAMOND rejects '.', e.g. from gffread -y without -S)."
                )
            logger.info(f"Processed {record_count} protein sequences to {self.output_faa_path}")

            if not self._verify_output(self.output_faa_path):
                return None

            # Ensure uniqueness (ProteinProcessor wrote directly; enforce uniqueness if duplicates created)
            try:
                _deduplicate_fasta_ids(self.output_faa_path)
            except Exception as ce:
                logger.warning(f"Post deduplicate skipped: {ce}")

            return str(self.output_faa_path)

        except Exception as e:
            logger.exception(f"Error processing protein sequences: {str(e)}")
            return None


def get_processor(config: GeneralConfig):
    """Factory function to get the appropriate processor based on mode."""
    mode = config.mode

    try:
        if mode == 'prok':
            return ProkProcessor(config)
        elif mode == 'meta':
            return MetaProcessor(config)
        elif mode == 'protein':
            return ProteinProcessor(config)
        else:
            raise ValueError(f"Unsupported mode: {mode}")
    except Exception as e:
        logger.error(f"Error creating processor for mode '{mode}': {str(e)}")
        raise

# ===================== ID Cleaning Utilities =====================

ID_MAX_LENGTH = 80
_prodigal_tail_pattern = re.compile(r'(?:_\d+){2,}$')  # kept (unused now, for possible future use)

def clean_sequence_id(raw_id: str, max_len: int = ID_MAX_LENGTH) -> str:
    """
    Clean a raw sequence ID.
    Rules:
      1. Keep only the first token (before first whitespace) - removes Prodigal comments.
      2. If contains 'JGI' (case-insensitive), replace every '|' with '-' to avoid delimiter conflicts.
      3. Truncate only if length exceeds max_len.
    Note: Preserves full Prodigal locus IDs like 'MGYG000290007_1_2' without modification.
    """
    if not raw_id:
        return raw_id
    
    # Split by whitespace and take first token only
    rid = raw_id.split()[0]
    
    # Handle JGI pipe conflicts
    if "JGI" in rid.upper():
        rid = rid.replace("|", "-")
    
    # Truncate if too long
    if len(rid) > max_len:
        rid = rid[:max_len]
    
    return rid

def _apply_id_cleaning(faa_path: Path, gff_path: Path | None = None):
    """
    Re-write FASTA (and optionally GFF) with cleaned IDs while preserving uniqueness.
    For GFF we replace protein_id=OLD with protein_id=NEW.
    Special handling for Pyrodigal output format.
    """
    faa_path = Path(faa_path)
    if not faa_path.exists():
        raise FileNotFoundError(f"FASTA not found for ID cleaning: {faa_path}")

    # First pass: collect mapping from full header to clean ID
    mapping = {}
    used = set()
    headers = []
    
    with faa_path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                # Get the full header (everything after >)
                full_hdr = line[1:].strip()
                headers.append(full_hdr)
    
    for full_hdr in headers:
        # Extract original ID - handle Pyrodigal format specially
        # Use the helper function to extract ID consistently
        original_id = _extract_id_from_header(full_hdr)
        
        # Clean the ID
        new_id = clean_sequence_id(original_id)
        
        # Handle collisions (rare, but ensure uniqueness)
        if new_id in used and new_id != original_id:
            suffix = 1
            base = new_id
            while f"{base}_{suffix}" in used:
                suffix += 1
            new_id = f"{base}_{suffix}"
        
        used.add(new_id)
        mapping[full_hdr] = new_id

    # Check if any changes are needed
    changes_needed = any(mapping[full_hdr] != _extract_id_from_header(full_hdr) for full_hdr in headers)
    if not changes_needed:
        logger.info("ID cleaning: no changes required (all IDs already clean)")
        return

    logger.info(f"ID cleaning: updating {sum(mapping[h] != _extract_id_from_header(h) for h in headers)} IDs (FASTA + GFF)")

    # Rewrite FASTA
    tmp_faa = faa_path.with_suffix(".tmp.clean")
    with faa_path.open() as fin, tmp_faa.open("w") as fout:
        for line in fin:
            if line.startswith(">"):
                full_hdr = line[1:].strip()
                new_id = mapping.get(full_hdr, _extract_id_from_header(full_hdr))
                fout.write(f">{new_id}\n")
            else:
                fout.write(line)
    shutil.move(tmp_faa, faa_path)

    # Rewrite GFF protein_id attributes if provided
    if gff_path:
        gff_path = Path(gff_path)
        if gff_path.exists():
            tmp_gff = gff_path.with_suffix(".tmp.clean")
            pid_pattern = re.compile(r'(protein_id=)([^;]+)')
            
            # Create reverse mapping from original_id to new_id
            id_mapping = {}
            for full_hdr, new_id in mapping.items():
                original_id = _extract_id_from_header(full_hdr)
                id_mapping[original_id] = new_id
            
            with gff_path.open() as fin, tmp_gff.open("w") as fout:
                for line in fin:
                    if line.startswith("#") or '\t' not in line:
                        fout.write(line)
                        continue
                    def repl(m):
                        old_id = m.group(2)
                        new_id = id_mapping.get(old_id, clean_sequence_id(old_id))
                        return f"{m.group(1)}{new_id}"
                    newline = pid_pattern.sub(repl, line)
                    fout.write(newline)
            shutil.move(tmp_gff, gff_path)

def _extract_id_from_header(full_header: str) -> str:
    """
    Extract protein ID from FASTA header, handling both Pyrodigal and standard formats.
    """
    if " # " in full_header:
        # Pyrodigal format: "MGYG000290007_1_2 # 1008 # 2369 # -1 # ID=18_2;..."
        return full_header.split(" # ")[0]
    else:
        # Standard format: take first token
        return full_header.split()[0]

def _deduplicate_fasta_ids(faa_path: Path):
    """
    Ensure FASTA IDs are unique (after direct write in ProteinProcessor).
    """
    faa_path = Path(faa_path)
    if not faa_path.exists():
        return
    tmp = faa_path.with_suffix(".tmp.dedup")
    seen = {}
    with faa_path.open() as fin, tmp.open("w") as fout:
        for line in fin:
            if line.startswith(">"):
                rid = line[1:].strip()
                base = rid
                if base in seen:
                    seen[base] += 1
                    rid = f"{base}_{seen[base]}"
                else:
                    seen[base] = 0
                fout.write(f">{rid}\n")
            else:
                fout.write(line)
    if any(v > 0 for v in seen.values()):
        shutil.move(tmp, faa_path)
    else:
        tmp.unlink(missing_ok=True)
