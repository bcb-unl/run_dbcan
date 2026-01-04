import logging
from pathlib import Path
from typing import Dict, Iterable, Set, List, Optional

import pandas as pd
from Bio import SeqIO

logger = logging.getLogger(__name__)
import dbcan.constants.process_utils_constants as PU
import dbcan.constants.base_constants as B


def _is_nonempty_file(path: Path) -> bool:
    return path.exists() and path.stat().st_size > 0

def _to_path(p) -> Path:
    return p if isinstance(p, Path) else Path(p)


def process_results(results: Optional[Iterable[List]], output_file, temp_hits_file: Optional[Path] = None) -> None:
    """Write HMMER hits to disk and filter highly-overlapping domains (per target).

    Accepts either an in-memory iterable of hit rows or a temporary hits TSV file.
    """
    out_path = _to_path(output_file)
    df = pd.DataFrame(columns=PU.HMMER_COLUMN_NAMES)

    try:
        if temp_hits_file and _is_nonempty_file(temp_hits_file):
            # Use chunked reading for large files (>100MB)
            file_size_mb = temp_hits_file.stat().st_size / (1024 * 1024)
            if file_size_mb > 100:
                logger.info(f"Large file detected ({file_size_mb:.1f}MB), using chunked reading")
                chunks = []
                for chunk in pd.read_csv(temp_hits_file, sep='\t', header=None, names=PU.HMMER_COLUMN_NAMES, chunksize=100000):
                    chunks.append(chunk)
                df = pd.concat(chunks, ignore_index=True)
            else:
                df = pd.read_csv(temp_hits_file, sep='\t', header=None, names=PU.HMMER_COLUMN_NAMES)
        elif results:
            df = pd.DataFrame(list(results), columns=PU.HMMER_COLUMN_NAMES)

        if df.empty:
            df.to_csv(out_path, index=False, sep='\t')
            return

        # enforce numeric types
        for col in [PU.TARGET_FROM_COLUMN, PU.TARGET_TO_COLUMN]:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors='coerce')
        if PU.IEVALUE_COLUMN in df.columns:
            df[PU.IEVALUE_COLUMN] = pd.to_numeric(df[PU.IEVALUE_COLUMN], errors='coerce')

        # basic cleaning
        df = df.dropna(subset=[PU.TARGET_NAME_COLUMN, PU.TARGET_FROM_COLUMN, PU.TARGET_TO_COLUMN])
        df.sort_values(by=[PU.TARGET_NAME_COLUMN, PU.TARGET_FROM_COLUMN, PU.TARGET_TO_COLUMN], inplace=True)
        df_filtered = filter_overlaps(df)
        df_filtered.to_csv(out_path, index=False, sep='\t')
    except Exception as e:
        logger.error(f"Error in process_results: {e}", exc_info=True)
        # write empty as fallback
        pd.DataFrame(columns=PU.HMMER_COLUMN_NAMES).to_csv(out_path, index=False, sep='\t')
    finally:
        if temp_hits_file:
            try:
                temp_hits_file.unlink(missing_ok=True)
            except Exception:
                logger.debug("Failed to remove temp hits file: %s", temp_hits_file)


def filter_overlaps(df: pd.DataFrame, overlap_ratio_threshold: float = PU.OVERLAP_RATIO_THRESHOLD) -> pd.DataFrame:
    """Within each Target Name, remove strongly-overlapping hits; keep better i-Evalue."""
    required = {PU.TARGET_NAME_COLUMN, PU.TARGET_FROM_COLUMN, PU.TARGET_TO_COLUMN}
    if not required.issubset(df.columns):
        missing = required - set(df.columns)
        logger.warning(f"filter_overlaps missing columns: {missing}, return original df")
        return df

    # types
    for col in [PU.TARGET_FROM_COLUMN, PU.TARGET_TO_COLUMN]:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    if PU.IEVALUE_COLUMN in df.columns:
        df[PU.IEVALUE_COLUMN] = pd.to_numeric(df[PU.IEVALUE_COLUMN], errors='coerce')

    df = df.dropna(subset=[PU.TARGET_NAME_COLUMN, PU.TARGET_FROM_COLUMN, PU.TARGET_TO_COLUMN])
    df = df.sort_values(by=[PU.TARGET_NAME_COLUMN, PU.TARGET_FROM_COLUMN, PU.TARGET_TO_COLUMN])

    kept_rows = []
    for _, group in df.groupby(PU.TARGET_NAME_COLUMN, sort=False):
        group = group.reset_index(drop=True)
        keep: List[pd.Series] = []
        for i in range(len(group)):
            if not keep:
                keep.append(group.iloc[i])
                continue

            last = keep[-1]
            cur = group.iloc[i]
            overlap = min(last[PU.TARGET_TO_COLUMN], cur[PU.TARGET_TO_COLUMN]) - max(last[PU.TARGET_FROM_COLUMN], cur[PU.TARGET_FROM_COLUMN])
            if overlap > 0:
                len_last = max(1, (last[PU.TARGET_TO_COLUMN] - last[PU.TARGET_FROM_COLUMN]))
                len_cur = max(1, (cur[PU.TARGET_TO_COLUMN] - cur[PU.TARGET_FROM_COLUMN]))
                r_last = overlap / len_last
                r_cur = overlap / len_cur
                if r_last > overlap_ratio_threshold or r_cur > overlap_ratio_threshold:
                    last_eval = last[PU.IEVALUE_COLUMN] if PU.IEVALUE_COLUMN in last else float('inf')
                    cur_eval = cur[PU.IEVALUE_COLUMN] if PU.IEVALUE_COLUMN in cur else float('inf')
                    try:
                        if float(last_eval) > float(cur_eval):
                            keep[-1] = cur
                    except Exception:
                        # keep earlier if not comparable
                        pass
                else:
                    keep.append(cur)
            else:
                keep.append(cur)
        kept_rows.extend(keep)

    return pd.DataFrame(kept_rows, columns=df.columns)


def process_cgc_sig_results(tc_config, tfdiamond_config, tf_config, stp_config, sulfatase_config, peptidase_config) -> None:
    """Combine TCDB, TF and STP results into one file and filter overlaps."""
    try:
        columns = PU.CGC_SIG_RESULT_COLUMN
        # resolve output_dir from any given config
        output_dir = getattr(tc_config, 'output_dir', None)
        out_dir_path = Path(output_dir)
        out_dir_path.mkdir(parents=True, exist_ok=True)

        # build output file paths
        output_files = {
            'TC': out_dir_path / B.TCDB_DIAMOND_OUTPUT,
            'STP': out_dir_path / B.STP_HMM_RESULT_FILE,
            'Sulfatase': out_dir_path / B.SULFATLAS_DIAMOND_OUTPUT,
            'Peptidase': out_dir_path / B.PEPTIDASE_DIAMOND_OUTPUT,
        }
        # add TF variants based on flags
        if tfdiamond_config is not None and getattr(tfdiamond_config, 'prokaryotic', True):
            output_files['TF_prok'] = out_dir_path / B.TF_DIAMOND_OUTPUT
        if tf_config is not None and getattr(tf_config, 'fungi', False):
            output_files['TF_fungi'] = out_dir_path / B.TF_HMM_RESULT_FILE

        frames = []
        for name, fpath in output_files.items():
            if not _is_nonempty_file(fpath):
                logger.warning(f"{name} output file not found or empty: {fpath}")
                continue
            try:
                # Use chunked reading for large files (>100MB)
                file_size_mb = fpath.stat().st_size / (1024 * 1024)
                if file_size_mb > 100:
                    logger.info(f"Large file detected ({file_size_mb:.1f}MB), using chunked reading for {name}")
                    chunks = []
                    for chunk in pd.read_csv(fpath, sep='\t', header=0, chunksize=100000):
                        chunks.append(chunk)
                    df = pd.concat(chunks, ignore_index=True)
                else:
                    df = pd.read_csv(fpath, sep='\t', header=0)
                
                # align columns
                missing = [c for c in columns if c not in df.columns]
                if missing:
                    if len(df.columns) >= len(columns):
                        df = df.iloc[:, :len(columns)]
                        df.columns = columns
                    else:
                        logger.warning(f"{name} columns missing {missing}, skip file: {fpath}")
                        continue
                else:
                    df = df[columns]

                for col in [PU.TARGET_FROM_COLUMN, PU.TARGET_TO_COLUMN, PU.ANNOTATE_FROM_COLUMN, PU.ANNOTATE_TO_COLUMN, PU.COVERAGE_COLUMN]:
                    if col in df.columns:
                        df[col] = pd.to_numeric(df[col], errors='coerce')

                if PU.IEVALUE_COLUMN in df.columns:
                    df[PU.IEVALUE_COLUMN] = pd.to_numeric(df[PU.IEVALUE_COLUMN], errors='coerce')
                df = df.dropna(subset=[PU.TARGET_NAME_COLUMN, PU.TARGET_FROM_COLUMN, PU.TARGET_TO_COLUMN])
                df.sort_values(by=[PU.TARGET_NAME_COLUMN, PU.TARGET_FROM_COLUMN, PU.TARGET_TO_COLUMN], inplace=True)

                df['Type'] = name
                frames.append(df)
                logger.info(f"Loaded {len(df)} {name} annotations from {fpath}")
            except Exception as e:
                logger.error(f"Error reading {name} output file: {e}", exc_info=True)

        out_file = out_dir_path / B.GFF_CGC_SIG_FILE
        if not frames:
            logger.warning("No valid CGC annotation data found")
            pd.DataFrame(columns=columns + ['Type']).to_csv(out_file, index=False, sep='\t')
            return

        # Optimize DataFrame concatenation for large datasets
        # Concatenate in chunks if total size is large
        total_rows = sum(len(f) for f in frames)
        if total_rows > 1000000:  # 1M rows threshold
            logger.info(f"Large dataset detected ({total_rows} rows), concatenating in chunks")
            # Concatenate in smaller groups to reduce memory peak
            chunk_size = 3  # Concatenate 3 dataframes at a time
            result_frames = []
            for i in range(0, len(frames), chunk_size):
                chunk = frames[i:i+chunk_size]
                result_frames.append(pd.concat(chunk, ignore_index=True))
            total_df = pd.concat(result_frames, ignore_index=True)
        else:
            total_df = pd.concat(frames, ignore_index=True)
        
        filtered_df = filter_overlaps(total_df)
        filtered_df.to_csv(out_file, index=False, sep='\t')
        logger.info(f"Saved {len(filtered_df)} CGC annotations to {out_file}")
    except Exception as e:
        logger.error(f"Error processing CGC signature results: {e}", exc_info=True)


def process_cgc_null_pfam_annotation(Pfam_config) -> None:
    """Post-process PFAM HMM null annotations for CGC pipeline outputs."""
    output_dir = Path(getattr(Pfam_config, 'output_dir', '.'))
    pfam_hmm_output = output_dir / B.PFAM_HMM_RESULT_FILE

    if not _is_nonempty_file(pfam_hmm_output):
        logger.warning(f"PFAM HMM output file not found or empty: {pfam_hmm_output}")
        pd.DataFrame(columns=PU.HMMER_COLUMN_NAMES).to_csv(pfam_hmm_output, index=False, sep='\t')
        return

    try:
        # Use chunked reading for large files (>100MB)
        file_size_mb = pfam_hmm_output.stat().st_size / (1024 * 1024)
        if file_size_mb > 100:
            logger.info(f"Large file detected ({file_size_mb:.1f}MB), using chunked reading")
            chunks = []
            for chunk in pd.read_csv(pfam_hmm_output, sep='\t', header=0, chunksize=100000):
                chunks.append(chunk)
            df = pd.concat(chunks, ignore_index=True)
        else:
            df = pd.read_csv(pfam_hmm_output, sep='\t', header=0)
        
        # align to HMMER_COLUMN_NAMES
        missing = [c for c in PU.HMMER_COLUMN_NAMES if c not in df.columns]
        if missing:
            if len(df.columns) >= len(PU.HMMER_COLUMN_NAMES):
                df = df.iloc[:, :len(PU.HMMER_COLUMN_NAMES)]
                df.columns = PU.HMMER_COLUMN_NAMES
            else:
                logger.warning(f"PFAM HMM output columns missing {missing}, skip filtering")
        for col in [PU.TARGET_FROM_COLUMN, PU.TARGET_TO_COLUMN]:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors='coerce')
        if PU.IEVALUE_COLUMN in df.columns:
            df[PU.IEVALUE_COLUMN] = pd.to_numeric(df[PU.IEVALUE_COLUMN], errors='coerce')
        df = df.dropna(subset=[PU.TARGET_NAME_COLUMN, PU.TARGET_FROM_COLUMN, PU.TARGET_TO_COLUMN])
        df.sort_values(by=[PU.TARGET_NAME_COLUMN, PU.TARGET_FROM_COLUMN, PU.TARGET_TO_COLUMN], inplace=True)

        df_filtered = filter_overlaps(df)
        df_filtered.to_csv(pfam_hmm_output, index=False, sep='\t')
        logger.info(f"Processed PFAM HMM results saved to: {pfam_hmm_output}")
    except Exception as e:
        logger.error(f"Error processing CGC null PFAM annotation: {e}", exc_info=True)


def extract_null_protein_ids_from_cgc(cgc_standard_out_file) -> Set[str]:
    """Extract protein IDs with 'null' type in cgc_standard_out.tsv."""
    path = _to_path(cgc_standard_out_file)
    if not path.exists():
        logger.warning(f"CGC standard out not found: {path}")
        return set()

    null_ids: Set[str] = set()
    with path.open() as f:
        header_seen = False
        for line in f:
            if not header_seen:
                header_seen = True
                continue
            if line.startswith("#") or not line.strip():
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 4 and fields[1].lower() == "null":
                null_ids.add(fields[3])
    return null_ids


def extract_fasta_by_protein_ids(input_faa, output_faa, protein_ids: Iterable[str]) -> int:
    """Extract sequences by IDs from input_faa to output_faa."""
    in_path = _to_path(input_faa)
    out_path = _to_path(output_faa)

    if not in_path.exists():
        logger.warning(f"Protein fasta not found: {in_path}")
        return 0

    count = 0
    with out_path.open("w") as out_handle:
        for record in SeqIO.parse(str(in_path), "fasta"):
            if record.id in protein_ids:
                SeqIO.write(record, out_handle, "fasta")
                count += 1
    return count


def extract_null_fasta_from_cgc(cgc_standard_out_file, input_faa, output_faa) -> None:
    """Extract FASTA for null protein IDs defined in cgc_standard_out.tsv."""
    null_protein_ids = extract_null_protein_ids_from_cgc(cgc_standard_out_file)
    count = extract_fasta_by_protein_ids(input_faa, output_faa, null_protein_ids)
    logger.info(f"Extracted {count} null protein sequences to {output_faa}")


def annotate_cgc_null_with_pfam_and_gff(cgc_standard_out_file, pfam_hmm_result_file, gff_file, output_cgc_file, output_gff_file) -> None:
    """Annotate CGC null entries with Pfam annotations in both cgc_standard_out.tsv and cgc.gff files."""
    pfam_path = _to_path(pfam_hmm_result_file)
    cgc_path = _to_path(cgc_standard_out_file)
    gff_path = _to_path(gff_file)
    out_cgc_path = _to_path(output_cgc_file)
    out_gff_path = _to_path(output_gff_file)

    if not pfam_path.exists():
        logger.warning(f"PFAM HMM result file not found: {pfam_path}")
        return
    if not cgc_path.exists():
        logger.warning(f"CGC standard out not found: {cgc_path}")
        return
    if not gff_path.exists():
        logger.warning(f"GFF file not found: {gff_path}")
        return

    # build mapping: protein_id -> pfam annotation
    pfam_map: Dict[str, str] = {}
    # Use chunked reading for large files (>100MB)
    file_size_mb = pfam_path.stat().st_size / (1024 * 1024)
    if file_size_mb > 100:
        logger.info(f"Large file detected ({file_size_mb:.1f}MB), using chunked reading for mapping")
        for chunk in pd.read_csv(pfam_path, sep='\t', header=0, chunksize=100000):
            if not {PU.TARGET_NAME_COLUMN, PU.HMM_NAME_COLUMN}.issubset(chunk.columns):
                logger.warning(f"PFAM result missing required columns in {pfam_path}")
                return
            for _, row in chunk.iterrows():
                protein_id = str(row[PU.TARGET_NAME_COLUMN])
                annot = str(row[PU.HMM_NAME_COLUMN])
                pfam_map[protein_id] = annot
    else:
        df = pd.read_csv(pfam_path, sep='\t', header=0)
        if not {PU.TARGET_NAME_COLUMN, PU.HMM_NAME_COLUMN}.issubset(df.columns):
            logger.warning(f"PFAM result missing required columns in {pfam_path}")
            return
        for _, row in df.iterrows():
            protein_id = str(row[PU.TARGET_NAME_COLUMN])
            annot = str(row[PU.HMM_NAME_COLUMN])
            pfam_map[protein_id] = annot

    # update cgc_standard_out.tsv
    with cgc_path.open() as fin, out_cgc_path.open('w') as fout:
        header = fin.readline()
        fout.write(header)
        for line in fin:
            if line.startswith("#") or not line.strip():
                fout.write(line)
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) >= 8 and fields[1].lower() == "null":
                protein_id = fields[3]
                if protein_id in pfam_map:
                    fields[1] = "Pfam"
                    fields[7] = pfam_map[protein_id]  # Gene Annotation column
            fout.write('\t'.join(fields) + '\n')

    # update cgc.gff
    with gff_path.open() as fin, out_gff_path.open('w') as fout:
        for line in fin:
            if line.startswith("#") or not line.strip():
                fout.write(line)
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9:
                fout.write(line)
                continue
            attr = fields[8]
            attr_dict: Dict[str, str] = {}
            for item in attr.split(';'):
                if '=' in item:
                    k, v = item.split('=', 1)
                    attr_dict[k] = v
            protein_id = attr_dict.get(B.GFF_PROTEIN_ID_COL, None)
            if protein_id and attr_dict.get(B.CGC_ANNOTATION_COLUMN, '').lower() == 'null' and protein_id in pfam_map:
                attr_dict[B.CGC_ANNOTATION_COLUMN] = f"Pfam|{pfam_map[protein_id]}"
            new_attr = ';'.join([f"{k}={v}" for k, v in attr_dict.items()])
            fields[8] = new_attr
            fout.write('\t'.join(fields) + '\n')


def extract_null_protein_ids_from_cgc_gff(cgc_gff_file) -> Set[str]:
    """Extract protein_id values from cgc.gff where CGC_annotation=null."""
    path = _to_path(cgc_gff_file)
    null_ids: Set[str] = set()
    if not path.exists():
        logger.warning(f"GFF file not found: {path}")
        return null_ids

    with path.open() as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            attrs: Dict[str, str] = {}
            for kv in parts[8].split(';'):
                if '=' in kv:
                    k, v = kv.split('=', 1)
                    attrs[k] = v
            if attrs.get(B.CGC_ANNOTATION_COLUMN, '').lower() == 'null':
                pid = attrs.get(B.GFF_PROTEIN_ID_COL)
                if pid:
                    null_ids.add(pid)
    return null_ids


def extract_null_fasta_from_gff(cgc_gff_file, input_faa, output_faa) -> None:
    """Extract FASTA for null protein IDs from cgc.gff."""
    null_protein_ids = extract_null_protein_ids_from_cgc_gff(cgc_gff_file)
    count = extract_fasta_by_protein_ids(input_faa, output_faa, null_protein_ids)
    logger.info(f"Extracted {count} null protein sequences (from GFF) to {output_faa}")
