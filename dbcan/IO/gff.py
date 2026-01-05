import logging
import sys
import multiprocessing
from pathlib import Path
from typing import Dict, Any

import pandas as pd
from Bio import SeqIO
from BCBio import GFF
from tqdm import tqdm

from dbcan.configs.base_config import GFFConfig  # keep as-is unless you migrate configs
import dbcan.constants.gff_constants as G  # for __all__

logger = logging.getLogger(__name__)

# Utility: read GFF to DataFrame
def read_gff_to_df(path: str, columns=None, drop_comments: bool = True) -> pd.DataFrame:
    """
    Read a GFF file into a DataFrame with the standard 9 columns.
    If columns is None, uses GFF_COLUMNS from constants.

    Important:
    - Do NOT treat '#' as inline comment, because attributes/protein_id may contain '#'.
    - We only drop lines that START WITH '#', which are true comment lines in GFF.
    """
    if columns is None:
        columns = G.GFF_COLUMNS
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"GFF file not found: {p}")

    # Never treat '#' as inline comment; some protein_id contain '#'
    # Use the python engine for robust parsing of mixed content.
    df = pd.read_csv(
        p,
        sep='\t',
        header=None,
        names=columns,
        comment=None,      # critical: keep inline '#'
        engine='python'
    )

    if drop_comments:
        # Only drop lines that start with '#' (true comment lines)
        first_col = columns[0]
        df = df[~df[first_col].astype(str).str.startswith('#')]

    # Sanity check: if many rows have 'protein_id=.*#' but no 'CGC_annotation=', warn the caller.
    attr_col = columns[-1]
    attrs_str = df[attr_col].astype(str)
    has_hash_pid = attrs_str.str.contains(r'protein_id=[^;]*#').any()
    has_cgc_ann = attrs_str.str.contains(r'\bCGC_annotation=').any()
    if has_hash_pid and not has_cgc_ann:
        logger.warning(
            "GFF appears to contain inline '#' in protein_id but no 'CGC_annotation=' was detected. "
            "Make sure no upstream reader uses comment='#'."
        )

    return df

def read_gff(path: str) -> pd.DataFrame:
    """Back-compat API: return DataFrame."""
    return read_gff_to_df(path, columns=G.GFF_COLUMNS, drop_comments=True)

def _sort_record_features(record):
    """Worker function to sort features within a single GFF record."""
    sorted_features = sorted(record.features, key=lambda f: (f.location.start, f.location.end))
    record.features = sorted_features
    return record


def _sort_features_in_record(record):
    """Return record with its features sorted by genomic start (then end)."""
    record.features = sorted(
        record.features,
        key=lambda f: (
            int(getattr(f.location, "start", 0)),
            int(getattr(f.location, "end", 0))
        )
    )
    return record


def _natural_key(s: str):
    """Natural sort helper: splits digits and text."""
    import re
    return [int(t) if t.isdigit() else t for t in re.findall(r'\d+|[^\d]+', s)]


class GFFProcessor:
    """Base GFF processor. Holds only config and derives paths via properties."""

    def __init__(self, config: GFFConfig):
        self.config = config
        self._validate_basic()

    # Properties (single source of truth: config)
    @property
    def output_dir_path(self) -> Path:
        return Path(self.config.output_dir)

    @property
    def input_gff_path(self) -> Path:
        return Path(self.config.input_gff)

    @property
    def input_total_faa_path(self) -> Path:
        return self.output_dir_path / G.GFF_INPUT_PROTEIN_FILE

    @property
    def cazyme_overview_path(self) -> Path:
        return self.output_dir_path / G.GFF_CAZYME_OVERVIEW_FILE

    @property
    def cgc_sig_path(self) -> Path:
        return self.output_dir_path / G.GFF_CGC_SIG_FILE

    @property
    def output_gff_path(self) -> Path:
        return self.output_dir_path / G.GFF_OUTPUT_FILE

    def _validate_basic(self):
        # Only minimal upfront checks
        self.output_dir_path.mkdir(parents=True, exist_ok=True)

    def load_cgc_type(self) -> Dict[str, Dict[str, Any]]:
        """
        Load CAZyme overview and CGC signature files and build a mapping:
        protein_id -> { GFF_CGC_ANNOTATION_COL: annotation }
        """
        try:
            if not self.cazyme_overview_path.exists():
                logger.error(f"CAZyme overview file not found: {self.cazyme_overview_path}")
                return {}

            if not self.cgc_sig_path.exists():
                logger.error(f"CGC signature file not found: {self.cgc_sig_path}")
                return {}

            logger.info(f"Loading CAZyme data from {self.cazyme_overview_path}")
            df = pd.read_csv(self.cazyme_overview_path, sep='\t')

            required_cols = [G.GFF_GENE_ID_COL, G.GFF_TOOLS_COUNT_COL, G.GFF_RECOMMEND_RESULTS_COL]
            if not all(col in df.columns for col in required_cols):
                logger.error(f"Missing required columns in CAZyme overview. Found: {df.columns.tolist()}")
                return {}

            df[G.GFF_GENE_ID_COL] = df[G.GFF_GENE_ID_COL].astype(str)
            df[G.GFF_RECOMMEND_RESULTS_COL] = df[G.GFF_RECOMMEND_RESULTS_COL].astype(str)
            df = df.rename(columns={G.GFF_GENE_ID_COL: G.GFF_PROTEIN_ID_COL, G.GFF_RECOMMEND_RESULTS_COL: G.GFF_CAZYME_COL})

            df[G.GFF_TOOLS_COUNT_COL] = pd.to_numeric(df[G.GFF_TOOLS_COUNT_COL], errors='coerce')
            overview_df = (
                df[df[G.GFF_TOOLS_COUNT_COL] >= G.GFF_MIN_TOOL_COUNT][[G.GFF_PROTEIN_ID_COL, G.GFF_CAZYME_COL]].drop_duplicates()
            )
            overview_df[G.GFF_CGC_ANNOTATION_COL] = G.GFF_CAZYME_PREFIX + overview_df[G.GFF_CAZYME_COL].astype(str)

            logger.info(f"Loading CGC signature data from {self.cgc_sig_path}")
            try:
                cgc_sig_chunks = pd.read_csv(
                    self.cgc_sig_path,
                    sep='\t',
                    usecols=G.GFF_CGC_SIG_COLUMNS,
                    header=None,
                    names=[G.GFF_FUNCTION_ANNOTATION_COL, G.GFF_PROTEIN_ID_COL, G.GFF_TYPE_COL],
                    chunksize=100_000
                )
                cgc_sig_df_list = []
                chunk_count = 0
                for chunk in cgc_sig_chunks:
                    chunk_count += 1
                    for col in (G.GFF_FUNCTION_ANNOTATION_COL, G.GFF_PROTEIN_ID_COL, G.GFF_TYPE_COL):
                        chunk[col] = chunk[col].astype(str)
                    before_rows = len(chunk)
                    chunk = chunk.dropna(subset=[G.GFF_PROTEIN_ID_COL])
                    chunk[G.GFF_CGC_ANNOTATION_COL] = chunk[G.GFF_TYPE_COL] + '|' + chunk[G.GFF_FUNCTION_ANNOTATION_COL]
                    cgc_sig_df_list.append(chunk[[G.GFF_PROTEIN_ID_COL, G.GFF_CGC_ANNOTATION_COL]])
                cgc_sig_df = pd.concat(cgc_sig_df_list, ignore_index=True) if cgc_sig_df_list else pd.DataFrame(
                    columns=[G.GFF_PROTEIN_ID_COL, G.GFF_CGC_ANNOTATION_COL]
                )
                logger.info(f"CGC signature chunks read: {chunk_count}, total rows: {len(cgc_sig_df)}; "
                            f"TC rows: {(cgc_sig_df[G.GFF_CGC_ANNOTATION_COL].str.startswith('TC|').sum() if not cgc_sig_df.empty else 0)}")
            except Exception as e:
                logger.error(f"Error reading CGC signature file: {e}")
                return overview_df.set_index(G.GFF_PROTEIN_ID_COL).to_dict('index')

            logger.info("Combining CAZyme and CGC signature data")
            combined_df = pd.concat(
                [
                    overview_df[[G.GFF_PROTEIN_ID_COL, G.GFF_CGC_ANNOTATION_COL]],
                    cgc_sig_df[[G.GFF_PROTEIN_ID_COL, G.GFF_CGC_ANNOTATION_COL]],
                ],
                ignore_index=True,
            ).drop_duplicates()

            combined_df = combined_df.groupby(G.GFF_PROTEIN_ID_COL, sort=False)[G.GFF_CGC_ANNOTATION_COL] \
                                     .apply(lambda x: '+'.join(list(x))) \
                                     .reset_index()
            lost_tc = (combined_df[G.GFF_CGC_ANNOTATION_COL].str.contains(r'TC\|', regex=True).sum())
            logger.info(f"Combined annotations: {len(combined_df)} proteins; with TC: {lost_tc}")
            return combined_df.set_index(G.GFF_PROTEIN_ID_COL).to_dict('index')

        except Exception as e:
            logger.error(f"Error loading CGC data: {e}")
            import traceback
            traceback.print_exc()
            return {}

    def _is_prodigal_mode(self) -> bool:
        """
        Return True only when gff_type equals the prodigal format constant.
        (User clarified: decision solely by gff_type, not mode.)
        """
        gff_type = str(getattr(self.config, "gff_type", "")).lower()
        prodigal_const = str(getattr(G, "GFF_FORMAT_PRODIGAL", "prodigal")).lower()
        return gff_type == prodigal_const

    def process_gff(self) -> bool:
        """
        Pipeline:
          - If gff_type == prodigal: skip sorting (preserve original contig & feature order).
          - Else: sort features within each contig, preserve contig order.
          - Load CAZyme / CGC annotations.
        """
        temp_output_file = str(self.output_gff_path) + G.GFF_TEMP_SUFFIX
        try:
            if not self.input_gff_path.exists():
                logger.error(f"Input GFF file not found: {self.input_gff_path}")
                return False

            do_sort = not self._is_prodigal_mode()
            if do_sort:
                logger.info("[GFF] Sorting enabled (non-prodigal). Feature-level sort only per contig.")
                self.sort_gff(self.input_gff_path, temp_output_file, natural_sort=False)
                working_input = temp_output_file
            else:
                logger.info("[GFF] Prodigal type detected: skip sorting, use raw order.")
                working_input = str(self.input_gff_path)

            logger.info("Loading annotation tables (CAZyme / CGC signatures)")
            cgc_data = self.load_cgc_type()
            processed_cgc_data = self._preprocess_cgc_data(cgc_data)

            logger.info(f"Annotating GFF from: {working_input}")
            self._process_gff_format(working_input, str(self.output_gff_path), processed_cgc_data)

            if do_sort:
                tp = Path(temp_output_file)
                if tp.exists():
                    tp.unlink()

            if not self.output_gff_path.exists():
                logger.error("Expected output GFF not created.")
                return False

            logger.info(f"GFF processing complete: {self.output_gff_path}")
            return True
        except Exception as e:
            logger.error(f"Error processing GFF: {e}")
            import traceback; traceback.print_exc()
            if Path(temp_output_file).exists():
                Path(temp_output_file).unlink()
            return False

    # Hooks for subclasses
    def _preprocess_cgc_data(self, cgc_data: Dict[str, Dict[str, Any]]) -> Dict[str, Dict[str, Any]]:
        """Subclasses may override this to normalize protein IDs etc."""
        return cgc_data

    def _process_gff_format(self, input_file: str, output_file: str, cgc_data: Dict[str, Dict[str, Any]]):
        """Must be implemented by subclasses for specific GFF flavors."""
        raise NotImplementedError

    # (sort_gff unchanged except you can update docstring to emphasize behavior)
    def sort_gff(
        self,
        input_gff: str,
        output_gff: str,
        natural_sort: bool = False
    ):
        """
        Sort features inside each record by (start,end).
        Preserve record (contig) order exactly as parsed unless natural_sort=True.
        No cross-record reordering by default.
        """
        in_path = Path(input_gff)
        out_path = Path(output_gff)
        if not in_path.exists():
            raise FileNotFoundError(f"Input GFF file not found: {in_path}")

        logger.info(
            f"Sorting GFF (preserve contig order={not natural_sort}; feature-level sort only). "
            f"Input={in_path}"
        )

        records = []
        try:
            with in_path.open() as fh:
                for idx, rec in enumerate(GFF.parse(fh)):
                    _sort_features_in_record(rec)
                    records.append((idx, rec))
        except Exception as e:
            logger.error(f"Error parsing GFF for sorting: {e}")
            raise

        if natural_sort:
            logger.info("Applying natural sort to record IDs (user-enabled).")
            records.sort(key=lambda x: _natural_key(x[1].id))
        else:
            # Keep the original parse order (idx ascending)
            records.sort(key=lambda x: x[0])

        with out_path.open('w') as out_fh:
            for _, rec in records:
                GFF.write([rec], out_fh)

        logger.info(
            f"Sorted GFF written: {out_path} (records={len(records)})"
        )

    def write_gff(self,
                  record,
                  feature,
                  protein_id: str,
                  cgc_annotation: str,
                  out_file):
        """
        Write a single annotated GFF line.
        - Skip 'remark' (any case)
        - Attributes: protein_id ; CGC_annotation
        """
        try:
            if not hasattr(feature, "location") or feature.location is None:
                logger.debug("Skip feature without location")
                return
            ftype_raw = getattr(feature, "type", "")
            if ftype_raw and ftype_raw.lower() == "remark":
                return

            start = int(feature.location.start) + 1  # GFF is 1-based inclusive
            end = int(feature.location.end)
            strand_val = getattr(feature.location, "strand", None)
            strand = '+' if strand_val == 1 else '-' if strand_val == -1 else '.'

            protein_id = protein_id if protein_id else G.GFF_UNKNOWN_ANNOTATION
            cgc_annotation = cgc_annotation if cgc_annotation else G.GFF_NULL_ANNOTATION

            line = (
                f"{record.id}\t.\t{feature.type}\t{start}\t{end}\t.\t{strand}\t.\t"
                f"{G.GFF_PROTEIN_ID_ATTR_NCBI}={protein_id};"
                f"{G.GFF_CGC_ANNOTATION_COL}={cgc_annotation}\n"
            )
            out_file.write(line)
        except Exception as e:
            logger.error(f"Error writing GFF line: {e}")


class NCBIEukProcessor(GFFProcessor):
    """Processor for NCBI Eukaryotic GFF format"""

    def _process_gff_format(self, input_file, output_file, cgc_data):
        try:
            with Path(input_file).open() as in_file, Path(output_file).open('w') as out_file:
                for record in tqdm(GFF.parse(in_file), desc="Generating NCBI Euk GFF"):
                    for feature in record.features:
                        if feature.type == G.GFF_GENE_FEATURE:
                            protein_id = G.GFF_UNKNOWN_ANNOTATION
                            cgc_annotation = G.GFF_UNKNOWN_ANNOTATION
                            non_mrna_found = False

                            for sub_feature in getattr(feature, "sub_features", []):
                                if sub_feature.type != G.GFF_MRNA_FEATURE:
                                    protein_id = G.GFF_NA_PROTEIN_ID
                                    cgc_annotation = G.GFF_OTHER_PREFIX + sub_feature.type
                                    non_mrna_found = True
                                    break

                            if non_mrna_found:
                                self.write_gff(record, feature, protein_id, cgc_annotation, out_file)
                                continue

                            for sub_feature in getattr(feature, "sub_features", []):
                                if sub_feature.type == G.GFF_MRNA_FEATURE:
                                    for sub_sub_feature in getattr(sub_feature, "sub_features", []):
                                        if sub_sub_feature.type == G.GFF_CDS_FEATURE:
                                            protein_id = sub_sub_feature.qualifiers.get(G.GFF_PROTEIN_ID_ATTR_NCBI, [G.GFF_UNKNOWN_ANNOTATION])[0]
                                            break
                                    if protein_id != G.GFF_UNKNOWN_ANNOTATION:
                                        break

                            cgc_annotation = cgc_data.get(protein_id, {}).get(G.GFF_CGC_ANNOTATION_COL, G.GFF_NULL_ANNOTATION)
                            self.write_gff(record, feature, protein_id, cgc_annotation, out_file)
        except Exception as e:
            logger.error(f"Error processing NCBI Eukaryotic GFF: {e}")
            raise


class NCBIProkProcessor(GFFProcessor):
    """Processor for NCBI Prokaryotic GFF format"""

    def _process_gff_format(self, input_file, output_file, cgc_data):
        try:
            with Path(input_file).open() as in_file, Path(output_file).open('w') as out_file:
                for record in tqdm(GFF.parse(in_file), desc="Generating NCBI Prok GFF"):
                    for feature in record.features:
                        if feature.type == G.GFF_GENE_FEATURE:
                            protein_id = G.GFF_UNKNOWN_ANNOTATION
                            cgc_annotation = G.GFF_UNKNOWN_ANNOTATION
                            non_cds_found = False

                            for sub_feature in getattr(feature, "sub_features", []):
                                if G.GFF_CDS_FEATURE not in sub_feature.type:
                                    protein_id = G.GFF_NA_PROTEIN_ID
                                    cgc_annotation = G.GFF_OTHER_PREFIX + sub_feature.type
                                    non_cds_found = True
                                    break

                            if non_cds_found:
                                self.write_gff(record, feature, protein_id, cgc_annotation, out_file)
                                continue

                            for sub_feature in getattr(feature, "sub_features", []):
                                if sub_feature.type == G.GFF_CDS_FEATURE:
                                    protein_id = sub_feature.qualifiers.get(G.GFF_PROTEIN_ID_ATTR_NCBI, [G.GFF_UNKNOWN_ANNOTATION])[0]
                                if protein_id != G.GFF_UNKNOWN_ANNOTATION:
                                    break

                            cgc_annotation = cgc_data.get(protein_id, {}).get(G.GFF_CGC_ANNOTATION_COL, G.GFF_NULL_ANNOTATION)
                            self.write_gff(record, feature, protein_id, cgc_annotation, out_file)
        except Exception as e:
            logger.error(f"Error processing NCBI Prokaryotic GFF: {e}")
            raise


class JGIProcessor(GFFProcessor):
    """Processor for JGI GFF format"""

    def _preprocess_cgc_data(self, cgc_data):
        # JGI uses a different protein ID format - extract the actual ID (3rd token by '-', generally it is split by '|' but we modify it to '-')
        return {k.split('-')[2] if '-' in k else k: v for k, v in cgc_data.items()}

    def _process_gff_format(self, input_file, output_file, cgc_data):
        try:
            # Build original ID mapping from the FASTA file (if present)
            original_id_mapping = {}
            if self.input_total_faa_path.exists():
                for record in SeqIO.parse(self.input_total_faa_path, 'fasta'):
                    original_id = record.id
                    simplified_id = original_id.split('-')[2] if '-' in original_id else original_id
                    original_id_mapping[simplified_id] = original_id
            else:
                logger.error(f"Input protein sequences file not found: {self.input_total_faa_path}")

            with Path(input_file).open() as in_file, Path(output_file).open('w') as out_file:
                for record in tqdm(GFF.parse(in_file), desc="Generating JGI GFF"):
                    for feature in record.features:
                        if feature.type == G.GFF_GENE_FEATURE:
                            protein_id = feature.qualifiers.get(G.GFF_JGI_PROTEIN_ID_ATTR, [G.GFF_UNKNOWN_ANNOTATION])[0]
                            simplified_id = protein_id.split('-')[2] if '-' in protein_id else protein_id
                            cgc_annotation = cgc_data.get(simplified_id, {}).get(G.GFF_CGC_ANNOTATION_COL, G.GFF_NULL_ANNOTATION)
                            original_protein_id = original_id_mapping.get(simplified_id, protein_id)
                            self.write_gff(record, feature, original_protein_id, cgc_annotation, out_file)
        except Exception as e:
            logger.error(f"Error processing JGI GFF: {e}")
            raise


class ProdigalProcessor(GFFProcessor):
    """Processor for Prodigal GFF format (features often CDS)."""

    def _process_gff_format(self, input_file, output_file, cgc_data):
        try:
            wrote = 0
            with Path(input_file).open() as in_file, Path(output_file).open('w') as out_file:
                for record in tqdm(GFF.parse(in_file), desc="Generating Prodigal GFF"):
                    for feature in record.features:
                        # Prodigal usually puts IDs in 'ID' qualifier
                        protein_id = feature.qualifiers.get(G.GFF_ID_ATTR_PRODIGAL, [G.GFF_UNKNOWN_ANNOTATION])[0]
                        cgc_annotation = cgc_data.get(protein_id, {}).get(
                            G.GFF_CGC_ANNOTATION_COL,
                            G.GFF_NULL_ANNOTATION
                        )
                        self.write_gff(record, feature, protein_id, cgc_annotation, out_file)
                        wrote += 1
            logger.info(f"[Prodigal] Written feature lines: {wrote}")
        except Exception as e:
            logger.error(f"Error processing Prodigal GFF: {e}")
            raise


def get_gff_processor(config: GFFConfig) -> GFFProcessor:
    """Factory: build a processor by gff_type."""
    gff_type = config.gff_type
    if gff_type == G.GFF_FORMAT_NCBI_EUK:
        return NCBIEukProcessor(config)
    elif gff_type == G.GFF_FORMAT_NCBI_PROK:
        return NCBIProkProcessor(config)
    elif gff_type == G.GFF_FORMAT_JGI:
        return JGIProcessor(config)
    elif gff_type == G.GFF_FORMAT_PRODIGAL:
        return ProdigalProcessor(config)
    else:
        raise ValueError(f"Unsupported GFF type: {gff_type}")
