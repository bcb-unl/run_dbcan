import os
import re
import logging
import pandas as pd
from dbcan.configs.cgcfinder_config import CGCFinderConfig
import dbcan.constants.cgcfinder_constants as C

logger = logging.getLogger(__name__)

# Prefer I/O module to parse GFF (compatible with both cases and function names)
from dbcan.IO.gff import read_gff_to_df as io_read_gff_df


class CGCFinder:
    """CGCFinder"""

    def __init__(self, config: CGCFinderConfig):
        self.config = config
        self._validate()

    # ---- properties (directly read from config) ----
    @property
    def output_dir(self) -> str:
        return self.config.output_dir

    @property
    def filename(self) -> str:
        return self.config.gff_file or os.path.join(self.output_dir, C.CGC_GFF_FILE)

    @property
    def num_null_gene(self) -> int:
        return self.config.num_null_gene

    @property
    def base_pair_distance(self) -> int:
        return self.config.base_pair_distance

    @property
    def use_null_genes(self) -> bool:
        return self.config.use_null_genes

    @property
    def use_distance(self) -> bool:
        return self.config.use_distance

    @property
    def additional_genes(self):
        return self.config.additional_genes or []

    @property
    def additional_logic(self) -> str:
        return self.config.additional_logic

    @property
    def min_core_cazyme(self) -> int:
        return self.config.min_core_cazyme

    @property
    def min_cluster_genes(self) -> int:
        return self.config.min_cluster_genes

    @property
    def additional_min_categories(self) -> int:
        return self.config.additional_min_categories

    @property
    def feature_types(self):
        return self.config.feature_types
    @property
    def extend_mode(self) -> str:
        return self.config.extend_mode

    @property
    def extend_bp(self) -> int:
        return self.config.extend_bp

    @property
    def extend_gene_count(self) -> int:
        return self.config.extend_gene_count

    # ---- parameter validation (single entry) ----
    def _validate(self):
        os.makedirs(self.output_dir, exist_ok=True)
        if self.extend_mode not in ('none', 'bp', 'gene'):
            raise ValueError("extend_mode must be 'none', 'bp' or 'gene'")
        if self.extend_mode == 'bp' and self.extend_bp < 0:
            raise ValueError("extend_bp must be >= 0")
        if self.extend_mode == 'gene' and self.extend_gene_count < 0:
            raise ValueError("extend_gene_count must be >= 0")
        if self.additional_logic not in ('all', 'any'):
            raise ValueError("additional_logic must be 'all' or 'any'")
        for v in (self.min_core_cazyme, self.min_cluster_genes, self.additional_min_categories):
            if v < 0:
                raise ValueError("min_core_cazyme / min_cluster_genes / additional_min_categories must be >= 0")

    def read_gff(self):
        """Read GFF using I/O module; fall back to simple parsing on failure"""
        try:
            if not os.path.exists(self.filename):
                logger.error(f"GFF file not found: {self.filename}")
                return False

            logger.info(f"Reading GFF file: {self.filename}")
            if io_read_gff_df:
                try:
                    self.df = io_read_gff_df(self.filename, columns=C.GFF_COLUMNS)
                except Exception as e:
                    logger.warning(f"Failed to read GFF with io_read_gff_df: {e}")
                    # Use chunked reading for large files (>100MB)
                    file_size_mb = Path(self.filename).stat().st_size / (1024 * 1024)
                    if file_size_mb > 100:
                        logger.info(f"Large GFF file detected ({file_size_mb:.1f}MB), using chunked reading")
                        chunks = []
                        for chunk in pd.read_csv(self.filename, sep='\t', names=C.GFF_COLUMNS, chunksize=100000):
                            chunks.append(chunk)
                        self.df = pd.concat(chunks, ignore_index=True)
                    else:
                        self.df = pd.read_csv(self.filename, sep='\t', names=C.GFF_COLUMNS)
            #print("test0",self.df  )

            # Filter comments and feature types
            # if C.CONTIG_ID_COLUMN in self.df.columns:
            #     self.df = self.df[~self.df[C.CONTIG_ID_COLUMN].astype(str).str.startswith('#')]
            # if 'type' in self.df.columns and self.feature_types:
            #     self.df = self.df[self.df['type'].isin(self.feature_types)]

            # Parse attributes column
            def parse_attributes(attr_str: str) -> dict:
                if not isinstance(attr_str, str) or not attr_str:
                    return {}
                items = {}
                for part in attr_str.split(';'):
                    if not part or '=' not in part:
                        continue
                    k, v = part.split('=', 1)
                    items[k.strip()] = v.strip()
                return items

            attrs = self.df[C.ATTRIBUTES_COLUMN].astype(str).map(parse_attributes)
            self.df[C.CGC_ANNOTATION_COLUMN] = attrs.map(lambda d: d.get(C.CGC_ANNOTATION_ATTR, 'null'))
            self.df[C.PROTEIN_ID_COLUMN] = attrs.map(lambda d: d.get(C.PROTEIN_ID_ATTR, ''))
            self.df = self.df[C.CGC_SELECTED_COLUMNS]

            logger.info(f"Loaded {len(self.df)} records from GFF file after filtering")
            return True
        except Exception as e:
            logger.error(f"Error reading GFF file: {str(e)}")
            import traceback; traceback.print_exc()
            return False

    def mark_signature_genes(self):
        """Mark core/additional/signature genes"""
        try:
            core_pattern = '|'.join(re.escape(s) for s in C.CGC_CORE_SIG_TYPES)
            add_pattern = '|'.join(re.escape(s) for s in (self.additional_genes or [])) if self.additional_genes else ''

            self.df[C.IS_CORE_COLUMN] = self.df[C.CGC_ANNOTATION_COLUMN].str.contains(core_pattern, na=False, regex=True)
            self.df[C.IS_ADDITIONAL_COLUMN] = self.df[C.CGC_ANNOTATION_COLUMN].str.contains(add_pattern, na=False, regex=True) if add_pattern else False
            self.df[C.IS_SIGNATURE_COLUMN] = self.df[C.IS_CORE_COLUMN] | self.df[C.IS_ADDITIONAL_COLUMN]

            logger.info(f"Marked signatures: core={int(self.df[C.IS_CORE_COLUMN].sum())}, "
                        f"additional={int(self.df[C.IS_ADDITIONAL_COLUMN].sum())}")
            return True
        except Exception as e:
            logger.error(f"Error marking signature genes: {str(e)}")
            import traceback; traceback.print_exc()
            return False

    def find_cgc_clusters(self):
        """Identify CGC clusters based on the defined criteria"""
        try:
            if not hasattr(self, 'df') or self.df.empty:
                logger.error("No GFF data loaded or no signature genes marked.")
                return []

            clusters = []
            cgc_id = 1

            logger.info(
                f"Finding CGC clusters using "
                f"{'distance' if self.use_distance else 'no distance'}, "
                f"{'null genes' if self.use_null_genes else 'no null genes'}; "
                f"max null genes: {self.num_null_gene}, bp distance: {self.base_pair_distance if self.use_distance else 'N/A'}; "
                f"extension: {self.extend_mode}"
                f"{(' (' + str(self.extend_bp) + ' bp)' ) if self.extend_mode=='bp' else ( ' (' + str(self.extend_gene_count) + ' genes)' if self.extend_mode=='gene' else '')}"
            )

            for contig, contig_df in self.df.groupby(C.CONTIG_ID_COLUMN):
                contig_df = contig_df.sort_values([C.START_COLUMN, C.END_COLUMN], kind="mergesort")
                sig_indices = contig_df[contig_df[C.IS_SIGNATURE_COLUMN]].index.to_numpy()

                if len(sig_indices) < 2:
                    continue

                starts = contig_df.loc[sig_indices, C.START_COLUMN].to_numpy()
                ends = contig_df.loc[sig_indices, C.END_COLUMN].to_numpy()
                idx_array = contig_df.index.to_numpy()
                pos_map = {label: pos for pos, label in enumerate(idx_array)}

                last_index = None
                start_index = None

                for i, sig_index in enumerate(sig_indices):
                    if last_index is None:
                        start_index = last_index = sig_index
                        continue

                    distance_valid = (starts[i] - ends[i - 1] <= self.base_pair_distance) if self.use_distance else True
                    pos_curr = pos_map[sig_index]
                    pos_last = pos_map[last_index]
                    null_gene_count = max(0, pos_curr - pos_last - 1)
                    null_gene_valid = (null_gene_count <= self.num_null_gene) if self.use_null_genes else True

                    if distance_valid and null_gene_valid:
                        last_index = sig_index
                    else:
                        base_cluster_df = self._get_base_cluster_window(contig_df, start_index, last_index)
                        cluster_df = self._extend_cluster(contig_df, start_index, last_index)
                        if self.validate_cluster(base_cluster_df):
                            clusters.append(self.process_cluster(cluster_df, cgc_id))
                            cgc_id += 1
                        start_index = last_index = sig_index

                if last_index is not None and start_index is not None:
                    base_cluster_df = self._get_base_cluster_window(contig_df, start_index, last_index)
                    cluster_df = self._extend_cluster(contig_df, start_index, last_index)
                    if self.validate_cluster(base_cluster_df):
                        clusters.append(self.process_cluster(cluster_df, cgc_id))
                        cgc_id += 1

            logger.info(f"Found {len(clusters)} CGC clusters")
            return clusters
        except Exception as e:
            logger.error(f"Error finding CGC clusters: {str(e)}")
            import traceback
            traceback.print_exc()
            return []

    # NEW: safe base-window slicer using positional indices
    def _get_base_cluster_window(self, contig_df: pd.DataFrame, start_index, end_index) -> pd.DataFrame:
        idx_array = contig_df.index.to_numpy()
        pos_map = {label: pos for pos, label in enumerate(idx_array)}
        if start_index in pos_map and end_index in pos_map:
            left_pos = pos_map[start_index]
            right_pos = pos_map[end_index]
            if left_pos > right_pos:
                left_pos, right_pos = right_pos, left_pos
            labels = idx_array[left_pos:right_pos + 1]
            return contig_df.loc[labels]
        if (start_index in contig_df.index) and (end_index in contig_df.index):
            try:
                return contig_df.loc[min(start_index, end_index):max(start_index, end_index)]
            except Exception:
                pass
        return contig_df.iloc[0:0]

    # NEW: extend cluster on both sides by bp or gene count
    def _extend_cluster(self, contig_df: pd.DataFrame, start_index, end_index) -> pd.DataFrame:
        idx_array = contig_df.index.to_numpy()
        pos_map = {label: pos for pos, label in enumerate(idx_array)}
        if start_index not in pos_map or end_index not in pos_map:
            return self._get_base_cluster_window(contig_df, start_index, end_index)

        left_pos = pos_map[start_index]
        right_pos = pos_map[end_index]
        if left_pos > right_pos:
            left_pos, right_pos = right_pos, left_pos

        if self.extend_mode == 'none':
            labels = idx_array[left_pos:right_pos + 1]
            return contig_df.loc[labels]

        if self.extend_mode == 'gene':
            k = max(0, int(self.extend_gene_count))
            new_left = max(0, left_pos - k)
            new_right = min(len(idx_array) - 1, right_pos + k)
            labels = idx_array[new_left:new_right + 1]
            return contig_df.loc[labels]

        if self.extend_mode == 'bp':
            bp = max(0, int(self.extend_bp))
            current_labels = idx_array[left_pos:right_pos + 1]
            min_start = int(contig_df.loc[current_labels, C.START_COLUMN].min())
            max_end = int(contig_df.loc[current_labels, C.END_COLUMN].max())
            left_boundary = min_start - bp
            right_boundary = max_end + bp

            new_left = left_pos
            while new_left > 0:
                prev_label = idx_array[new_left - 1]
                prev_end = int(contig_df.at[prev_label, C.END_COLUMN])
                if prev_end >= left_boundary:
                    new_left -= 1
                else:
                    break

            new_right = right_pos
            last_idx = len(idx_array) - 1
            while new_right < last_idx:
                next_label = idx_array[new_right + 1]
                next_start = int(contig_df.at[next_label, C.START_COLUMN])
                if next_start <= right_boundary:
                    new_right += 1
                else:
                    break

            labels = idx_array[new_left:new_right + 1]
            return contig_df.loc[labels]

        labels = idx_array[left_pos:right_pos + 1]
        return contig_df.loc[labels]

    def validate_cluster(self, cluster_df):
        if len(cluster_df) < max(1, self.min_cluster_genes):
            return False
        cazyme_count = int(cluster_df[C.IS_CORE_COLUMN].sum())
        has_core = cazyme_count >= self.min_core_cazyme

        # fix: compare with string "CAZyme" instead of list C.CGC_CORE_SIG_TYPES
        if len(self.additional_genes) == 1 and self.additional_genes[0] == "CAZyme":
            return cazyme_count >= max(2, self.min_core_cazyme)

        additional_annotations = set()
        if self.additional_genes:
            for annotation in cluster_df[cluster_df[C.IS_ADDITIONAL_COLUMN]][C.CGC_ANNOTATION_COLUMN]:
                for gene_type in self.additional_genes:
                    if gene_type in annotation:
                        additional_annotations.add(gene_type)

        if not self.additional_genes:
            return has_core

        if self.additional_logic == 'all':
            has_required = set(self.additional_genes).issubset(additional_annotations)
        else:
            has_required = len(additional_annotations) >= max(1, self.additional_min_categories)

        return has_core and has_required

    @staticmethod
    def get_gene_type(annotation_str):
        PRIORITY = C.PRIORITY
        types = [ann.split('|')[0] for ann in annotation_str.split('+')]
        return sorted(types, key=lambda t: PRIORITY.get(t, 99))[0] if types else C.NULL_GENE_TYPE

    def process_cluster(self, cluster_df, cgc_id):
        return [{
            C.CGC_ID_FIELD: f'CGC{cgc_id}',
            C.GENE_TYPE_FIELD: self.get_gene_type(gene[C.CGC_ANNOTATION_COLUMN]),
            C.CONTIG_ID_COLUMN: gene[C.CONTIG_ID_COLUMN],
            C.CGC_PROTEIN_ID_FIELD: gene[C.PROTEIN_ID_COLUMN],
            C.GENE_START_FIELD: gene[C.START_COLUMN],
            C.GENE_STOP_FIELD: gene[C.END_COLUMN],
            C.GENE_STRAND_FIELD: gene[C.STRAND_COLUMN],
            C.GENE_ANNOTATION_FIELD: gene[C.CGC_ANNOTATION_COLUMN]
        } for _, gene in cluster_df.iterrows()]

    def output_clusters(self, clusters):
        try:
            output_path = os.path.join(self.output_dir, C.CGC_RESULT_FILE)
            if not clusters:
                logger.warning("No CGC clusters found to output")
                pd.DataFrame(columns=C.CGC_OUTPUT_COLUMNS).to_csv(output_path, sep='\t', index=False)

                summary_path = output_path.replace(".tsv", "_summary.tsv")
                pd.DataFrame(columns=["CGC#", "Contig ID", "Cluster Start", "Cluster End", "Genes",
                                     "CAZymes", "TC", "TF", "STP", "Sulfatase", "Peptidase",
                                     "Signatures", "Length (bp)"]).to_csv(summary_path, sep='\t', index=False)
                logger.info(f"Empty CGC output file created at {output_path}")
                return

            # flatten clusters
            rows = []
            for cluster in clusters:
                rows.extend(cluster)
            df_output = pd.DataFrame(rows)

            # write standard output with ONLY the original 8 columns
            df_output.to_csv(output_path, sep='\t', index=False)
            logger.info(f"CGC clusters have been written to {output_path}")

            # build summary with detailed additional gene type columns
            def _summarize(group: pd.DataFrame):
                contig = group[C.CONTIG_ID_COLUMN].iloc[0]
                start = int(group[C.GENE_START_FIELD].min())
                end = int(group[C.GENE_STOP_FIELD].max())
                genes = int(len(group))
                cazy = int((group[C.GENE_TYPE_FIELD] == "CAZyme").sum())

                # Count specific additional gene types
                tc = int((group[C.GENE_TYPE_FIELD] == "TC").sum())
                tf = int((group[C.GENE_TYPE_FIELD] == "TF").sum())
                stp = int((group[C.GENE_TYPE_FIELD] == "STP").sum())
                sulfatase = int((group[C.GENE_TYPE_FIELD] == "Sulfatase").sum())
                peptidase = int((group[C.GENE_TYPE_FIELD] == "Peptidase").sum())

                sigs = int(cazy + tc + tf + stp + sulfatase + peptidase)
                span = end - start + 1
                return pd.Series({
                    "Contig ID": contig,
                    "Cluster Start": start,
                    "Cluster End": end,
                    "Genes": genes,
                    "CAZymes": cazy,
                    "TC": tc,
                    "TF": tf,
                    "STP": stp,
                    "Sulfatase": sulfatase,
                    "Peptidase": peptidase,
                    "Signatures": sigs,
                    "Length (bp)": span
                })

            summary_df = df_output.groupby(C.CGC_ID_FIELD, sort=False).apply(_summarize).reset_index()
            summary_path = output_path.replace(".tsv", "_summary.tsv")
            summary_df.to_csv(summary_path, sep='\t', index=False)
            logger.info(f"CGC summary has been written to {summary_path}")
        except Exception as e:
            logger.error(f"Error outputting CGC clusters: {str(e)}")
            import traceback; traceback.print_exc()

    def run(self):
        if not self.read_gff():
            return False
        if not self.mark_signature_genes():
            return False
        clusters = self.find_cgc_clusters()
        self.output_clusters(clusters)
        logger.info("CGCFinder run completed")
        return True
