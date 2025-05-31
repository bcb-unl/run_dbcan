import pandas as pd
import numpy as np
import os
import logging
from dbcan.parameter import CGCFinderConfig
from dbcan.constants import (GFF_COLUMNS, CGC_CORE_SIG_TYPES, CGC_DEFAULT_NULL_GENE,
                              CGC_DEFAULT_BP_DISTANCE, CGC_DEFAULT_USE_NULL_GENES,
                              CGC_DEFAULT_USE_DISTANCE, CGC_DEFAULT_ADDITIONAL_GENES,
                              CGC_GFF_FILE, CGC_RESULT_FILE, CGC_OUTPUT_COLUMNS,
                              CGC_ANNOTATION_COLUMN, ATTRIBUTES_COLUMN,
                              PROTEIN_ID_COLUMN, PROTEIN_ID_ATTR, CGC_ANNOTATION_ATTR,
                              CGC_SELECTED_COLUMNS, IS_CORE_COLUMN,
                              IS_ADDITIONAL_COLUMN, IS_SIGNATURE_COLUMN,
                              CONTIG_ID_COLUMN, START_COLUMN, END_COLUMN,
                              STRAND_COLUMN, CGC_ID_FIELD, CGC_PROTEIN_ID_FIELD,
                              GENE_TYPE_FIELD, GENE_START_FIELD,
                              GENE_STOP_FIELD, GENE_STRAND_FIELD,
                              GENE_ANNOTATION_FIELD, NULL_GENE_TYPE)

class CGCFinder:
    """CGCFinder"""

    def __init__(self, config: CGCFinderConfig):
        """Initialize the CGCFinder with configuration."""
        self.config = config
        self._setup_processor()

    def _setup_processor(self):
        """setup the processor with derived attributes"""
        # basic attributes
        self.output_dir = self._derive_output_dir()
        self.filename = self._derive_filename()

        # attributes for CGC cluster identification
        self.num_null_gene = self._derive_num_null_gene()
        self.base_pair_distance = self._derive_base_pair_distance()
        self.use_null_genes = self._derive_use_null_genes()
        self.use_distance = self._derive_use_distance()
        self.additional_genes = self._derive_additional_genes()

        # verify attributes
        self._validate_attributes()

    def _validate_attributes(self):
        """check if required attributes are set and output directory exists"""
        required_attrs = ['output_dir', 'filename']

        for attr in required_attrs:
            if getattr(self, attr, None) is None:
                raise ValueError(f"Required attribute '{attr}' was not properly set")

        # ensure output directory exists
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir, exist_ok=True)
            logging.info(f"Created output directory: {self.output_dir}")

    def _derive_output_dir(self):
        """generate output directory from config"""
        return self.config.output_dir

    def _derive_filename(self):
        """generate filename from config"""
        return os.path.join(self.output_dir, CGC_GFF_FILE)

    def _derive_num_null_gene(self):
        """generate the maximum number of null genes allowed between signature genes"""
        return getattr(self.config, 'num_null_gene', CGC_DEFAULT_NULL_GENE)

    def _derive_base_pair_distance(self):
        """generate the maximum base pair distance between signature genes"""
        return getattr(self.config, 'base_pair_distance', CGC_DEFAULT_BP_DISTANCE)

    def _derive_use_null_genes(self):
        """control whether to consider null genes in the distance calculation"""
        return getattr(self.config, 'use_null_genes', CGC_DEFAULT_USE_NULL_GENES)

    def _derive_use_distance(self):
        """consider distance between signature genes when identifying clusters"""
        return getattr(self.config, 'use_distance', CGC_DEFAULT_USE_DISTANCE)

    def _derive_additional_genes(self):
        """generate additional genes to be considered as signature genes"""
        return getattr(self.config, 'additional_genes', CGC_DEFAULT_ADDITIONAL_GENES)

    def read_gff(self):
        """read GFF file and extract relevant information"""
        try:
            if not os.path.exists(self.filename):
                logging.error(f"GFF file not found: {self.filename}")
                return False

            logging.info(f"Reading GFF file: {self.filename}")
            self.df = pd.read_csv(self.filename, sep='\t', names=GFF_COLUMNS, comment='#')

            # extract relevant columns
            self.df[CGC_ANNOTATION_COLUMN] = self.df[ATTRIBUTES_COLUMN].apply(
                lambda x: dict(item.split('=') for item in x.split(';') if '=' in item).get(CGC_ANNOTATION_ATTR, '')
            )
            self.df[PROTEIN_ID_COLUMN] = self.df[ATTRIBUTES_COLUMN].apply(
                lambda x: dict(item.split('=') for item in x.split(';') if '=' in item).get(PROTEIN_ID_ATTR, '')
            )
            self.df = self.df[CGC_SELECTED_COLUMNS]
            logging.info(f"Loaded {len(self.df)} records from GFF file")
            return True
        except Exception as e:
            logging.error(f"Error reading GFF file: {str(e)}")
            import traceback
            traceback.print_exc()
            return False

    def mark_signature_genes(self):
        """annotate signature genes in the dataframe"""
        try:
            if not hasattr(self, 'df') or self.df.empty:
                logging.error("No GFF data loaded. Run read_gff() first.")
                return False

            core_sig_types = CGC_CORE_SIG_TYPES
            self.df[IS_CORE_COLUMN] = self.df[CGC_ANNOTATION_COLUMN].str.contains('|'.join(core_sig_types), na=False)
            self.df[IS_ADDITIONAL_COLUMN] = self.df[CGC_ANNOTATION_COLUMN].str.contains('|'.join(self.additional_genes), na=False)
            self.df[IS_SIGNATURE_COLUMN] = self.df[IS_CORE_COLUMN] | self.df[IS_ADDITIONAL_COLUMN]

            sig_gene_count = self.df[IS_SIGNATURE_COLUMN].sum()
            logging.info(f"Marked {sig_gene_count} signature genes ({self.df[IS_CORE_COLUMN].sum()} core, {self.df[IS_ADDITIONAL_COLUMN].sum()} additional)")
            return True
        except Exception as e:
            logging.error(f"Error marking signature genes: {str(e)}")
            import traceback
            traceback.print_exc()
            return False

    def find_cgc_clusters(self):
        """identify CGC clusters based on the defined criteria"""
        try:
            if not hasattr(self, 'df') or self.df.empty:
                logging.error("No GFF data loaded or no signature genes marked.")
                return []

            clusters = []
            cgc_id = 1

            logging.info(
                f"Finding CGC clusters using "
                f"{'distance' if self.use_distance else 'no distance'}, "
                f"{'null genes' if self.use_null_genes else 'no null genes'}; "
                f"max null genes: {self.num_null_gene}, bp distance: {self.base_pair_distance if self.use_distance else 'N/A'} "
            )

            for contig, contig_df in self.df.groupby(CONTIG_ID_COLUMN):
                sig_indices = contig_df[contig_df[IS_SIGNATURE_COLUMN]].index.to_numpy()

                if len(sig_indices) < 2:
                    continue  # need at least 2 signature genes to form a cluster

                starts = contig_df.loc[sig_indices, START_COLUMN].to_numpy()
                ends = contig_df.loc[sig_indices, END_COLUMN].to_numpy()

                last_index = None
                start_index = None

                for i, sig_index in enumerate(sig_indices):
                    if last_index is None:
                        start_index = last_index = sig_index
                        continue

                    distance_valid = (starts[i] - ends[i - 1] <= self.base_pair_distance) if self.use_distance else True
                    null_gene_count = (sig_index - last_index - 1)
                    null_gene_valid = (null_gene_count <= self.num_null_gene) if self.use_null_genes else True

                    if distance_valid and null_gene_valid:
                        last_index = sig_index
                    else:
                        cluster_df = contig_df.loc[start_index:last_index]
                        if self.validate_cluster(cluster_df):
                            clusters.append(self.process_cluster(cluster_df, cgc_id))
                            cgc_id += 1
                        start_index = last_index = sig_index

                # process the last cluster if it exists
                if last_index is not None and start_index is not None:
                    cluster_df = contig_df.loc[start_index:last_index]
                    if self.validate_cluster(cluster_df):
                        clusters.append(self.process_cluster(cluster_df, cgc_id))
                        cgc_id += 1

            logging.info(f"Found {len(clusters)} CGC clusters")
            return clusters
        except Exception as e:
            logging.error(f"Error finding CGC clusters: {str(e)}")
            import traceback
            traceback.print_exc()
            return []

    def validate_cluster(self, cluster_df):
        """validate if a cluster meets the criteria for being a CGC"""
        if len(cluster_df) < 2:
            return False
        has_core = cluster_df[IS_CORE_COLUMN].any()

        # special case: CAZyme-only mode, need at least 2 CAZymes
        if len(self.additional_genes) == 1 and self.additional_genes[0] == "CAZyme":
            cazyme_count = cluster_df[IS_CORE_COLUMN].sum()
            return has_core and cazyme_count >= 2

        # nomal case: check if all additional genes are present
        additional_annotations = set()
        for annotation in cluster_df[cluster_df[IS_ADDITIONAL_COLUMN]][CGC_ANNOTATION_COLUMN]:
            for gene_type in self.additional_genes:
                if gene_type in annotation:
                    additional_annotations.add(gene_type)

        has_all_additional = set(self.additional_genes).issubset(additional_annotations)
        return (has_core and has_all_additional)


    @staticmethod
    def get_gene_type(annotation_str):
        PRIORITY = {'CAZyme': 0, 'TC': 1, 'TF': 2, 'STP': 3, 'SULFATLAS':4, 'PEPTIDASE':5}
        types = [ann.split('|')[0] for ann in annotation_str.split('+')]
        return sorted(types, key=lambda t: PRIORITY.get(t, 99))[0] if types else NULL_GENE_TYPE

    def process_cluster(self, cluster_df, cgc_id):
        """format a cluster for output"""
        return [{
            CGC_ID_FIELD: f'CGC{cgc_id}',
            GENE_TYPE_FIELD: self.get_gene_type(gene[CGC_ANNOTATION_COLUMN]),
            CONTIG_ID_COLUMN: gene[CONTIG_ID_COLUMN],
            CGC_PROTEIN_ID_FIELD: gene[PROTEIN_ID_COLUMN],
            GENE_START_FIELD: gene[START_COLUMN],
            GENE_STOP_FIELD: gene[END_COLUMN],
            GENE_STRAND_FIELD: gene[STRAND_COLUMN],
            GENE_ANNOTATION_FIELD: gene[CGC_ANNOTATION_COLUMN]
        } for _, gene in cluster_df.iterrows()]

    def output_clusters(self, clusters):
        """export identified CGC clusters to a TSV file"""
        try:
            if not clusters:
                logging.warning("No CGC clusters found to output")
                # generate empty file
                empty_df = pd.DataFrame(columns=CGC_OUTPUT_COLUMNS)
                output_path = os.path.join(self.output_dir, CGC_RESULT_FILE)
                empty_df.to_csv(output_path, sep='\t', index=False)
                logging.info(f"Empty CGC output file created at {output_path}")
                return

            rows = []
            for cluster in clusters:
                rows.extend(cluster)

            df_output = pd.DataFrame(rows)
            output_path = os.path.join(self.output_dir, CGC_RESULT_FILE)
            df_output.to_csv(output_path, sep='\t', index=False)
            logging.info(f"CGC clusters have been written to {output_path}")
        except Exception as e:
            logging.error(f"Error outputting CGC clusters: {str(e)}")
            import traceback
            traceback.print_exc()

    def run(self):
        """run the CGCFinder"""
        if not self.read_gff():
            return False
        if not self.mark_signature_genes():
            return False
        clusters = self.find_cgc_clusters()
        self.output_clusters(clusters)
        logging.info("CGCFinder run completed")
        return True
