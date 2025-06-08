import pandas as pd
import os
import re
import logging
from dbcan.parameter import OverviewGeneratorConfig
from Bio import SeqIO
from dbcan.constants import (OVERVIEW_FILE, NON_CAZYME_FAA_FILE, INPUT_PROTEIN_NAME,
                           DIAMOND_RESULT_FILE, DBCAN_SUB_RESULT_FILE, DBCAN_HMM_RESULT_FILE,
                           OVERVIEW_COLUMNS, DIAMOND_COLUMN_NAMES_OVERVIEW, DBCAN_SUB_COLUMN_NAMES_OVERVIEW,
                           DBCAN_HMM_COLUMN_NAMES_OVERVIEW, GENE_ID_FIELD, EC_FIELD, DBCAN_HMM_FIELD,
                           DBCAN_SUB_FIELD, DIAMOND_FIELD, TOOLS_COUNT_FIELD,
                           RECOMMEND_RESULTS_FIELD, EMPTY_RESULT_PLACEHOLDER,
                           SUBFAMILY_NAME_FIELD, HMM_NAME_FIELD, TARGET_NAME_FIELD,
                           TARGET_FROM_FIELD, TARGET_TO_FIELD, I_EVALUE_FIELD,
                           CAZY_ID_FIELD, SUBFAMILY_EC_FIELD, OVERVIEW_OVERLAP_THRESHOLD,
                           MIN_TOOLS_FOR_RECOMMENDATION, CAZY_ID_PATTERN,
                           RESULT_SEPARATOR, EC_SEPARATOR)

class OverviewGenerator:
    """Generate overview of CAZyme annotations using template method pattern"""

    def __init__(self, config: OverviewGeneratorConfig):
        """Initialize with configuration"""
        self.config = config
        self._setup_processor()

    def _setup_processor(self):
        """Set up processor attributes using template method pattern"""
        self.output_dir = self.config.output_dir
        self.file_paths = self._derive_file_paths()
        self.column_names = self._derive_column_names()
        self.overview_columns = OVERVIEW_COLUMNS
        self.overlap_threshold = OVERVIEW_OVERLAP_THRESHOLD

        # For non-CAZyme FAA generation
        self.cazyme_overview = os.path.join(self.output_dir, OVERVIEW_FILE)
        self.input_total_faa = self._derive_input_total_faa()

        # Validate required attributes
        self._validate_attributes()

    def _validate_attributes(self):
        """Validate that all required attributes are properly set"""
        if not os.path.exists(self.output_dir):
            raise FileNotFoundError(f"Output directory not found: {self.output_dir}")

        # Check that at least one result file exists
        found_files = False
        for key, file_path in self.file_paths.items():
            if os.path.exists(file_path):
                found_files = True
                logging.info(f"Found {key} results at {file_path}")
            else:
                logging.warning(f"{key} results not found at {file_path}")

        if not found_files:
            logging.warning("No CAZyme annotation results found. Overview will be empty.")

    def _derive_file_paths(self):
        """Derive file paths for annotation results"""
        return {
            'diamond': os.path.join(self.output_dir, DIAMOND_RESULT_FILE),
            'dbcan_sub': os.path.join(self.output_dir, DBCAN_SUB_RESULT_FILE),
            'dbcan_hmm': os.path.join(self.output_dir, DBCAN_HMM_RESULT_FILE)
        }

    def _derive_column_names(self):
        """Derive column names for annotation results"""
        return {
            'diamond': DIAMOND_COLUMN_NAMES_OVERVIEW,
            'dbcan_sub': DBCAN_SUB_COLUMN_NAMES_OVERVIEW,
            'dbcan_hmm': DBCAN_HMM_COLUMN_NAMES_OVERVIEW
        }

    def _derive_input_total_faa(self):
        """Derive path to total input FAA file"""
        return os.path.join(self.output_dir, INPUT_PROTEIN_NAME)

    def load_data(self):
        """Load data from annotation result files"""
        data = {}
        for key, file_path in self.file_paths.items():
            if os.path.exists(file_path):
                try:
                    df = pd.read_csv(file_path, sep='\t')

                    # Ensure all required columns exist
                    if not all(col in df.columns for col in self.column_names[key]):
                        logging.warning(f"Missing columns in {file_path}. Expected: {self.column_names[key]}, Found: {df.columns}")
                        continue

                    # Filter to only required columns
                    df = df[self.column_names[key]]

                    # Process specific to each result type
                    if key == 'diamond':
                        df[CAZY_ID_FIELD] = df[CAZY_ID_FIELD].apply(self.extract_cazy_id)
                    elif key in ['dbcan_hmm', 'dbcan_sub']:
                        # Handle HMM name formatting
                        hmm_col = HMM_NAME_FIELD if key == 'dbcan_hmm' else SUBFAMILY_NAME_FIELD
                        df[hmm_col] = df[hmm_col].apply(
                            lambda x: x.split('.hmm')[0] if isinstance(x, str) and '.hmm' in x else x
                        )

                    data[key] = df
                    logging.info(f"Loaded {len(df)} rows from {key} results")
                except Exception as e:
                    logging.error(f"Error loading {key} results: {e}")

        return data

    @staticmethod
    def extract_cazy_id(cazy_id):
        """Extract CAZy ID from DIAMOND result"""
        if not isinstance(cazy_id, str):
            return cazy_id

        parts = cazy_id.split('|')
        for part in parts:
            if re.match(CAZY_ID_PATTERN, part):
                return RESULT_SEPARATOR.join(parts[parts.index(part):])
        return cazy_id

    def calculate_overlap(self, start1, end1, start2, end2):
        """Calculate overlap between two regions"""
        start_max = max(start1, start2)
        end_min = min(end1, end2)
        overlap = max(0, end_min - start_max + 1)
        length1 = end1 - start1 + 1
        length2 = end2 - start2 + 1
        return overlap / min(length1, length2) > self.overlap_threshold

    def select_best_result(self, group):
        """Choose the best result from a group of overlapping annotations"""
        # First priority: HMM results that contain "_" (CAZy subfamily)
        for entry in group:
            if "_" in entry[0] and entry[3] == 'hmm':
                return entry

        # Second priority: subfamily (sub) results
        sub_results = [entry for entry in group if entry[3] == 'sub']
        if sub_results:
            return sub_results[0]

        # Third priority: remaining HMM results
        hmm_results = [entry for entry in group if entry[3] == 'hmm']
        if hmm_results:
            return hmm_results[0]

        # Fallback: first result in group
        return group[0]

    def determine_best_result(self, gene_id, data):
        """Determine best result for a gene"""
        results = {EC_FIELD: EMPTY_RESULT_PLACEHOLDER, DBCAN_HMM_FIELD: EMPTY_RESULT_PLACEHOLDER, DBCAN_SUB_FIELD: EMPTY_RESULT_PLACEHOLDER, DIAMOND_FIELD: EMPTY_RESULT_PLACEHOLDER, TOOLS_COUNT_FIELD: 0, RECOMMEND_RESULTS_FIELD: EMPTY_RESULT_PLACEHOLDER}

        # Process HMMER results
        if 'dbcan_hmm' in data and not data['dbcan_hmm'].empty:
            hmm_results = data['dbcan_hmm'][data['dbcan_hmm'][TARGET_NAME_FIELD] == gene_id]
            if not hmm_results.empty:
                results[DBCAN_HMM_FIELD] = RESULT_SEPARATOR.join([f"{row[HMM_NAME_FIELD]}({row[TARGET_FROM_FIELD]}-{row[TARGET_TO_FIELD]})" for _, row in hmm_results.iterrows()])
                results[TOOLS_COUNT_FIELD] += 1

        # Process dbCAN-sub results
        if 'dbcan_sub' in data and not data['dbcan_sub'].empty:
            sub_results = data['dbcan_sub'][data['dbcan_sub'][TARGET_NAME_FIELD] == gene_id]
            if not sub_results.empty:
                results[DBCAN_SUB_FIELD] = RESULT_SEPARATOR.join([f"{row[SUBFAMILY_NAME_FIELD]}({row[TARGET_FROM_FIELD]}-{row[TARGET_TO_FIELD]})" for _, row in sub_results.iterrows()])
                results[EC_FIELD] = EC_SEPARATOR.join([str(ec) if ec is not None else EMPTY_RESULT_PLACEHOLDER for ec in sub_results[SUBFAMILY_EC_FIELD].fillna(EMPTY_RESULT_PLACEHOLDER).tolist()])
                results[TOOLS_COUNT_FIELD] += 1

        # Process DIAMOND results
        if 'diamond' in data and not data['diamond'].empty:
            diamond_results = data['diamond'][data['diamond'][GENE_ID_FIELD] == gene_id]
            if not diamond_results.empty:
                results[DIAMOND_FIELD] = RESULT_SEPARATOR.join(diamond_results[CAZY_ID_FIELD].tolist())
                results[TOOLS_COUNT_FIELD] += 1

        # Only add Recommend Results if at least 2 tools detected the gene
        if results[TOOLS_COUNT_FIELD] >= MIN_TOOLS_FOR_RECOMMENDATION:

            if results[DBCAN_HMM_FIELD] != EMPTY_RESULT_PLACEHOLDER and results[DBCAN_SUB_FIELD] != EMPTY_RESULT_PLACEHOLDER:
                # combine HMM and subfamily results
                all_results = []
                for _, hr in hmm_results.iterrows():
                    all_results.append((hr[HMM_NAME_FIELD], hr[TARGET_FROM_FIELD], hr[TARGET_TO_FIELD], 'hmm'))
                for _, sr in sub_results.iterrows():
                    all_results.append((sr[SUBFAMILY_NAME_FIELD], sr[TARGET_FROM_FIELD], sr[TARGET_TO_FIELD], 'sub'))

                overlap_results = self.graph_based_grouping(all_results)

                sorted_results = sorted(overlap_results, key=lambda x: x[1])
                domain_with_range = [f"{res[0]}({res[1]}-{res[2]})" for res in sorted_results]
                domain_with_range = list(dict.fromkeys(domain_with_range))
                domain_names = [d.split('(')[0] for d in domain_with_range]
                #domain_names = domain_with_range
                results[RECOMMEND_RESULTS_FIELD] = EC_SEPARATOR.join(domain_names)

                #results[RECOMMEND_RESULTS_FIELD] = EC_SEPARATOR.join([str(res[0]) for res in sorted_results])
            elif results[DBCAN_HMM_FIELD] != EMPTY_RESULT_PLACEHOLDER:
                results[RECOMMEND_RESULTS_FIELD] = EC_SEPARATOR.join([name.split('(')[0] for name in results[DBCAN_HMM_FIELD].split(RESULT_SEPARATOR)])
            elif results[DBCAN_SUB_FIELD] != EMPTY_RESULT_PLACEHOLDER:
                results[RECOMMEND_RESULTS_FIELD] = EC_SEPARATOR.join([name.split('(')[0] for name in results[DBCAN_SUB_FIELD].split(RESULT_SEPARATOR)])

        return results

    def graph_based_grouping(self, all_results):
        """Process annotation results using improved domain-aware grouping algorithm"""
        if not all_results:
            return []

        # First separate results by type and sort by position
        hmm_results = sorted([r for r in all_results if r[3] == 'hmm'], key=lambda x: x[1])
        sub_results = sorted([r for r in all_results if r[3] == 'sub'], key=lambda x: x[1])

        # Check for special case 1: one subfamily spans multiple HMMs
        if len(sub_results) == 1 and len(hmm_results) > 1:
            sub = sub_results[0]
            # Check if this subfamily spans all HMMs
            spans_all = True
            for hmm in hmm_results:
                if not self.calculate_overlap(sub[1], sub[2], hmm[1], hmm[2]):
                    spans_all = False
                    break

            if spans_all:
                # If one subfamily spans all HMMs, prioritize it according to rules
                all_annotations = sub_results + hmm_results
                best = self.select_best_result(all_annotations)
                return [best]

        # Check for special case 2: one HMM spans multiple subfamilies of same name
        if len(hmm_results) == 1 and len(sub_results) > 1:
            hmm = hmm_results[0]
            # Check if all subfamilies have the same name
            if "_" in hmm[0]:
                sub_names = set(sub[0] for sub in sub_results)
                if len(sub_names) == 1:
                    # Check if HMM spans all these subfamilies
                    spans_all = True
                    for sub in sub_results:
                        if not self.calculate_overlap(hmm[1], hmm[2], sub[1], sub[2]):
                            spans_all = False
                            break

                    if spans_all:
                        # If one HMM spans all same-name subfamilies, use selection rules
                        all_annotations = [hmm] + sub_results
                        best = self.select_best_result(all_annotations)
                        return [best]

        # Create domain groups based on position overlap
        domain_groups = []
        processed_hmms = set()
        processed_subs = set()

        # Step 1: Group directly overlapping subfamilies
        sub_groups = []
        for i, sub1 in enumerate(sub_results):
            if i in processed_subs:
                continue

            current_group = [sub1]
            processed_subs.add(i)

            # Use more robust algorithm to find all connected subfamilies
            has_new_addition = True
            while has_new_addition:
                has_new_addition = False
                for j, sub2 in enumerate(sub_results):
                    if j in processed_subs:
                        continue

                    # Check if sub2 overlaps with any subfamily in the current group
                    if any(self.calculate_overlap(group_sub[1], group_sub[2], sub2[1], sub2[2]) for group_sub in current_group):
                        current_group.append(sub2)
                        processed_subs.add(j)
                        has_new_addition = True

            sub_groups.append(current_group)

        # Reset processed_subs for the main processing
        processed_subs = set()

        # Step 2: Create comprehensive groups that account for all overlaps
        # First, create initial groups from each annotation
        all_groups = []

        # Start with subfamilies
        for sub_group in sub_groups:
            group = sub_group.copy()

            # Find all HMMs that overlap with any subfamily in this group
            for i, hmm in enumerate(hmm_results):
                if any(self.calculate_overlap(sub[1], sub[2], hmm[1], hmm[2]) for sub in sub_group):
                    group.append(hmm)
                    processed_hmms.add(i)

            all_groups.append(group)

            # Mark all subfamilies as processed
            for sub in sub_group:
                for j, sub_j in enumerate(sub_results):
                    if sub == sub_j:
                        processed_subs.add(j)
                        break

        # Process any remaining HMMs that don't overlap with any subfamily
        for i, hmm in enumerate(hmm_results):
            if i in processed_hmms:
                continue

            # Check if this HMM overlaps with any other unprocessed HMM
            group = [hmm]
            processed_hmms.add(i)

            for j, hmm2 in enumerate(hmm_results):
                if j in processed_hmms or j == i:
                    continue

                if self.calculate_overlap(hmm[1], hmm[2], hmm2[1], hmm2[2]):
                    group.append(hmm2)
                    processed_hmms.add(j)

            all_groups.append(group)

        # Step 3: Select best result from each group
        final_results = []
        for group in all_groups:
            if not group:
                continue
            # Select best result according to priority rules
            best = self.select_best_result(group)
            final_results.append(best)

        return final_results

    def aggregate_data(self, gene_ids, data):
        """Aggregate data for all genes"""
        aggregated_results = []
        for gene_id in sorted(gene_ids):
            result = self.determine_best_result(gene_id, data)
            aggregated_results.append([gene_id] + list(result.values()))
        return pd.DataFrame(aggregated_results, columns=self.overview_columns)

    def generate_non_cazyme_faa(self):
        """Generate FAA file with non-CAZyme sequences"""
        try:
            # Only generate if input and output can be determined
            if not hasattr(self, 'input_total_faa') or self.input_total_faa is None:
                logging.error("Cannot generate non-CAZyme FAA: input_total_faa not set")
                return

            if not os.path.exists(self.input_total_faa):
                logging.error(f"Cannot generate non-CAZyme FAA: input file not found: {self.input_total_faa}")
                return

            if not os.path.exists(self.cazyme_overview):
                logging.error(f"Cannot generate non-CAZyme FAA: overview file not found: {self.cazyme_overview}")
                return

            # Read overview and get CAZyme IDs with >= 2 tools
            df = pd.read_csv(self.cazyme_overview, sep='\t')
            filtered_df = df[df[TOOLS_COUNT_FIELD] >= MIN_TOOLS_FOR_RECOMMENDATION]
            cazyme_ids = set(filtered_df[GENE_ID_FIELD].tolist())

            # Write non-CAZymes to output
            output_path = os.path.join(self.output_dir, NON_CAZYME_FAA_FILE)
            count = 0

            with open(self.input_total_faa, 'r') as infile, open(output_path, 'w') as outfile:
                for record in SeqIO.parse(infile, 'fasta'):
                    header_id = record.id.split()[0]
                    if header_id not in cazyme_ids:
                        SeqIO.write(record, outfile, 'fasta')
                        count += 1

            logging.info(f"Non-CAZyme FAA file generated with {count} sequences at {output_path}")

        except Exception as e:
            logging.error(f"Failed to generate non-CAZyme FAA: {str(e)}")

    def run(self):
        """Run overview generation"""
        try:
            # Load data from result files
            loaded_data = self.load_data()

            # If no data was loaded, create empty overview
            if not loaded_data:
                logging.warning("No annotation results found. Creating empty overview.")
                empty_df = pd.DataFrame(columns=self.overview_columns)
                output_path = os.path.join(self.output_dir, OVERVIEW_FILE)
                empty_df.to_csv(output_path, sep='\t', index=False)
                print(f"Empty overview saved to: {output_path}")
                return

            # Collect all gene IDs from all datasets
            gene_ids = set()
            for key, dataset in loaded_data.items():
                id_col = TARGET_NAME_FIELD if key in ['dbcan_hmm', 'dbcan_sub'] else GENE_ID_FIELD
                if id_col in dataset.columns:
                    gene_ids.update(dataset[id_col].unique())

            # Aggregate data for all genes
            aggregated_results = self.aggregate_data(gene_ids, loaded_data)

            # Save overview to file
            output_path = os.path.join(self.output_dir, OVERVIEW_FILE)
            aggregated_results.to_csv(output_path, sep='\t', index=False)
            print(f"Aggregated results saved to: {output_path}")

            # Generate non-CAZyme FAA file
            self.generate_non_cazyme_faa()

        except Exception as e:
            logging.error(f"Error generating overview: {str(e)}")
            import traceback
            traceback.print_exc()
