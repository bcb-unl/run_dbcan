import logging
import re
from pathlib import Path
from typing import Dict, List, Iterable, Tuple, Any

import pandas as pd
from Bio import SeqIO

from dbcan.configs.base_config import OverviewGeneratorConfig
import dbcan.constants.OverviewGenerator_constants as O
# from dbcan.constants import (
#     OVERVIEW_FILE, INPUT_PROTEIN_NAME,
#     DIAMOND_RESULT_FILE, DBCAN_SUB_RESULT_FILE, DBCAN_HMM_RESULT_FILE,
#     OVERVIEW_COLUMNS, DIAMOND_COLUMN_NAMES_OVERVIEW, DBCAN_SUB_COLUMN_NAMES_OVERVIEW,
#     DBCAN_HMM_COLUMN_NAMES_OVERVIEW, GENE_ID_FIELD, EC_FIELD, DBCAN_HMM_FIELD,
#     DBCAN_SUB_FIELD, DIAMOND_FIELD, TOOLS_COUNT_FIELD,
#     RECOMMEND_RESULTS_FIELD, EMPTY_RESULT_PLACEHOLDER,
#     SUBFAMILY_NAME_FIELD, HMM_NAME_FIELD, TARGET_NAME_FIELD,
#     TARGET_FROM_FIELD, TARGET_TO_FIELD, I_EVALUE_FIELD,
#     CAZY_ID_FIELD, SUBFAMILY_EC_FIELD, OVERVIEW_OVERLAP_THRESHOLD,
#     MIN_TOOLS_FOR_RECOMMENDATION, CAZY_ID_PATTERN,
#     RESULT_SEPARATOR, EC_SEPARATOR, DBCAN_SUB_SUBSTRATE_COLUMN
# )

logger = logging.getLogger(__name__)

class OverviewGenerator:
    """Generate overview of CAZyme annotations."""

    def __init__(self, config: OverviewGeneratorConfig):
        self.config = config
        self._validate_basic()

    # ---------------------------
    # Properties (single source)
    # ---------------------------
    @property
    def output_dir_path(self) -> Path:
        return Path(self.config.output_dir)

    @property
    def diamond_result_path(self) -> Path:
        return self.output_dir_path / O.DIAMOND_RESULT_FILE

    @property
    def dbcan_sub_result_path(self) -> Path:
        return self.output_dir_path / O.DBCAN_SUB_RESULT_FILE

    @property
    def dbcan_hmm_result_path(self) -> Path:
        return self.output_dir_path / O.DBCAN_HMM_RESULT_FILE

    @property
    def overview_output_path(self) -> Path:
        return self.output_dir_path / O.OVERVIEW_FILE

    @property
    def input_total_faa_path(self) -> Path:
        return self.output_dir_path / O.INPUT_PROTEIN_NAME

    @property
    def result_file_map(self) -> Dict[str, Path]:
        return {
            "diamond": self.diamond_result_path,
            "dbcan_sub": self.dbcan_sub_result_path,
            "dbcan_hmm": self.dbcan_hmm_result_path,
        }

    @property
    def column_names_map(self) -> Dict[str, List[str]]:
        cols = {
            "diamond": list(O.DIAMOND_COLUMN_NAMES_OVERVIEW),
            "dbcan_sub": list(O.DBCAN_SUB_COLUMN_NAMES_OVERVIEW),
            "dbcan_hmm": list(O.DBCAN_HMM_COLUMN_NAMES_OVERVIEW),
        }
        # Ensure substrate column exists for dbcan_sub
        if O.DBCAN_SUB_SUBSTRATE_COLUMN not in cols["dbcan_sub"]:
            cols["dbcan_sub"].append(O.DBCAN_SUB_SUBSTRATE_COLUMN)
        return cols

    @property
    def overview_columns(self) -> List[str]:
        cols = list(O.OVERVIEW_COLUMNS)
        if O.SUBSTRATE_FIELD not in cols:
            cols.append(O.SUBSTRATE_FIELD)
        return cols

    @property
    def overlap_threshold(self) -> float:
        return O.OVERVIEW_OVERLAP_THRESHOLD

    # ---------------------------
    # Validation
    # ---------------------------
    def _validate_basic(self):
        self.output_dir_path.mkdir(parents=True, exist_ok=True)

    # ---------------------------
    # Data loading
    # ---------------------------
    def load_data(self) -> Dict[str, pd.DataFrame]:
        data: Dict[str, pd.DataFrame] = {}
        for key, path in self.result_file_map.items():
            if not path.exists():
                logger.warning(f"{key} results not found at {path}")
                continue
            try:
                # Use chunked reading for large files (>100MB)
                file_size_mb = path.stat().st_size / (1024 * 1024)
                if file_size_mb > 100:
                    logger.info(f"Large file detected ({file_size_mb:.1f}MB), using chunked reading for {key}")
                    chunks = []
                    for chunk in pd.read_csv(path, sep='\t', chunksize=100000):
                        chunks.append(chunk)
                    df = pd.concat(chunks, ignore_index=True)
                else:
                    df = pd.read_csv(path, sep='\t')
                
                required = list(self.column_names_map[key])
                missing = [c for c in required if c not in df.columns]
                if missing:
                    # Allow auto-add substrate for dbcan_sub
                    if key == "dbcan_sub" and O.DBCAN_SUB_SUBSTRATE_COLUMN in missing:
                        df[O.DBCAN_SUB_SUBSTRATE_COLUMN] = "-"
                        missing = [m for m in missing if m != O.DBCAN_SUB_SUBSTRATE_COLUMN]
                    if missing:
                        logger.warning(
                            f"Missing columns in {path}. Expected subset: {required}. "
                            f"Found: {df.columns.tolist()}. Missing: {missing}. Skipping dataset."
                        )
                        continue

                df = df[self.column_names_map[key]]

                if key == "diamond":
                    df[O.CAZY_ID_FIELD] = df[O.CAZY_ID_FIELD].apply(self.extract_cazy_id)
                elif key in ("dbcan_hmm", "dbcan_sub"):
                    hmm_col = O.HMM_NAME_FIELD if key == "dbcan_hmm" else O.SUBFAMILY_NAME_FIELD
                    df[hmm_col] = df[hmm_col].apply(
                        lambda x: x.split(".hmm")[0] if isinstance(x, str) and ".hmm" in x else x
                    )
                    if key == "dbcan_sub" and O.DBCAN_SUB_SUBSTRATE_COLUMN in df.columns:
                        df[O.DBCAN_SUB_SUBSTRATE_COLUMN] = (
                            df[O.DBCAN_SUB_SUBSTRATE_COLUMN].fillna("-").astype(str)
                        )

                data[key] = df
                logger.info(f"Loaded {len(df)} rows from {key} results")
            except Exception as e:
                logger.error(f"Error loading {key} results: {e}")
        return data

    # ---------------------------
    # Helper / extraction logic
    # ---------------------------
    @staticmethod
    def extract_cazy_id(cazy_id):
        """
        Extract CAZy family ids from a '|' separated string.
        Rules:
          - If any token contains 'fasta' (case-insensitive), truncate at the first pure-digit token
            (drop that token and everything after it).
          - Find the first token matching O.CAZY_ID_PATTERN; return all matching tokens from that
            position onward, joined by O.RESULT_SEPARATOR (keeps multi-domain like GH5|GH3|GH6).
          - If truncation happened but no CAZy token found, return the first remaining token.
          - Otherwise, return the original value.
        """
        if not isinstance(cazy_id, str) or not cazy_id:
            return cazy_id

        parts = [p.strip() for p in cazy_id.split('|')]
        truncated = False

        # Special handling when filename present
        if any("fasta" in p.lower() for p in parts if isinstance(p, str)):
            for idx, tok in enumerate(parts):
                if tok.isdigit():
                    parts = parts[:idx]
                    truncated = True
                    break

        # Find first CAZy token index
        first_idx = None
        for i, t in enumerate(parts):
            if re.match(O.CAZY_ID_PATTERN, t or ""):
                first_idx = i
                break

        if first_idx is not None:
            # Collect all CAZy tokens from first match to the end (preserve multi-domain)
            matches = [t for t in parts[first_idx:] if re.match(O.CAZY_ID_PATTERN, t or "")]
            if matches:
                return O.RESULT_SEPARATOR.join(matches)

        # Fallbacks
        if truncated and parts:
            return parts[0]
        return cazy_id

    def calculate_overlap(self, start1: int, end1: int, start2: int, end2: int) -> bool:
        start_max = max(start1, start2)
        end_min = min(end1, end2)
        overlap = max(0, end_min - start_max + 1)
        length1 = end1 - start1 + 1
        length2 = end2 - start2 + 1
        return overlap / min(length1, length2) > self.overlap_threshold

    def select_best_result(self, group: List[Tuple]) -> Tuple:
        # Priority 1: HMM results containing "_" (subfamily pattern in HMM result)
        for entry in group:
            if "_" in entry[0] and entry[3] == 'hmm':
                return entry
        # Priority 2: subfamily results
        subs = [e for e in group if e[3] == 'sub']
        if subs:
            return subs[0]
        # Priority 3: remaining HMM
        hmms = [e for e in group if e[3] == 'hmm']
        if hmms:
            return hmms[0]
        return group[0]

    def determine_best_result(self, gene_id: str, data: Dict[str, pd.DataFrame]) -> Dict[str, Any]:
        results = {
            O.EC_FIELD: O.EMPTY_RESULT_PLACEHOLDER,
            O.DBCAN_HMM_FIELD: O.EMPTY_RESULT_PLACEHOLDER,
            O.DBCAN_SUB_FIELD: O.EMPTY_RESULT_PLACEHOLDER,
            O.DIAMOND_FIELD: O.EMPTY_RESULT_PLACEHOLDER,
            O.TOOLS_COUNT_FIELD: 0,
            O.RECOMMEND_RESULTS_FIELD: O.EMPTY_RESULT_PLACEHOLDER,
            O.SUBSTRATE_FIELD: O.EMPTY_RESULT_PLACEHOLDER,
        }

        # HMM results
        hmm_results = pd.DataFrame()
        if "dbcan_hmm" in data and not data["dbcan_hmm"].empty:
            hmm_results = data["dbcan_hmm"][data["dbcan_hmm"][O.TARGET_NAME_FIELD] == gene_id]
            if not hmm_results.empty:
                results[O.DBCAN_HMM_FIELD] = O.RESULT_SEPARATOR.join(
                    [
                        f"{row[O.HMM_NAME_FIELD]}({row[O.TARGET_FROM_FIELD]}-{row[O.TARGET_TO_FIELD]})"
                        for _, row in hmm_results.iterrows()
                    ]
                )
                results[O.TOOLS_COUNT_FIELD] += 1

        # Subfamily results
        sub_results = pd.DataFrame()
        if "dbcan_sub" in data and not data["dbcan_sub"].empty:
            sub_results = data["dbcan_sub"][data["dbcan_sub"][O.TARGET_NAME_FIELD] == gene_id]
            if not sub_results.empty:
                results[O.DBCAN_SUB_FIELD] = O.RESULT_SEPARATOR.join(
                    [
                        f"{row[O.SUBFAMILY_NAME_FIELD]}({row[O.TARGET_FROM_FIELD]}-{row[O.TARGET_TO_FIELD]})"
                        for _, row in sub_results.iterrows()
                    ]
                )
                results[O.EC_FIELD] = O.EC_SEPARATOR.join(
                    [
                        str(ec) if ec is not None else O.EMPTY_RESULT_PLACEHOLDER
                        for ec in sub_results[O.SUBFAMILY_EC_FIELD]
                        .fillna(O.EMPTY_RESULT_PLACEHOLDER)
                        .tolist()
                    ]
                )
                results[O.TOOLS_COUNT_FIELD] += 1
                # Substrate aggregation
                if O.DBCAN_SUB_SUBSTRATE_COLUMN in sub_results.columns:
                    subs_raw = [
                        s
                        for s in sub_results[O.DBCAN_SUB_SUBSTRATE_COLUMN]
                        .astype(str)
                        .tolist()
                        if s and s != '-' and s.lower() != 'nan'
                    ]
                    if subs_raw:
                        flat: List[str] = []
                        for s in subs_raw:
                            flat.extend([p for p in re.split(r'[;,]', s) if p])
                        uniq = sorted(
                            {p.strip() for p in flat if p and p.strip()}
                        )
                        if uniq:
                            results[O.SUBSTRATE_FIELD] = O.SUB_SEPARATOR.join(uniq)

        # DIAMOND results
        if "diamond" in data and not data["diamond"].empty:
            diamond_results = data["diamond"][data["diamond"][O.GENE_ID_FIELD] == gene_id]
            if not diamond_results.empty:
                results[O.DIAMOND_FIELD] = O.RESULT_SEPARATOR.join(
                    diamond_results[O.CAZY_ID_FIELD].tolist()
                )
                results[O.TOOLS_COUNT_FIELD] += 1

        # Recommend results if >=2 tools
        if results[O.TOOLS_COUNT_FIELD] >= O.MIN_TOOLS_FOR_RECOMMENDATION:
            if (
                results[O.DBCAN_HMM_FIELD] != O.EMPTY_RESULT_PLACEHOLDER
                and results[O.DBCAN_SUB_FIELD] != O.EMPTY_RESULT_PLACEHOLDER
                and not hmm_results.empty
                and not sub_results.empty
            ):
                all_ann = []
                for _, hr in hmm_results.iterrows():
                    all_ann.append(
                        (hr[O.HMM_NAME_FIELD], hr[O.TARGET_FROM_FIELD], hr[O.TARGET_TO_FIELD], 'hmm')
                    )
                for _, sr in sub_results.iterrows():
                    all_ann.append(
                        (
                            sr[O.SUBFAMILY_NAME_FIELD],
                            sr[O.TARGET_FROM_FIELD],
                            sr[O.TARGET_TO_FIELD],
                            'sub',
                        )
                    )
                grouped = self.graph_based_grouping(all_ann)
                sorted_results = sorted(grouped, key=lambda x: x[1])
                domain_with_range = [
                    f"{res[0]}({res[1]}-{res[2]})" for res in sorted_results
                ]
                domain_with_range = list(dict.fromkeys(domain_with_range))
                domain_names = [d.split('(')[0] for d in domain_with_range]
                results[O.RECOMMEND_RESULTS_FIELD] = O.EC_SEPARATOR.join(domain_names)
            elif results[O.DBCAN_HMM_FIELD] != O.EMPTY_RESULT_PLACEHOLDER:
                results[O.RECOMMEND_RESULTS_FIELD] = O.EC_SEPARATOR.join(
                    [name.split('(')[0] for name in results[O.DBCAN_HMM_FIELD].split(O.RESULT_SEPARATOR)]
                )
            elif results[O.DBCAN_SUB_FIELD] != O.EMPTY_RESULT_PLACEHOLDER:
                results[O.RECOMMEND_RESULTS_FIELD] = O.EC_SEPARATOR.join(
                    [name.split('(')[0] for name in results[O.DBCAN_SUB_FIELD].split(O.RESULT_SEPARATOR)]
                )

        return results

    def graph_based_grouping(self, all_results: List[Tuple]) -> List[Tuple]:
        if not all_results:
            return []

        hmm_results = sorted([r for r in all_results if r[3] == 'hmm'], key=lambda x: x[1])
        sub_results = sorted([r for r in all_results if r[3] == 'sub'], key=lambda x: x[1])

        # Special case 1
        if len(sub_results) == 1 and len(hmm_results) > 1:
            sub = sub_results[0]
            if all(self.calculate_overlap(sub[1], sub[2], hmm[1], hmm[2]) for hmm in hmm_results):
                best = self.select_best_result(sub_results + hmm_results)
                return [best]

        # Special case 2
        if len(hmm_results) == 1 and len(sub_results) > 1:
            hmm = hmm_results[0]
            if "_" in hmm[0]:
                sub_names = {s[0] for s in sub_results}
                if len(sub_names) == 1:
                    if all(self.calculate_overlap(hmm[1], hmm[2], sub[1], sub[2]) for sub in sub_results):
                        best = self.select_best_result([hmm] + sub_results)
                        return [best]

        # Group subfamilies by overlap
        processed_subs = set()
        sub_groups: List[List[Tuple]] = []
        for i, sub1 in enumerate(sub_results):
            if i in processed_subs:
                continue
            group = [sub1]
            processed_subs.add(i)
            changed = True
            while changed:
                changed = False
                for j, sub2 in enumerate(sub_results):
                    if j in processed_subs:
                        continue
                    if any(self.calculate_overlap(g[1], g[2], sub2[1], sub2[2]) for g in group):
                        group.append(sub2)
                        processed_subs.add(j)
                        changed = True
            sub_groups.append(group)

        processed_hmms = set()
        all_groups: List[List[Tuple]] = []

        for sub_group in sub_groups:
            group = sub_group.copy()
            for i, hmm in enumerate(hmm_results):
                if any(self.calculate_overlap(sub[1], sub[2], hmm[1], hmm[2]) for sub in sub_group):
                    group.append(hmm)
                    processed_hmms.add(i)
            all_groups.append(group)

        for i, hmm in enumerate(hmm_results):
            if i in processed_hmms:
                continue
            group = [hmm]
            processed_hmms.add(i)
            for j, hmm2 in enumerate(hmm_results):
                if j in processed_hmms or j == i:
                    continue
                if self.calculate_overlap(hmm[1], hmm[2], hmm2[1], hmm2[2]):
                    group.append(hmm2)
                    processed_hmms.add(j)
            all_groups.append(group)

        final_results = []
        for group in all_groups:
            if group:
                final_results.append(self.select_best_result(group))
        return final_results

    def aggregate_data(self, gene_ids: Iterable[str], data: Dict[str, pd.DataFrame]) -> pd.DataFrame:
        aggregated = []
        for gene_id in sorted(gene_ids):
            result = self.determine_best_result(gene_id, data)
            row = []
            for col in self.overview_columns:
                if col == self.overview_columns[0]:
                    row.append(gene_id)
                else:
                    row.append(result.get(col, O.EMPTY_RESULT_PLACEHOLDER))
            aggregated.append(row)
        return pd.DataFrame(aggregated, columns=self.overview_columns)

    # ---------------------------
    # Orchestration
    # ---------------------------
    def run(self):
        try:
            loaded_data = self.load_data()
            if not loaded_data:
                logger.warning("No annotation result files found. Creating empty overview.")
                empty_df = pd.DataFrame(columns=self.overview_columns)
                empty_df.to_csv(self.overview_output_path, sep='\t', index=False)
                logger.info(f"Empty overview saved to: {self.overview_output_path}")
                return

            gene_ids = set()
            for key, df in loaded_data.items():
                id_col = O.TARGET_NAME_FIELD if key in ("dbcan_hmm", "dbcan_sub") else O.GENE_ID_FIELD
                if id_col in df.columns:
                    gene_ids.update(df[id_col].unique())

            aggregated_df = self.aggregate_data(gene_ids, loaded_data)
            aggregated_df.to_csv(self.overview_output_path, sep='\t', index=False)
            logger.info(f"Aggregated overview saved to: {self.overview_output_path}")

        except Exception as e:
            logger.error(f"Error generating overview: {e}", exc_info=True)
