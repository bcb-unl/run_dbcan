import logging
import re
from pathlib import Path
from typing import Dict, List, Iterable, Tuple, Any

import pandas as pd
from Bio import SeqIO

from dbcan.configs.base_config import OverviewGeneratorConfig
from dbcan.process.domain_resolution import (
    intervals_length,
    intervals_overlap,
    parse_intervals,
)
from dbcan.process.process_utils import alignment_details_path
import dbcan.constants.OverviewGenerator_constants as O
import dbcan.constants.process_utils_constants as PU
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
        # (source, target, ali_from, ali_to) -> match-state intervals, loaded
        # from the optional alignment-details sidecars written next to the HMM
        # result tables. Empty when a run used legacy ali-span resolution.
        self._matched_intervals: Dict[Tuple[str, str, int, int], Tuple[Tuple[int, int], ...]] = {}
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
    def alignment_details_map(self) -> Dict[str, Path]:
        """Optional match-state sidecars for the two HMM-based sources."""
        return {
            "dbcan_hmm": alignment_details_path(self.dbcan_hmm_result_path),
            "dbcan_sub": alignment_details_path(self.dbcan_sub_result_path),
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
    def load_alignment_details(self) -> None:
        """Load the resolved-domain sidecars, if the run produced them.

        The public result tables keep the historical column schema, so the
        match-state evidence used to resolve the architecture is carried in a
        separate ``*_alignment_details.tsv`` file. Consuming it here means the
        overview works from already-resolved domains instead of re-inferring the
        architecture from the continuous alignment span. When the sidecar is
        absent (legacy ali-span run, or a hand-assembled output directory) the
        overview transparently falls back to span overlap.
        """
        self._matched_intervals = {}
        for source, path in self.alignment_details_map.items():
            if not path.exists():
                continue
            try:
                df = pd.read_csv(path, sep='\t')
            except Exception as e:
                logger.warning(f"Could not read alignment details {path}: {e}")
                continue
            required = {
                PU.TARGET_NAME_COLUMN,
                PU.TARGET_FROM_COLUMN,
                PU.TARGET_TO_COLUMN,
                PU.MATCHED_INTERVALS_COLUMN,
            }
            if not required.issubset(df.columns):
                logger.warning(f"Alignment details {path} missing columns {required - set(df.columns)}")
                continue
            for _, row in df.iterrows():
                try:
                    ali_from = int(row[PU.TARGET_FROM_COLUMN])
                    ali_to = int(row[PU.TARGET_TO_COLUMN])
                except (TypeError, ValueError):
                    continue
                key = (source, str(row[PU.TARGET_NAME_COLUMN]), ali_from, ali_to)
                # Two resolved hits can only share an identical span if a
                # coexistence policy allowed them, in which case they describe
                # the same region anyway; first one wins.
                self._matched_intervals.setdefault(
                    key,
                    parse_intervals(row[PU.MATCHED_INTERVALS_COLUMN], ali_from, ali_to),
                )
            logger.info(f"Loaded {len(df)} {source} alignment details from {path}")

    def _lookup_intervals(self, source: str, target: str, ali_from, ali_to):
        try:
            key = (source, str(target), int(ali_from), int(ali_to))
        except (TypeError, ValueError):
            return None
        return self._matched_intervals.get(key)

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
        """Alignment-span overlap: the fallback when no match-state data exists."""
        start_max = max(start1, start2)
        end_min = min(end1, end2)
        overlap = max(0, end_min - start_max + 1)
        length1 = end1 - start1 + 1
        length2 = end2 - start2 + 1
        return overlap / min(length1, length2) > self.overlap_threshold

    def annotations_overlap(self, ann1: Tuple, ann2: Tuple) -> bool:
        """Do two annotations cover the same conserved region?

        Uses match-state intervals when both annotations have them, so a
        catalytic domain does not appear to cover an accessory domain that only
        sits inside its insert-state region. Falls back to the alignment span
        otherwise.
        """
        iv1 = ann1[4] if len(ann1) > 4 else None
        iv2 = ann2[4] if len(ann2) > 4 else None
        if iv1 and iv2:
            overlap = intervals_overlap(iv1, iv2)
            len1 = max(1, intervals_length(iv1))
            len2 = max(1, intervals_length(iv2))
            return overlap / min(len1, len2) > self.overlap_threshold
        return self.calculate_overlap(ann1[1], ann1[2], ann2[1], ann2[2])

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

    def determine_best_result(self, gene_id: str, indexed: Dict[str, Dict[str, pd.DataFrame]]) -> Dict[str, Any]:
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
        hmm_results = indexed.get("dbcan_hmm", {}).get(gene_id)
        if hmm_results is not None and not hmm_results.empty:
            results[O.DBCAN_HMM_FIELD] = O.RESULT_SEPARATOR.join(
                [
                    f"{row[O.HMM_NAME_FIELD]}({row[O.TARGET_FROM_FIELD]}-{row[O.TARGET_TO_FIELD]})"
                    for _, row in hmm_results.iterrows()
                ]
            )
            results[O.TOOLS_COUNT_FIELD] += 1

        # Subfamily results
        sub_results = indexed.get("dbcan_sub", {}).get(gene_id)
        if sub_results is not None and not sub_results.empty:
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
        diamond_results = indexed.get("diamond", {}).get(gene_id)
        if diamond_results is not None and not diamond_results.empty:
            results[O.DIAMOND_FIELD] = O.RESULT_SEPARATOR.join(
                diamond_results[O.CAZY_ID_FIELD].tolist()
            )
            results[O.TOOLS_COUNT_FIELD] += 1

        # Recommend results if >=2 tools
        if results[O.TOOLS_COUNT_FIELD] >= O.MIN_TOOLS_FOR_RECOMMENDATION:
            if (
                results[O.DBCAN_HMM_FIELD] != O.EMPTY_RESULT_PLACEHOLDER
                and results[O.DBCAN_SUB_FIELD] != O.EMPTY_RESULT_PLACEHOLDER
                and hmm_results is not None
                and sub_results is not None
                and not hmm_results.empty
                and not sub_results.empty
            ):
                all_ann = []
                for _, hr in hmm_results.iterrows():
                    all_ann.append(
                        (
                            hr[O.HMM_NAME_FIELD],
                            hr[O.TARGET_FROM_FIELD],
                            hr[O.TARGET_TO_FIELD],
                            'hmm',
                            self._lookup_intervals(
                                "dbcan_hmm", gene_id, hr[O.TARGET_FROM_FIELD], hr[O.TARGET_TO_FIELD]
                            ),
                        )
                    )
                for _, sr in sub_results.iterrows():
                    all_ann.append(
                        (
                            sr[O.SUBFAMILY_NAME_FIELD],
                            sr[O.TARGET_FROM_FIELD],
                            sr[O.TARGET_TO_FIELD],
                            'sub',
                            self._lookup_intervals(
                                "dbcan_sub", gene_id, sr[O.TARGET_FROM_FIELD], sr[O.TARGET_TO_FIELD]
                            ),
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
            elif results[O.DIAMOND_FIELD] != O.EMPTY_RESULT_PLACEHOLDER:
                results[O.RECOMMEND_RESULTS_FIELD] = O.EC_SEPARATOR.join(
                    results[O.DIAMOND_FIELD].split(O.RESULT_SEPARATOR)
                )

        return results

    def graph_based_grouping(self, all_results: List[Tuple]) -> List[Tuple]:
        """Pick one representative annotation per conserved region.

        This groups *alternate annotations of the same region* (a dbCAN-HMM
        family and the dbCAN-sub subfamily that call it) so the recommendation
        does not list the same domain twice. It is not a second domain-conflict
        algorithm: dbCAN-HMM hits reaching this point already survived
        match-state-aware resolution, so two of them are distinct domains (or an
        explicitly allowed coexisting pair) and are never merged with each other.
        """
        if not all_results:
            return []

        hmm_results = sorted([r for r in all_results if r[3] == 'hmm'], key=lambda x: x[1])
        sub_results = sorted([r for r in all_results if r[3] == 'sub'], key=lambda x: x[1])

        # Special case 1: a single subfamily call spanning every HMM domain.
        if len(sub_results) == 1 and len(hmm_results) > 1:
            sub = sub_results[0]
            if all(self.annotations_overlap(sub, hmm) for hmm in hmm_results):
                best = self.select_best_result(sub_results + hmm_results)
                return [best]

        # Special case 2: one subfamily-style HMM call vs repeats of one subfamily.
        if len(hmm_results) == 1 and len(sub_results) > 1:
            hmm = hmm_results[0]
            if "_" in hmm[0]:
                sub_names = {s[0] for s in sub_results}
                if len(sub_names) == 1:
                    if all(self.annotations_overlap(hmm, sub) for sub in sub_results):
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
                    if any(self.annotations_overlap(g, sub2) for g in group):
                        group.append(sub2)
                        processed_subs.add(j)
                        changed = True
            sub_groups.append(group)

        processed_hmms = set()
        all_groups: List[List[Tuple]] = []

        for sub_group in sub_groups:
            group = sub_group.copy()
            for i, hmm in enumerate(hmm_results):
                if any(self.annotations_overlap(sub, hmm) for sub in sub_group):
                    group.append(hmm)
                    processed_hmms.add(i)
            all_groups.append(group)

        # Every remaining HMM hit is its own region. Grouping HMM hits with each
        # other here would re-collapse an architecture that domain resolution
        # already accepted — the exact bug that made a CBM sitting inside a
        # catalytic domain's insertion disappear from Recommend Results even
        # though dbCAN_hmm_results.tsv kept it.
        for i, hmm in enumerate(hmm_results):
            if i in processed_hmms:
                continue
            processed_hmms.add(i)
            all_groups.append([hmm])

        final_results = []
        for group in all_groups:
            if group:
                final_results.append(self.select_best_result(group))
        return final_results

    def aggregate_data(self, gene_ids: Iterable[str], indexed: Dict[str, Dict[str, pd.DataFrame]]) -> pd.DataFrame:
        cols = self.overview_columns
        aggregated: List[List[Any]] = []
        for gene_id in sorted(gene_ids):
            result = self.determine_best_result(gene_id, indexed)
            row: List[Any] = [gene_id]
            for col in cols[1:]:
                row.append(result.get(col, O.EMPTY_RESULT_PLACEHOLDER))
            aggregated.append(row)
        return pd.DataFrame(aggregated, columns=cols)

    # ---------------------------
    # Orchestration
    # ---------------------------
    def run(self):
        try:
            self.load_alignment_details()
            loaded_data = self.load_data()
            if not loaded_data:
                logger.warning("No annotation result files found. Creating empty overview.")
                empty_df = pd.DataFrame(columns=self.overview_columns)
                empty_df.to_csv(self.overview_output_path, sep='\t', index=False)
                logger.info(f"Empty overview saved to: {self.overview_output_path}")
                return

            gene_ids = set()
            indexed: Dict[str, Dict[str, pd.DataFrame]] = {}
            for key, df in loaded_data.items():
                if key in ("dbcan_hmm", "dbcan_sub"):
                    id_col = O.TARGET_NAME_FIELD
                else:
                    id_col = O.GENE_ID_FIELD
                if id_col in df.columns:
                    unique_ids = df[id_col].unique()
                    gene_ids.update(unique_ids)
                    # Build per-gene index to avoid repeated full-table scans
                    groups: Dict[str, pd.DataFrame] = {}
                    for gid, sub_df in df.groupby(id_col, sort=False):
                        groups[str(gid)] = sub_df
                    indexed[key] = groups

            aggregated_df = self.aggregate_data(gene_ids, indexed)
            aggregated_df.to_csv(self.overview_output_path, sep='\t', index=False)
            logger.info(f"Aggregated overview saved to: {self.overview_output_path}")

        except Exception as e:
            logger.error(f"Error generating overview: {e}", exc_info=True)
