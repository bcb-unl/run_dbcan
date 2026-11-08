import logging
import re
from pathlib import Path
from typing import Dict, List, Iterable, Tuple, Any

import pandas as pd
from dbcan.IO.table import read_tsv
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

_CAZY_FAMILY_RE = re.compile(O.CAZY_ID_PATTERN)

# CAZyme_annotation method name (--methods, lower-cased) -> result key used below
METHOD_TO_RESULT_KEY = {"diamond": "diamond", "hmm": "dbcan_hmm", "dbcansub": "dbcan_sub"}
RESULT_KEY_TO_METHOD = {"diamond": "diamond", "dbcan_hmm": "hmm", "dbcan_sub": "dbCANsub"}

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
    def selected_result_keys(self) -> set:
        """Result tables to merge: only the methods run in this CAZyme_annotation call.

        Files of other methods may be left in output_dir by an earlier run; merging them would
        report stale hits and inflate #ofTools. ``methods=None`` (direct API use) keeps all three.
        """
        methods = getattr(self.config, "methods", None)
        if methods is None:
            return set(self.result_file_map)
        if isinstance(methods, str):
            methods = [m.strip() for m in methods.split(",") if m.strip()]
        keys = set()
        for m in methods:
            key = METHOD_TO_RESULT_KEY.get(str(m).lower())
            if key is None:
                raise ValueError(f"Unknown CAZyme annotation method: {m!r} (expected diamond, hmm, dbCANsub)")
            keys.add(key)
        return keys

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
        selected = self.selected_result_keys
        for key, path in self.result_file_map.items():
            if key not in selected:
                if path.exists():
                    logger.warning(
                        f"Ignoring {path.name}: method '{RESULT_KEY_TO_METHOD[key]}' was not selected with --methods in this run "
                        f"(the file is left over from an earlier run)."
                    )
                continue
            if not path.exists():
                logger.warning(f"{key} results not found at {path}")
                continue
            try:
                # Use chunked reading for large files (>100MB)
                file_size_mb = path.stat().st_size / (1024 * 1024)
                if file_size_mb > 100:
                    logger.info(f"Large file detected ({file_size_mb:.1f}MB), using chunked reading for {key}")
                    chunks = []
                    for chunk in read_tsv(path, chunksize=100000):
                        chunks.append(chunk)
                    df = pd.concat(chunks, ignore_index=True)
                else:
                    df = read_tsv(path)
                
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
        CAZy families of one CAZy DIAMOND subject header (sseqid), joined with O.RESULT_SEPARATOR.

        The header is '|' separated: an accession followed by families and EC numbers
        ('CAD7272370.1|GH33', '1803517A|GH13_6|GH13_6|3.2.1.1') or, for JGI entries, families
        followed by a numeric protein ID and the source file
        ('GH15|CBM20|716895|Lacbi81306_1_GeneCatalog_proteins_20140617.aa.fasta').

        Only whole family tokens (O.CAZY_ID_PATTERN, e.g. GH13 or GH13_6) are returned, in order
        and without repeats, so accessions such as 'AAK49401.1' (starts with 'AA'), EC numbers,
        protein numbers and file names are never reported as families. If no token is a family
        the value is returned unchanged.
        """
        if not isinstance(cazy_id, str) or not cazy_id:
            return cazy_id
        families = [t for t in (p.strip() for p in cazy_id.split('|')) if _CAZY_FAMILY_RE.match(t)]
        if not families:
            return cazy_id
        return O.RESULT_SEPARATOR.join(dict.fromkeys(families))

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
                    # one entry per EC token; '-' marks an EC without a substrate
                    uniq = sorted(
                        {p.strip() for p in flat if p.strip() and p.strip() != O.EMPTY_RESULT_PLACEHOLDER}
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
        algorithm: dbCAN-HMM hits reaching this point already survived domain
        resolution, so two of them are distinct domains (or an explicitly allowed
        coexisting pair) and are never merged with each other.
        """
        if not all_results:
            return []

        hmm_results = sorted([r for r in all_results if r[3] == 'hmm'], key=lambda x: x[1])
        sub_results = sorted([r for r in all_results if r[3] == 'sub'], key=lambda x: x[1])

        # Special case 1: a single subfamily call spanning every HMM domain.
        if len(sub_results) == 1 and len(hmm_results) > 1:
            sub = sub_results[0]
            if all(self.calculate_overlap(sub[1], sub[2], hmm[1], hmm[2]) for hmm in hmm_results):
                best = self.select_best_result(sub_results + hmm_results)
                return [best]

        # Special case 2: one subfamily-style HMM call vs repeats of one subfamily.
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

        # Every remaining HMM hit is its own region: grouping HMM hits with each
        # other here would drop a domain that dbCAN_hmm_results.tsv kept (e.g. a
        # pair allowed by a coexistence policy).
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
                    # IDs are opaque strings (issue #82): the same str() key is used for
                    # the gene list and for the per-gene index.
                    df = df.assign(**{id_col: df[id_col].astype(str)})
                    gene_ids.update(df[id_col].unique())
                    # Build per-gene index to avoid repeated full-table scans
                    groups: Dict[str, pd.DataFrame] = {}
                    for gid, sub_df in df.groupby(id_col, sort=False):
                        groups[gid] = sub_df
                    indexed[key] = groups

            aggregated_df = self.aggregate_data(gene_ids, indexed)
            aggregated_df.to_csv(self.overview_output_path, sep='\t', index=False)
            logger.info(f"Aggregated overview saved to: {self.overview_output_path}")

        except Exception as e:
            logger.error(f"Error generating overview: {e}", exc_info=True)
            raise
