import logging
from pathlib import Path
from typing import Dict, Tuple, List

import pandas as pd
from dbcan.IO.table import read_tsv
from dbcan.configs.pyhmmer_config import DBCANSUBConfig
import dbcan.constants.process_dbcan_sub_constants as P

logger = logging.getLogger(__name__)

class DBCANSUBProcessor:
    """Process dbCAN-sub results (add substrate info). Config is the single source of truth."""

    def __init__(self, config: DBCANSUBConfig):
        self.config = config

    @property
    def input_file_path(self) -> Path:
        # Raw HMM output produced by hmmsearch
        return Path(self.config.output_dir) / P.DBCAN_SUB_HMM_RAW_FILE

    @property
    def output_file_path(self) -> Path:
        # Processed results path
        return Path(self.config.output_dir) / P.DBCAN_SUB_HMM_RESULT_FILE

    @property
    def mapping_file_path(self) -> Path:
        # Substrate mapping table path
        return Path(self.config.db_dir) / P.SUBSTRATE_MAPPING_FILE

    def _validate_for_run(self) -> bool:
        ok = True
        if not self.input_file_path.exists():
            logger.warning(f"dbCAN-sub raw file not found: {self.input_file_path}")
            ok = False
        elif self.input_file_path.stat().st_size == 0:
            logger.warning(f"dbCAN-sub raw file is empty: {self.input_file_path}")
            ok = False
        if not self.mapping_file_path.exists():
            logger.warning(f"Substrate mapping file not found: {self.mapping_file_path}. Substrate will be '-'.")
            # Mapping file missing is not fatal, we can still process with empty substrate column
        return ok

    def load_substrate_mapping(self) -> Dict[Tuple[str, str], List[str]]:
        """
        Load mapping: (family, EC) -> list of substrates (unique, case-insensitively sorted).
        Keep all candidates so that for INCOMPLETE EC (a.b.c.-) we can deterministically choose one.

        Supports both layouts of fam-substrate-mapping.tsv:
        - 5 columns: Substrate_high_level, Substrate_curated, Family, Name, EC_Number;
        - the curated web release (e.g. fam-substrate-mapping-08262025.tsv): leading '#' comment
          line(s) and two extra columns, new_Substrate_high_level (the current substrate) and
          type (NEW / CHANGED / REMOVED). REMOVED rows are dropped and new_Substrate_high_level
          is used when it is filled.
        Rows without a substrate are ignored. A file that exists but cannot be parsed raises.
        """
        if not self.mapping_file_path.exists():
            logger.warning(f"Can't find substrate mapping file: {self.mapping_file_path}")
            return {}
        df = read_tsv(self.mapping_file_path, id_columns=(), comment=None, dtype=str, quoting=0)
        # drop leading comment lines that pandas took as the header
        while len(df.columns) and str(df.columns[0]).startswith('#'):
            df.columns = [str(c).strip() for c in df.iloc[0]]
            df = df.iloc[1:].reset_index(drop=True)
        df.columns = [str(c).strip() for c in df.columns]
        required = {P.MAPPING_SUBSTRATE_COLUMN, P.MAPPING_FAMILY_COLUMN, P.MAPPING_EC_COLUMN}
        if not required.issubset(df.columns):
            raise ValueError(
                f"{self.mapping_file_path} lacks columns {sorted(required - set(df.columns))}; "
                f"found {list(df.columns)}"
            )
        df = df.fillna('')
        if P.MAPPING_TYPE_COLUMN in df.columns:
            df = df[df[P.MAPPING_TYPE_COLUMN].str.strip().str.upper() != 'REMOVED']
        substrate = df[P.MAPPING_SUBSTRATE_COLUMN].str.strip()
        if P.MAPPING_NEW_SUBSTRATE_COLUMN in df.columns:
            new_sub = df[P.MAPPING_NEW_SUBSTRATE_COLUMN].str.strip()
            substrate = new_sub.where(new_sub != '', substrate)
        family = df[P.MAPPING_FAMILY_COLUMN].str.strip()
        ec = df[P.MAPPING_EC_COLUMN].str.strip().replace({'': '-', 'NA': '-', 'nan': '-'})
        out: Dict[Tuple[str, str], List[str]] = {}
        for fam, e, sub in zip(family, ec, substrate):
            if not fam or not sub or sub == '-':
                continue
            out.setdefault((fam, e), [])
            if sub not in out[(fam, e)]:
                out[(fam, e)].append(sub)
        for k in out:
            out[k] = sorted(out[k], key=str.lower)
        return out

    def process_dbcan_sub(self) -> None:
        logger.info(f"Starting dbCAN-sub processing. Input: {self.input_file_path}, Output: {self.output_file_path}")
        # Get output columns for creating empty file if needed
        out_cols = getattr(P, "DBCAN_SUB_COLUMN_NAMES", None)
        if out_cols is None:
            logger.error(f"DBCAN_SUB_COLUMN_NAMES not found in constants. Cannot create empty file with proper headers.")
            # Still try to create a basic empty file
            try:
                basic_cols = ['Subfam Name', 'Subfam Composition', 'Subfam EC', 'Substrate', 
                             'HMM Length', 'Target Name', 'Target Length', 'i-Evalue',
                             'HMM From', 'HMM To', 'Target From', 'Target To', 'Coverage', 'HMM File Name']
                empty_df = pd.DataFrame(columns=basic_cols)
                self.output_file_path.parent.mkdir(parents=True, exist_ok=True)
                empty_df.to_csv(self.output_file_path, sep='\t', index=False)
                logger.info(f"Created empty dbCAN-sub results file with basic headers -> {self.output_file_path}")
            except Exception as e:
                logger.error(f"Failed to create basic empty file: {e}", exc_info=True)
            return
        else:
            logger.info(f"DBCAN_SUB_COLUMN_NAMES found: {len(out_cols)} columns: {out_cols[:5]}...")
        
        if not self._validate_for_run():
            logger.info(f"Validation failed. Input file exists: {self.input_file_path.exists()}, "
                       f"size: {self.input_file_path.stat().st_size if self.input_file_path.exists() else 0}")
            # Create empty file with headers when input file doesn't exist or is empty
            if out_cols:
                try:
                    empty_df = pd.DataFrame(columns=out_cols)
                    self.output_file_path.parent.mkdir(parents=True, exist_ok=True)
                    empty_df.to_csv(self.output_file_path, sep='\t', index=False)
                    logger.info(f"Created empty dbCAN-sub results file with headers -> {self.output_file_path}")
                except Exception as e:
                    logger.error(f"Failed to create empty file: {e}", exc_info=True)
            else:
                logger.warning(f"Cannot create empty file: DBCAN_SUB_COLUMN_NAMES is None")
            return

        subs_dict = self.load_substrate_mapping()
        subfam_table = self.load_subfamily_substrate_table()
        if not subs_dict:
            logger.warning("No substrate mapping data loaded. Substrate annotation will be '-'.")

        try:            
            df = read_tsv(self.input_file_path)
            if df.empty:
                logger.warning("No dbCAN-sub results to process (DataFrame is empty)")
                # Create empty file with headers when DataFrame is empty
                if out_cols:
                    try:
                        empty_df = pd.DataFrame(columns=out_cols)
                        self.output_file_path.parent.mkdir(parents=True, exist_ok=True)
                        empty_df.to_csv(self.output_file_path, sep='\t', index=False)
                        logger.info(f"Created empty dbCAN-sub results file with headers (no data rows) -> {self.output_file_path}")
                    except Exception as e:
                        logger.error(f"Failed to create empty file: {e}", exc_info=True)
                else:
                    logger.warning(f"Cannot create empty file: DBCAN_SUB_COLUMN_NAMES is None")
                return

            hmm_name_col = P.DBCAN_SUB_HMM_NAME_COLUMN
            subfamily_name_col = P.DBCAN_SUB_SUBFAMILY_NAME_COLUMN
            subfamily_comp_col = P.DBCAN_SUB_SUBFAMILY_COMP_COLUMN
            subfamily_ec_col = P.DBCAN_SUB_SUBFAMILY_EC_COLUMN
            substrate_col = P.DBCAN_SUB_SUBSTRATE_COLUMN

            if hmm_name_col not in df.columns:
                logger.warning(f"Column '{hmm_name_col}' not found in raw. Filling derived columns with '-' and writing back.")
                for col in (subfamily_name_col, subfamily_comp_col, subfamily_ec_col, substrate_col):
                    if col not in df.columns:
                        df[col] = '-'
            else:
                # Derive columns from HMM Name
                df[subfamily_name_col] = df[hmm_name_col].apply(self._extract_subfamily_names)
                df[subfamily_comp_col] = df[hmm_name_col].apply(self._extract_subfamily_components)
                df[subfamily_ec_col] = df[hmm_name_col].apply(self._extract_subfamily_ecs)
                df[substrate_col] = df[hmm_name_col].apply(lambda x: self.get_substrates(str(x), subs_dict, subfam_table))

                # Drop original HMM Name
                df.drop(columns=[hmm_name_col], inplace=True, errors='ignore')

            # Normalize derived columns
            for col in (subfamily_name_col, subfamily_comp_col, subfamily_ec_col, substrate_col):
                if col not in df.columns:
                    df[col] = '-'
                df[col] = df[col].fillna('-').astype(str)

            # Reorder/ensure final columns if constant provided
            if out_cols:
                for c in out_cols:
                    if c not in df.columns:
                        df[c] = '-'
                df = df[out_cols]

            self.output_file_path.parent.mkdir(parents=True, exist_ok=True)
            df.to_csv(self.output_file_path, sep='\t', index=False)
            logger.info(f"Successfully processed dbCAN-sub results ({len(df)} entries) -> {self.output_file_path.name}")
        except Exception as e:
            logger.error(f"Error processing dbCAN-sub results: {e}", exc_info=True)
            # Create empty file with headers even when there's an error
            if out_cols:
                try:
                    empty_df = pd.DataFrame(columns=out_cols)
                    self.output_file_path.parent.mkdir(parents=True, exist_ok=True)
                    empty_df.to_csv(self.output_file_path, sep='\t', index=False)
                    logger.info(f"Created empty dbCAN-sub results file with headers after error -> {self.output_file_path}")
                except Exception as write_error:
                    logger.error(f"Failed to create empty file: {write_error}", exc_info=True)
        
        # Final check: ensure output file exists
        if not self.output_file_path.exists():
            logger.warning(f"Output file does not exist after processing: {self.output_file_path}. Creating empty file as fallback.")
            if out_cols:
                try:
                    empty_df = pd.DataFrame(columns=out_cols)
                    self.output_file_path.parent.mkdir(parents=True, exist_ok=True)
                    empty_df.to_csv(self.output_file_path, sep='\t', index=False)
                    logger.info(f"Created empty dbCAN-sub results file as final fallback -> {self.output_file_path}")
                except Exception as write_error:
                    logger.error(f"Failed to create empty file in final fallback: {write_error}", exc_info=True)
            else:
                logger.error(f"Cannot create empty file: DBCAN_SUB_COLUMN_NAMES is None")
        else:
            logger.info(f"Output file exists: {self.output_file_path} (size: {self.output_file_path.stat().st_size} bytes)")

    @staticmethod
    def _extract_subfamily_names(hmm_name: str) -> str:
        parts = str(hmm_name).split('|')
        names = [p.split('.hmm')[0] for p in parts if p.endswith('.hmm')]
        return ';'.join(names) if names else '-'

    @staticmethod
    def _extract_subfamily_components(hmm_name: str) -> str:
        parts = str(hmm_name).split('|')
        # Keep segments that are not .hmm names and do not look like EC (not 4-part dotted)
        comps = [p for p in parts if not p.endswith('.hmm') and len(p.split('.')) != 4]
        return ';'.join(comps) if comps else '-'

    @staticmethod
    def _extract_subfamily_ecs(hmm_name: str) -> str:
        """EC tokens as 'EC:count', in the same order as the Substrate column.

        eHMM names omit the count when only one characterized protein has the EC
        (e.g. 'AA1_e1.hmm|AA1:57|AA1_3:52|1.10.3.2'); the dbCAN-sub web server lists such
        an EC with count 1, so it is written as '1.10.3.2:1'.
        """
        ecs_sorted = DBCANSUBProcessor._parse_sorted_ec_tokens(hmm_name)
        ecs = [t if ':' in t else f"{t}:1" for t in ecs_sorted]
        return ';'.join(ecs) if ecs else '-'

    @staticmethod
    def _parse_sorted_ec_tokens(hmm_name: str) -> List[str]:
        """
        Extract all EC tokens from an eHMM name and sort them deterministically.

        A token is ``a.b.c.d`` optionally followed by ``:count`` (e.g. ``3.2.1.4:18`` or
        ``1.10.3.2``); each of b, c, d may be ``-`` and d may be a preliminary ``nX`` number.
        - Sort by EC core numerically (a, b, c, d), with '-' treated as infinity (placed after numbers).
        - Stable within equal keys (preserves original order when keys tie).
        Returns the sorted list of original tokens (including ':count' when present).
        """
        parts = str(hmm_name).split('|')
        indexed = []
        for idx, p in enumerate(parts):
            p = p.strip()
            if p.endswith('.hmm'):
                continue
            ec_core = p.split(':', 1)[0]
            if not P.EC_CORE_PATTERN.match(ec_core):
                continue
            indexed.append((idx, ec_core, p))

        def _conv(s):
            if s == '-':
                return 10**9  # push incomplete EC (d='-') to the end
            try:
                return int(s)
            except Exception:
                return 10**8  # e.g. preliminary 'n1' numbers: after numbers, before '-'

        def _key(item):
            _, core, _ = item
            a, b, c, d = core.split('.')
            return (_conv(a), _conv(b), _conv(c), _conv(d))

        indexed_sorted = sorted(indexed, key=lambda it: (_key(it), it[0]))
        return [t for _, __, t in indexed_sorted]

    def get_substrates(self, profile_info: str, subs_dict: Dict[Tuple[str, str], List[str]],
                       subfam_table: Dict[Tuple[str, str], str] | None = None) -> str:
        """
        Map each EC token of an eHMM to a substrate.

        The family is taken from the eHMM name (``GH71_e32.hmm`` -> ``GH71``). EC tokens are
        ``a.b.c.d`` or ``a.b.c.-``, with or without ``:count``.

        - COMPLETE EC (d != '-'):
            1) subs_dict[(family, EC)]; if several substrates, the last one alphabetically
            2) fallback: subs_dict[(family, '-')] (family-level row of the mapping file)
            3) fallback: the family's only substrate, if the family has exactly one
            otherwise '-'
        - INCOMPLETE EC (d == '-'): only subs_dict[(family, a.b.c.-)], otherwise '-'.
        - No EC token: subs_dict[(family, '-')], else the family's only substrate, else '-'.

        An EC number is looked up under every family of the eHMM composition (first the
        eHMM's own family, then the other families by member count): the EC of a CBM or
        multi-domain eHMM usually belongs to its catalytic family (e.g. CBM1_e3 with
        GH45:49 and 3.2.1.4). If ``subfam_table`` (the dbCAN-sub web server's per-subfamily
        substrate table, see DBCANSUB_SUBSTRATE_TABLE_FILE) has an entry for
        (subfamily, EC), that entry is used as is, so results match the web server.

        Returns one entry per EC token, joined with ';' in the same order as the
        ``Subfam EC`` column. A mapping value that itself lists several substrates
        (``beta-glucan;beta-mannan``) is written with ',' inside its entry so that the
        entries stay aligned with the EC tokens.
        """
        def _pick(v) -> str | None:
            """Pick a deterministic substrate from a list-like or single value."""
            if v is None:
                return None
            if isinstance(v, (list, set, tuple)):
                candidates = [x for x in v if x and x != '-']
                if not candidates:
                    return None
                return sorted(candidates, key=str.lower)[-1]
            # Legacy single value
            return v if isinstance(v, str) and v and v != '-' else None

        if not profile_info or not isinstance(profile_info, str):
            return '-'
        parts = profile_info.split('|')
        if not parts:
            return '-'

        try:
            # Family name inferred from the first HMM token (strip '.hmm' and optional suffix parts)
            family = parts[0].split('.hmm')[0].split("_")[0]
        except (IndexError, AttributeError):
            return '-'

        # Family-level unique substrate set (used by the fallbacks below)
        family_all_subs = {s for (fam, _ec), vals in subs_dict.items() if fam == family
                           for s in (vals if isinstance(vals, (list, set, tuple)) else [vals])}
        family_all_subs_clean = {s for s in family_all_subs if s and s != '-'}

        def _family_fallback() -> str | None:
            v = _pick(subs_dict.get((family, '-')))
            if v:
                return v
            if len(family_all_subs_clean) == 1:
                return next(iter(family_all_subs_clean))
            return None

        def _entry(v) -> str:
            return v.replace(';', ',') if v else '-'

        ec_tokens = self._parse_sorted_ec_tokens(profile_info)
        if not ec_tokens:
            return _entry(_family_fallback())

        subfamily = parts[0].split('.hmm')[0]
        lookup_families = [family] + [f for f in self._composition_families(profile_info) if f != family]

        substrates_mapped = []
        for token in ec_tokens:
            ec_core = token.split(':', 1)[0]
            if subfam_table is not None and (subfamily, ec_core) in subfam_table:
                substrates_mapped.append(_entry(subfam_table[(subfamily, ec_core)]))
                continue
            v = None
            for fam in lookup_families:
                v = _pick(subs_dict.get((fam, ec_core)))
                if v:
                    break
            if not v and ec_core.split('.')[-1] != '-':
                v = _family_fallback()
            substrates_mapped.append(_entry(v))

        return ';'.join(substrates_mapped)

    @staticmethod
    def _composition_families(hmm_name: str) -> List[str]:
        """CAZy families of an eHMM composition ('GH7:184', 'GH5_5:4' -> GH5), most members first."""
        counts: Dict[str, int] = {}
        for p in str(hmm_name).split('|')[1:]:
            key, _, cnt = p.partition(':')
            m = P.FAMILY_TOKEN_PATTERN.match(key.strip())
            if not m:
                continue
            fam = m.group(1)
            counts[fam] = counts.get(fam, 0) + (int(cnt) if cnt.isdigit() else 0)
        return [f for f, _ in sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))]

    def load_subfamily_substrate_table(self) -> Dict[Tuple[str, str], str] | None:
        """Optional per-subfamily table exported from the dbCAN-sub web server.

        Columns: Subfamily, EC, Count, Substrate (one row per subfamily and EC; Substrate may be
        empty or '-'). Returns None when the file is not in db_dir.
        """
        path = Path(self.config.db_dir) / P.DBCANSUB_SUBSTRATE_TABLE_FILE
        if not path.exists():
            return None
        df = read_tsv(path, id_columns=(), dtype=str).fillna('')
        missing = {"Subfamily", "EC", "Substrate"} - set(df.columns)
        if missing:
            raise ValueError(f"{path} lacks columns {sorted(missing)}")
        table = {}
        for sub, ec, substrate in zip(df["Subfamily"].str.strip(), df["EC"].str.strip(), df["Substrate"].str.strip()):
            table[(sub, ec)] = substrate if substrate and substrate != '-' else ''
        logger.info(f"Loaded per-subfamily substrate table: {len(table)} (subfamily, EC) entries from {path.name}")
        return table
