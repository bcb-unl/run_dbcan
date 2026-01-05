import logging
from pathlib import Path
from typing import Dict, Tuple, List

import pandas as pd
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
        return ok

    def load_substrate_mapping(self) -> Dict[Tuple[str, str], List[str]]:
        """
        Load mapping: (family, EC) -> list of substrates (unique, case-insensitively sorted).
        Keep all candidates so that for INCOMPLETE EC (a.b.c.-) we can deterministically choose one.
        """
        try:
            df = pd.read_csv(self.mapping_file_path, sep='\t', header=None, skiprows=1, usecols=[2, 4, 0])
            df[2] = df[2].astype(str).str.strip()                      # family
            df[4] = df[4].astype(str).str.strip().replace({'NA': '-', 'nan': '-'}).fillna('-')  # EC
            df[0] = df[0].astype(str).str.strip()                      # substrate
            # Group to unique, sorted list per (family, EC)
            grp = df.groupby([2, 4])[0].apply(
                lambda s: sorted({v for v in s if v and v != '-'}, key=str.lower)
            )
            return {k: v for k, v in grp.items()}
        except FileNotFoundError:
            logger.warning(f"Can't find substrate mapping file: {self.mapping_file_path}")
            return {}
        except Exception as e:
            logger.error(f"Error loading substrate mapping: {e}", exc_info=True)
            return {}

    def process_dbcan_sub(self) -> None:
        if not self._validate_for_run():
            return

        subs_dict = self.load_substrate_mapping()
        if not subs_dict:
            logger.warning("No substrate mapping data loaded. Substrate annotation will be '-'.")

        try:
            df = pd.read_csv(self.input_file_path, sep='\t')
            if df.empty:
                logger.warning("No dbCAN-sub results to process")
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
                df[substrate_col] = df[hmm_name_col].apply(lambda x: self.get_substrates(str(x), subs_dict))

                # Sanity check: lengths must match per row; pad substrates if needed (no reordering)
                def _pad_match_len(ec_s: str, sub_s: str) -> str:
                    ecs = [t for t in (ec_s.split(';') if isinstance(ec_s, str) else []) if t != '']
                    subs = [t for t in (sub_s.split(';') if isinstance(sub_s, str) else [])]
                    if len(subs) < len(ecs):
                        subs = subs + ['-'] * (len(ecs) - len(subs))
                    elif len(subs) > len(ecs):
                        subs = subs[:len(ecs)]
                    return ';'.join(subs)
                df[substrate_col] = [
                    _pad_match_len(e, s) for e, s in zip(df[subfamily_ec_col].tolist(),
                                                         df[substrate_col].tolist())
                ]

                # Drop original HMM Name
                df.drop(columns=[hmm_name_col], inplace=True, errors='ignore')

            # Normalize derived columns
            for col in (subfamily_name_col, subfamily_comp_col, subfamily_ec_col, substrate_col):
                if col not in df.columns:
                    df[col] = '-'
                df[col] = df[col].fillna('-').astype(str)

            # Reorder/ensure final columns if constant provided
            out_cols = getattr(P, "DBCAN_SUB_COLUMN_NAMES", None)
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
        # Use the same sorted EC token list as substrate mapping to keep consistent order
        ecs_sorted = DBCANSUBProcessor._parse_sorted_ec_tokens(hmm_name)
        return ';'.join(ecs_sorted) if ecs_sorted else '-'

    @staticmethod
    def _parse_sorted_ec_tokens(hmm_name: str) -> List[str]:
        """
        Extract all EC tokens (format 'a.b.c.d:count') from HMM name and sort them deterministically:
        - Sort by EC core numerically (a, b, c, d), with '-' in d treated as infinity (placed after numbers).
        - Stable within equal keys (preserves original order when keys tie).
        Returns the sorted list of original tokens (including ':count').
        """
        parts = str(hmm_name).split('|')
        indexed = []
        for idx, p in enumerate(parts):
            if ':' not in p:
                continue
            ec_core = p.split(':')[0]
            segs = ec_core.split('.')
            if len(segs) != 4:
                continue
            indexed.append((idx, ec_core, p))

        def _conv(s):
            if s == '-':
                return 10**9  # push incomplete EC (d='-') to the end
            try:
                return int(s)
            except Exception:
                return s

        def _key(item):
            _, core, _ = item
            a, b, c, d = core.split('.')
            return (_conv(a), _conv(b), _conv(c), _conv(d))

        indexed_sorted = sorted(indexed, key=lambda it: (_key(it), it[0]))
        return [t for _, __, t in indexed_sorted]

    def get_substrates(self, profile_info: str, subs_dict: Dict[Tuple[str, str], List[str]]) -> str:
        """
        Map each EC token to a substrate individually.

        Rules (change applies to INCOMPLETE EC only):
        - Token format: a.b.c.d:count (we use a.b.c.d as EC core).
        - COMPLETE EC (d != '-'):
            1) pick from subs_dict[(family, EC)] — if multiple, choose the last one alphabetically
            2) fallback to subs_dict[(family, '-')] — same picking rule
            3) if the family has exactly one unique substrate overall -> use it; otherwise '-'
        - INCOMPLETE EC (d == '-'):
            Only look at subs_dict[(family, a.b.c.-)]; if multiple, choose the last one alphabetically.
            If not found, return '-' (NO fallback to subs_dict[(family, '-')]).
        - No EC tokens:
            1) subs_dict[(family, '-')]
            2) if the family has exactly one unique substrate overall -> use it; otherwise '-'
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

        # Family-level unique substrate set (used only for COMPLETE EC step-3 fallback)
        family_all_subs = {s for (fam, _ec), vals in subs_dict.items() if fam == family
                           for s in (vals if isinstance(vals, (list, set, tuple)) else [vals])}
        family_all_subs_clean = {s for s in family_all_subs if s and s != '-'}

        # Extract EC tokens in a deterministic order identical to _extract_subfamily_ecs
        ec_tokens = self._parse_sorted_ec_tokens(profile_info)

        # No EC tokens: family-level fallback
        if not ec_tokens:
            v = _pick(subs_dict.get((family, '-')))
            if v:
                return v
            if len(family_all_subs_clean) == 1:
                return next(iter(family_all_subs_clean))
            return '-'

        substrates_mapped = []
        for token in ec_tokens:
            ec_core = token.split(':')[0]
            segs = ec_core.split('.')
            is_incomplete = (segs[-1] == '-')

            if is_incomplete:
                # Only consider candidates under the exact key (family, a.b.c.-); pick last alphabetically if multiple
                v = _pick(subs_dict.get((family, ec_core)))
                if v:
                    substrates_mapped.append(v)
                    continue
                # No fallback to (family, '-'); strictly return '-'
                substrates_mapped.append('-')
                continue

            # COMPLETE EC handling (same logic, supports multi-candidate by picking last alphabetically)
            v = _pick(subs_dict.get((family, ec_core)))
            if v:
                substrates_mapped.append(v)
                continue
            v = _pick(subs_dict.get((family, '-')))
            if v:
                substrates_mapped.append(v)
                continue
            if len(family_all_subs_clean) == 1:
                substrates_mapped.append(next(iter(family_all_subs_clean)))
            else:
                substrates_mapped.append('-')

        return ';'.join(substrates_mapped) if substrates_mapped else '-'
