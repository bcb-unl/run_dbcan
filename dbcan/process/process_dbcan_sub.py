import logging
from pathlib import Path
from typing import Dict, Tuple

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
        # raw hmm output produced by hmmsearch
        return Path(self.config.output_dir) / P.DBCAN_SUB_HMM_RAW_FILE

    @property
    def output_file_path(self) -> Path:
        # processed results
        return Path(self.config.output_dir) / P.DBCAN_SUB_HMM_RESULT_FILE

    @property
    def mapping_file_path(self) -> Path:
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

    def load_substrate_mapping(self) -> Dict[Tuple[str, str], str]:
        try:
            df = pd.read_csv(self.mapping_file_path, sep='\t', header=None, skiprows=1, usecols=[2, 4, 0])
            df[2] = df[2].astype(str).str.strip()
            df[4] = df[4].astype(str).str.strip().replace({'NA': '-', 'nan': '-'}).fillna('-')
            df[0] = df[0].astype(str).str.strip()
            df['key'] = df.apply(lambda x: (x[2], x[4]), axis=1)
            ser = pd.Series(df[0].values, index=pd.MultiIndex.from_tuples(df['key']))
            return ser.groupby(level=[0, 1]).last().to_dict()
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
                # Drop original HMM Name
                df.drop(columns=[hmm_name_col], inplace=True, errors='ignore')

            # Normalize derived columns
            for col in (subfamily_name_col, subfamily_comp_col, subfamily_ec_col, substrate_col):
                if col not in df.columns:
                    df[col] = '-'
                df[col] = df[col].fillna('-').astype(str)

            # Reorder/ensure final columns
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
        comps = [p for p in parts if not p.endswith('.hmm') and len(p.split('.')) != 4]
        return ';'.join(comps) if comps else '-'

    @staticmethod
    def _extract_subfamily_ecs(hmm_name: str) -> str:
        parts = str(hmm_name).split('|')
        ecs = [p for p in parts if ':' in p and len(p.split(':')[0].split('.')) == 4]
        return ';'.join(ecs) if ecs else '-'

    def get_substrates(self, profile_info: str, subs_dict: Dict[Tuple[str, str], str]) -> str:
        """
        Map EACH EC token to a substrate individually.

        Updated rules refinement:
        - Token format: a.b.c.d:count (we use a.b.c.d as EC core).
        - For COMPLETE EC (d != '-'):
            1) (family, EC) -> substrate
            2) fallback (family, '-') if exists
            3) if family only has exactly one unique substrate (excluding '-') -> that substrate
            4) else '-'
        - For INCOMPLETE EC (d == '-'):
            Attempt mapping instead of unconditional 'unknown':
            1) (family, EC) direct lookup (e.g. (GH17, 2.4.1.-) )
            2) fallback (family, '-')
            3) if family only has one unique substrate -> that substrate
            4) else 'unknown'
          (Reason: families like GH17 may have only one substrate; incomplete EC still informative)
        - No EC tokens at all:
            1) (family, '-')
            2) if family only one unique substrate -> that substrate
            3) else '-'
        - Output preserves EC token order; one substrate per EC token separated by ';'.
        """
        if not profile_info or not isinstance(profile_info, str):
            return '-'
        parts = profile_info.split('|')
        if not parts:
            return '-'

        try:
            family = parts[0].split('.hmm')[0].split("_")[0]
        except (IndexError, AttributeError):
            return '-'

        # Collect all substrates mapped to this family (for uniqueness fallback)
        family_all_subs = {v for (fam, _ec), v in subs_dict.items() if fam == family}
        family_all_subs_clean = {s for s in family_all_subs if s and s != '-'}

        # Extract ordered EC tokens
        ec_tokens = []
        for p in parts:
            if ':' in p:
                ec_core = p.split(':')[0]
                if len(ec_core.split('.')) == 4:
                    ec_tokens.append(p)

        # No EC tokens: family-level fallback
        if not ec_tokens:
            if (family, '-') in subs_dict:
                return subs_dict[(family, '-')]
            if len(family_all_subs_clean) == 1:
                return next(iter(family_all_subs_clean))
            return '-'

        substrates_mapped = []
        for token in ec_tokens:
            ec_core = token.split(':')[0]
            segs = ec_core.split('.')
            is_incomplete = (segs[-1] == '-')

            if is_incomplete:
                # Incomplete EC fallback chain
                direct = subs_dict.get((family, ec_core))
                if direct:
                    substrates_mapped.append(direct)
                    continue
                fallback = subs_dict.get((family, '-'))
                if fallback:
                    substrates_mapped.append(fallback)
                    continue
                if len(family_all_subs_clean) == 1:
                    substrates_mapped.append(next(iter(family_all_subs_clean)))
                else:
                    substrates_mapped.append('unknown')
                continue

            # COMPLETE EC handling
            direct = subs_dict.get((family, ec_core))
            if direct:
                substrates_mapped.append(direct)
                continue
            fallback = subs_dict.get((family, '-'))
            if fallback:
                substrates_mapped.append(fallback)
                continue
            if len(family_all_subs_clean) == 1:
                substrates_mapped.append(next(iter(family_all_subs_clean)))
            else:
                substrates_mapped.append('-')

        return ';'.join(substrates_mapped) if substrates_mapped else '-'
