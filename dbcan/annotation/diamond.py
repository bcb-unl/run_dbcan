from pathlib import Path
import subprocess
import logging
import pandas as pd
from abc import ABC
import time

from dbcan.configs.diamond_config import (
    DiamondConfig, DiamondCAZyConfig,DiamondTCConfig, DiamondSulfataseConfig, DiamondPeptidaseConfig, DiamondTFConfig
)
import dbcan.constants.diamond_constants as D

logger = logging.getLogger(__name__)

class DiamondProcessor(ABC):
    def __init__(self, config: DiamondConfig):
        self.config = config
        self._validate()

    @property
    def diamond_db(self) -> str:
        return str(Path(self.config.db_dir) / self.config.db_file)

    @property
    def input_faa(self) -> str:
        return str(Path(self.config.output_dir) / self.config.input_faa)

    @property
    def output_file(self) -> str:
        return str(Path(self.config.output_dir) / self.config.output_file)

    @property
    def e_value_threshold(self) -> str:
        return str(self.config.e_value_threshold)

    @property
    def coverage_threshold(self) -> str:
        return str(self.config.coverage_threshold)

    @property
    def verbose_option(self) -> bool:
        return self.config.verbose_option

    def _validate(self):
        if self.config.e_value_threshold is None:
            raise ValueError("e_value_threshold must be set in the configuration")
        if not Path(self.diamond_db).exists():
            raise FileNotFoundError(f"Database file not found: {self.diamond_db}")
        if not Path(self.input_faa).exists():
            raise FileNotFoundError(f"Input file not found: {self.input_faa}")
        Path(self.output_file).parent.mkdir(parents=True, exist_ok=True)

    def run(self):
        outfmt_args = [D.DIAMOND_CMD_OUTFMT, D.DIAMOND_DEFAULT_OUTFMT]
        if self.config.outfmt_fields:
            outfmt_args.extend(self.config.outfmt_fields)

        cmd = [
            D.DIAMOND_CMD, D.DIAMOND_BLASTP_CMD,
            D.DIAMOND_CMD_DB, self.diamond_db,
            D.DIAMOND_CMD_QUERY, self.input_faa,
            D.DIAMOND_CMD_OUT, self.output_file,
            *outfmt_args,
            D.DIAMOND_CMD_EVALUE, str(self.config.e_value_threshold),
            D.DIAMOND_CMD_MAX_TARGET, str(D.DIAMOND_MAX_TARGET_SEQS),
            D.DIAMOND_CMD_THREADS, str(self.config.threads),
            D.DIAMOND_CMD_VERBOSE if self.verbose_option else D.DIAMOND_CMD_QUIET
        ]
        if self.config.coverage_threshold is not None:
            cmd.extend([D.DIAMOND_CMD_QUERY_COVER, str(self.config.coverage_threshold)])

        logger.info("Running DIAMOND BLASTp: %s", " ".join(cmd))
        start = time.time()
        try:
            result = subprocess.run(cmd, check=False, capture_output=True, text=True)
            elapsed = time.time() - start
            if result.returncode != 0:
                logger.error("DIAMOND BLASTp failed (exit %s)", result.returncode)
                if result.stderr:
                    logger.error("DIAMOND stderr: %s", result.stderr.strip())
                raise subprocess.CalledProcessError(result.returncode, cmd, output=result.stdout, stderr=result.stderr)

            logger.info("DIAMOND BLASTp completed in %.2fs", elapsed)
            self.format_results()
        except subprocess.CalledProcessError as e:
            logger.error("DIAMOND BLASTp failed: %s", e)
            raise

    def format_results(self):
        of = Path(self.output_file)
        if not of.exists() or of.stat().st_size == 0:
            # Create a header-only TSV to avoid downstream warnings
            cols = list(self.config.column_names) if getattr(self.config, "column_names", None) else []
            # Keep the 'Database' column convention if label is provided
            if self.config.label and 'Database' not in cols:
                cols.append('Database')
            header_df = pd.DataFrame(columns=cols)
            of.parent.mkdir(parents=True, exist_ok=True)
            header_df.to_csv(of, sep='\t', index=False)
            logger.info(f"No DIAMOND hits. Wrote header-only file: {self.output_file}")
            return
        try:
            df = pd.read_csv(of, sep='\t', header=None, names=self.config.column_names)
            if self.config.id_column:
                self._postprocess_ids(df, self.config.id_column)
            if self.config.label:
                df['Database'] = self.config.label
            df.to_csv(of, sep='\t', index=False)
            logger.info(f"Results formatted and saved to {self.output_file}")
        except Exception as e:
            logger.error(f"Error formatting results: {e}")
            raise

    def _postprocess_ids(self, df: pd.DataFrame, id_col: str):
        pass
class CAZYDiamondProcessor(DiamondProcessor):
    def __init__(self, config: DiamondCAZyConfig):
        super().__init__(config)

    def _postprocess_ids(self, df: pd.DataFrame, id_col: str):
        df[id_col] = df[id_col].apply(lambda x: x.split(' ')[0].split('|')[-1] if isinstance(x, str) else x)



class TCDBDiamondProcessor(DiamondProcessor):
    def __init__(self, config: DiamondTCConfig):
        super().__init__(config)

    def _postprocess_ids(self, df: pd.DataFrame, id_col: str):
        df[id_col] = df[id_col].apply(lambda x: x.split(' ')[0].split('|')[-1] if isinstance(x, str) else x)

class SulfatlasDiamondProcessor(DiamondProcessor):
    def __init__(self, config: DiamondSulfataseConfig):
        super().__init__(config)

    def _postprocess_ids(self, df: pd.DataFrame, id_col: str):
        df[id_col] = df[id_col].apply(
            lambda x: "_".join(x.split('|')[1].split('_')[1:]) if isinstance(x, str) and '|' in x else "unknown"
        )

class PeptidaseDiamondProcessor(DiamondProcessor):
    def __init__(self, config: DiamondPeptidaseConfig):
        super().__init__(config)

    def _postprocess_ids(self, df: pd.DataFrame, id_col: str):
        df[id_col] = df[id_col].apply(lambda x: x.split('|')[1] if isinstance(x, str) and "|" in x else "unknown")

class TFDiamondProcessor(DiamondProcessor):
    def __init__(self, config: DiamondTFConfig):
        super().__init__(config)
    def _postprocess_ids(self, df: pd.DataFrame, id_col: str):
        df[id_col] = df[id_col].apply(lambda x: x.split('|')[2] if isinstance(x, str) and "|" in x else "unknown")
