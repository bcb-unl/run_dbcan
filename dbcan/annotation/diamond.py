import os
import subprocess
import logging
import pandas as pd

from dbcan.parameter import DiamondConfig, DiamondTCConfig, DiamondSulfataseConfig, DiamondPeptidaseConfig
from dbcan.constants import (TC,TCDB_ID_COLUMN,SULFATLAS_ID_COLUMN,PEPTIDASE_ID_COLUMN,

    CAZY_COLUMN_NAMES, TCDB_COLUMN_NAMES, SULFATLAS_COLUMN_NAMES, PEPTIDASE_COLUMN_NAMES,
                             CAZY_DIAMOND_DB, TCDB_DIAMOND_DB,
                             CAZY_DIAMOND_OUTPUT, TCDB_DIAMOND_OUTPUT,DIAMOND_CMD,DIAMOND_BLASTP_CMD,
                             DIAMOND_CAZY_EVALUE_DEFAULT, DIAMOND_TCDB_EVALUE_DEFAULT,
                             DIAMOND_TCDB_COVERAGE_DEFAULT, DIAMOND_MAX_TARGET_SEQS,
                             DIAMOND_DEFAULT_OUTFMT,
                             DIAMOND_CMD_DB, DIAMOND_CMD_QUERY, DIAMOND_CMD_OUT,
                             DIAMOND_CMD_OUTFMT, DIAMOND_CMD_EVALUE,
                             DIAMOND_CMD_MAX_TARGET, DIAMOND_CMD_THREADS,
                             DIAMOND_CMD_VERBOSE, DIAMOND_CMD_QUIET,
                             DIAMOND_CMD_QUERY_COVER, TCDB_DIAMOND_OUTFMT_FIELDS, SULFATLAS_DIAMOND_OUTFMT_FIELDS, PEPTIDASE_DIAMOND_OUTFMT_FIELDS,
                                SULFATLAS, SULFATLAS_DIAMOND_DB,
                                SULFATLAS_DIAMOND_OUTPUT, SULFATLAS_COLUMN_NAMES,
                                DIAMOND_SULFATLAS_EVALUE_DEFAULT, DIAMOND_SULFATLAS_COVERAGE_DEFAULT,
                                PEPTIDASE, PEPTIDASE_DIAMOND_DB,
                                PEPTIDASE_DIAMOND_OUTPUT, PEPTIDASE_COLUMN_NAMES,
                                DIAMOND_PEPTIDASE_EVALUE_DEFAULT, DIAMOND_PEPTIDASE_COVERAGE_DEFAULT)

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class DiamondProcessor:
    """Base Diamond processor class using template method pattern"""

    def __init__(self, config):
        """Initialize with configuration"""
        self.config = config
        self._setup_processor()

    def _setup_processor(self):
        """Set up processor attributes using template method pattern"""
        self.diamond_db = self._derive_diamond_db()
        self.input_faa = self._derive_input_faa()
        self.output_file = self._derive_output_file()
        self.e_value_threshold = self._derive_e_value_threshold()
        self.threads = self._derive_threads()
        self.verbose_option = self._derive_verbose_option()

        # Validate required attributes
        self._validate_attributes()

    def _validate_attributes(self):
        """Validate that all required attributes are properly set"""
        required_attrs = ['diamond_db', 'input_faa', 'output_file',
                            'e_value_threshold', 'threads']

        for attr in required_attrs:
            if getattr(self, attr, None) is None:
                raise ValueError(f"Required attribute '{attr}' was not properly set")

        # Also validate file existence
        if not os.path.exists(self.diamond_db):
            raise FileNotFoundError(f"Database file not found: {self.diamond_db}")

        if not os.path.exists(self.input_faa):
            raise FileNotFoundError(f"Input file not found: {self.input_faa}")

        # Create output directory if it doesn't exist
        os.makedirs(os.path.dirname(os.path.abspath(self.output_file)), exist_ok=True)

    def _derive_diamond_db(self):
        """Derive DIAMOND database path - to be implemented by subclasses"""
        raise NotImplementedError("Subclasses must implement _derive_diamond_db()")

    def _derive_input_faa(self):
        """Derive input protein sequence file path - to be implemented by subclasses"""
        raise NotImplementedError("Subclasses must implement _derive_input_faa()")

    def _derive_output_file(self):
        """Derive output file path - to be implemented by subclasses"""
        raise NotImplementedError("Subclasses must implement _derive_output_file()")

    def _derive_e_value_threshold(self):
        """Derive E-value threshold - to be implemented by subclasses"""
        raise NotImplementedError("Subclasses must implement _derive_e_value_threshold()")

    def _derive_threads(self):
        """Derive number of threads to use"""
        return self.config.threads

    def _derive_verbose_option(self):
        """Derive verbose option flag"""
        return getattr(self.config, 'verbose_option', False)

    def run_diamond(self, outfmt=DIAMOND_DEFAULT_OUTFMT, extra_args=None):
        """Run DIAMOND BLASTP"""
        cmd = [DIAMOND_CMD,DIAMOND_BLASTP_CMD,
            DIAMOND_CMD_DB, self.diamond_db,
            DIAMOND_CMD_QUERY, self.input_faa,
            DIAMOND_CMD_OUT, self.output_file,
            DIAMOND_CMD_OUTFMT, outfmt,
            DIAMOND_CMD_EVALUE, str(self.e_value_threshold),
            DIAMOND_CMD_MAX_TARGET, str(DIAMOND_MAX_TARGET_SEQS),
            DIAMOND_CMD_THREADS, str(self.threads),
            DIAMOND_CMD_VERBOSE if self.verbose_option else DIAMOND_CMD_QUIET
        ]

        if extra_args:
            cmd.extend(extra_args)

        logging.info(f"Running DIAMOND BLASTp with {os.path.basename(self.diamond_db)}...")
        try:
            subprocess.run(cmd, check=True)
            logging.info("DIAMOND BLASTp completed")
        except subprocess.CalledProcessError as e:
            logging.error(f"DIAMOND BLASTp failed: {e}")
            raise

    def format_results(self, column_names, extra_processing=None):
        """Format results"""
        if not os.path.exists(self.output_file) or os.stat(self.output_file).st_size == 0:
            logging.warning(f"No results to format: {self.output_file} is empty or missing")
            return

        try:
            filtered_df = pd.read_csv(self.output_file, sep='\t', header=None, names=column_names)

            if extra_processing:
                extra_processing(filtered_df)

            filtered_df.to_csv(self.output_file, sep='\t', index=False)
            logging.info(f"Results formatted and saved to {self.output_file}")
        except Exception as e:
            logging.error(f"Error formatting results: {e}")
            raise


class CAZyDiamondProcessor(DiamondProcessor):
    """CAZyme DIAMOND processor"""

    def _derive_diamond_db(self):
        """Get CAZyme DIAMOND database path"""
        return os.path.join(self.config.db_dir, CAZY_DIAMOND_DB)

    def _derive_input_faa(self):
        """Get input protein sequence file path"""
        return os.path.join(self.config.output_dir, "uniInput.faa")

    def _derive_output_file(self):
        """Get output file path"""
        return os.path.join(self.config.output_dir, CAZY_DIAMOND_OUTPUT)

    def _derive_e_value_threshold(self):
        """Get E-value threshold for CAZyme searches"""
        # Use the value from config or default to DIAMOND_CAZY_EVALUE_DEFAULT
        return getattr(self.config, 'e_value_threshold', DIAMOND_CAZY_EVALUE_DEFAULT)

    def run(self):
        """Run CAZyme DIAMOND search"""
        self.run_diamond()

    def format_results(self):
        """Format CAZyme DIAMOND results"""
        super().format_results(CAZY_COLUMN_NAMES)


class TCDBDiamondProcessor(DiamondProcessor):
    """TCDB DIAMOND processor"""

    def _derive_diamond_db(self):
        """Get TCDB DIAMOND database path"""
        return os.path.join(self.config.db_dir, TCDB_DIAMOND_DB)

    def _derive_input_faa(self):
        """Get input protein sequence file path"""
        return os.path.join(self.config.output_dir, "uniInput.faa")

    def _derive_output_file(self):
        """Get output file path"""
        return os.path.join(self.config.output_dir, TCDB_DIAMOND_OUTPUT)

    def _derive_e_value_threshold(self):
        """Get E-value threshold for TCDB searches"""
        return getattr(self.config, 'e_value_threshold_tc', DIAMOND_TCDB_EVALUE_DEFAULT)

    def _derive_coverage_threshold(self):
        """Get coverage threshold for TCDB searches"""
        return getattr(self.config, 'coverage_threshold_tc', DIAMOND_TCDB_COVERAGE_DEFAULT)

    def run(self):
        """Run TCDB DIAMOND search"""
        # Get coverage threshold
        coverage_threshold = self._derive_coverage_threshold()

        # Set additional parameters
        extra_args = [
            DIAMOND_CMD_OUTFMT, DIAMOND_DEFAULT_OUTFMT, *TCDB_DIAMOND_OUTFMT_FIELDS,
            DIAMOND_CMD_QUERY_COVER, str(coverage_threshold)
        ]

        self.run_diamond(outfmt=DIAMOND_DEFAULT_OUTFMT, extra_args=extra_args)

    def format_results(self):
        """Format TCDB DIAMOND results"""
        def extra_processing(df):
            """Additional processing for TCDB results"""
            if TCDB_ID_COLUMN in df.columns:
                df[TCDB_ID_COLUMN] = df[TCDB_ID_COLUMN].apply(lambda x: x.split(' ')[0].split('|')[-1] if isinstance(x, str) else x)
            df['Database'] = TC

        super().format_results(TCDB_COLUMN_NAMES, extra_processing)


class SulfatlasDiamondProcessor(DiamondProcessor):
    """Sulfatlas DIAMOND processor"""

    def _derive_diamond_db(self):
        """Get Sulfatlas DIAMOND database path"""
        return os.path.join(self.config.db_dir, SULFATLAS_DIAMOND_DB)

    def _derive_input_faa(self):
        """Get input protein sequence file path"""
        return os.path.join(self.config.output_dir, "uniInput.faa")

    def _derive_output_file(self):
        """Get output file path"""
        return os.path.join(self.config.output_dir, SULFATLAS_DIAMOND_OUTPUT)

    def _derive_e_value_threshold(self):
        """Get E-value threshold for Sulfatlas searches"""
        return getattr(self.config, 'e_value_threshold_sulfatlas', DIAMOND_SULFATLAS_EVALUE_DEFAULT)

    def _derive_coverage_threshold(self):
        """Get coverage threshold for Sulfatlas searches"""
        return getattr(self.config, 'coverage_threshold_sulfatlas', DIAMOND_SULFATLAS_COVERAGE_DEFAULT)

    def run(self):
        """Run Sulfatlas DIAMOND search"""
        # Get coverage threshold
        coverage_threshold = self._derive_coverage_threshold()

        # Set additional parameters
        extra_args = [
            DIAMOND_CMD_OUTFMT, DIAMOND_DEFAULT_OUTFMT, *SULFATLAS_DIAMOND_OUTFMT_FIELDS,
            DIAMOND_CMD_QUERY_COVER, str(coverage_threshold)
        ]

        self.run_diamond(outfmt=DIAMOND_DEFAULT_OUTFMT, extra_args=extra_args)

    def format_results(self):
        """Format Sulfatlas DIAMOND results"""
        def extra_processing(df):
            """Additional processing for Sulfatlas results - extract S1_4 format IDs"""


            if SULFATLAS_ID_COLUMN in df.columns:
                df[SULFATLAS_ID_COLUMN] = df[SULFATLAS_ID_COLUMN].apply(
                    lambda x: "_".join(x.split('|')[1].split('_')[1:]) if isinstance(x, str) and '|' in x else "unknown"
                    )
            df['Database'] = SULFATLAS

        super().format_results(SULFATLAS_COLUMN_NAMES, extra_processing)


class PeptidaseDiamondProcessor(DiamondProcessor):
    """Peptidase DIAMOND processor"""

    def _derive_diamond_db(self):
        """Get Peptidase DIAMOND database path"""
        return os.path.join(self.config.db_dir, PEPTIDASE_DIAMOND_DB)

    def _derive_input_faa(self):
        """Get input protein sequence file path"""
        return os.path.join(self.config.output_dir, "uniInput.faa")

    def _derive_output_file(self):
        """Get output file path"""
        return os.path.join(self.config.output_dir, PEPTIDASE_DIAMOND_OUTPUT)

    def _derive_e_value_threshold(self):
        """Get E-value threshold for Peptidase searches"""
        return getattr(self.config, 'e_value_threshold_peptidase', DIAMOND_PEPTIDASE_EVALUE_DEFAULT)

    def _derive_coverage_threshold(self):
        """Get coverage threshold for Peptidase searches"""
        return getattr(self.config, 'coverage_threshold_peptidase', DIAMOND_PEPTIDASE_COVERAGE_DEFAULT)

    def run(self):
        """Run Peptidase DIAMOND search"""
        # Get coverage threshold
        coverage_threshold = self._derive_coverage_threshold()

        # Set additional parameters
        extra_args = [
            DIAMOND_CMD_OUTFMT, DIAMOND_DEFAULT_OUTFMT, *PEPTIDASE_DIAMOND_OUTFMT_FIELDS,
            DIAMOND_CMD_QUERY_COVER, str(coverage_threshold)
        ]

        self.run_diamond(outfmt=DIAMOND_DEFAULT_OUTFMT, extra_args=extra_args)

    def format_results(self):
        """Format Peptidase DIAMOND results"""
        def extra_processing(df):
            """Additional processing for Peptidase results - extract S01.001 format IDs"""
            import re
            # Extract S01.001 pattern from string like ">MER0000002|S01.001"
            if PEPTIDASE_ID_COLUMN in df.columns:
                df[PEPTIDASE_ID_COLUMN] = df[PEPTIDASE_ID_COLUMN].apply(
                lambda x: x.split('|')[1] if isinstance(x, str) and "|"  in x else "unknown"
                )

            df['Database'] = PEPTIDASE
        super().format_results(PEPTIDASE_COLUMN_NAMES, extra_processing)
