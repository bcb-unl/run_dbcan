import pandas as pd
import os
import logging
from dbcan.parameter import DBCANSUBProcessorConfig
from dbcan.constants import (DBCAN_SUB_COLUMN_NAMES, DBCAN_SUB_RESULT_FILE,
                            SUBSTRATE_MAPPING_FILE,
                            DBCAN_SUB_HMM_NAME_COLUMN, DBCAN_SUB_SUBFAMILY_NAME_COLUMN,
                            DBCAN_SUB_SUBFAMILY_COMP_COLUMN, DBCAN_SUB_SUBFAMILY_EC_COLUMN,
                            DBCAN_SUB_SUBSTRATE_COLUMN,
                            DBCAN_SUB_HMM_SUFFIX, DBCAN_SUB_SEPARATOR
                            )

class DBCANSUBProcessor:
    """Process dbCAN-sub results using template method pattern"""

    def __init__(self, config: DBCANSUBProcessorConfig):
        """Initialize with configuration"""
        self.config = config

        # Use template method pattern to set up all required attributes
        self._setup_processor()

    def _setup_processor(self):
        """Set up processor attributes using template method pattern"""
        self.output_file = self._derive_output_file()
        self.mapping_file = self._derive_mapping_file()

        # Validate required attributes
        self._validate_attributes()

    def _validate_attributes(self):
        """Validate that all required attributes are properly set"""
        required_attrs = ['output_file', 'mapping_file']

        for attr in required_attrs:
            if getattr(self, attr, None) is None:
                raise ValueError(f"Required attribute '{attr}' was not properly set")

        # Check if mapping file exists
        if not os.path.exists(self.mapping_file):
            raise FileNotFoundError(f"Substrate mapping file not found: {self.mapping_file}")

        # Check if output file exists (may not exist if hmmsearch failed)
        if not os.path.exists(self.output_file):
            logging.warning(f"dbCAN-sub results file not found: {self.output_file}. No substrate processing will be performed.")

    def _derive_output_file(self):
        """Derive output file path"""
        return getattr(self.config, 'output_file', None) or os.path.join(self.config.output_dir, DBCAN_SUB_RESULT_FILE)

    def _derive_mapping_file(self):
        """Derive mapping file path"""
        return getattr(self.config, 'mapping_file', None) or os.path.join(self.config.db_dir, SUBSTRATE_MAPPING_FILE)

    def load_substrate_mapping(self):
        """
        Load substrate mapping from file

        Example of mapping file:
        Substrate_high_level	 Substrate_curated	Family	Name	        EC_Number
        lignin	                 lignin	            AA1	    ferroxidase	    1.10.3.2
        chitin                    chitin            CBM14   Long description    NA
        """
        try:
            df = pd.read_csv(self.mapping_file, sep='\t', header=None, skiprows=1, usecols=[2, 4, 0])
            df[4] = df[4].str.strip().fillna('-')
            df['key'] = df.apply(lambda x: (x[2], x[4]), axis=1)
            return pd.Series(df[0].values, index=pd.MultiIndex.from_tuples(df['key'])).to_dict()
        except FileNotFoundError:
            logging.error(f"Can't find substrate mapping file: {self.mapping_file}")
            return {}
        except Exception as e:
            logging.error(f"Error loading substrate mapping: {e}")
            return {}

    def process_dbcan_sub(self):
        """Process dbCAN-sub results to add substrate information"""
        # Skip if output file doesn't exist
        if not os.path.exists(self.output_file):
            logging.warning(f"Output file not found: {self.output_file}. Skipping substrate processing.")
            return

        # Skip if file is empty
        if os.path.getsize(self.output_file) == 0:
            logging.warning(f"Output file is empty: {self.output_file}. Skipping substrate processing.")
            return

        subs_dict = self.load_substrate_mapping()
        if not subs_dict:
            logging.warning("No substrate mapping data loaded. Substrate annotation will be empty.")

        try:
            df = pd.read_csv(self.output_file, sep='\t')
            if df.empty:
                logging.warning("No dbCAN-sub results to process")
                return

            logging.info(f"Processing {len(df)} dbCAN-sub results")

            #  extract information from HMM Name : PL25_e0.hmm|PL25:38|PL0:1|3.2.1.122:13
            df[DBCAN_SUB_SUBFAMILY_NAME_COLUMN] = df[DBCAN_SUB_HMM_NAME_COLUMN].apply(
                lambda x: DBCAN_SUB_SEPARATOR.join(p.split('.')[0] for p in x.split(DBCAN_SUB_SEPARATOR) if DBCAN_SUB_HMM_SUFFIX in p))
            df[DBCAN_SUB_SUBFAMILY_COMP_COLUMN] = df[DBCAN_SUB_HMM_NAME_COLUMN].apply(
                lambda x: DBCAN_SUB_SEPARATOR.join(p for p in x.split(DBCAN_SUB_SEPARATOR) if DBCAN_SUB_HMM_SUFFIX not in p and len(p.split('.')) != 4))
            df[DBCAN_SUB_SUBFAMILY_EC_COLUMN] = df[DBCAN_SUB_HMM_NAME_COLUMN].apply(
                lambda x: DBCAN_SUB_SEPARATOR.join(p for p in x.split(DBCAN_SUB_SEPARATOR) if len(p.split('.')) == 4))
            df[DBCAN_SUB_SUBSTRATE_COLUMN] = df[DBCAN_SUB_HMM_NAME_COLUMN].apply(
                lambda x: self.get_substrates(x, subs_dict))
            df.drop(DBCAN_SUB_HMM_NAME_COLUMN, axis=1, inplace=True)

            # Ensure columns are in the right order
            for col in DBCAN_SUB_COLUMN_NAMES:
                if col not in df.columns:
                    df[col] = '-'  # Add missing columns with default value
            df = df[DBCAN_SUB_COLUMN_NAMES]

            df.to_csv(self.output_file, sep='\t', index=False)
            logging.info(f"Successfully processed dbCAN-sub results with substrate information ({len(df)} entries)")
        except Exception as e:
            logging.error(f"Error processing dbCAN-sub results: {str(e)}")
            import traceback
            traceback.print_exc()

    def get_substrates(self, profile_info, subs_dict):
        """Extract substrate information from profile info"""
        if not profile_info or not isinstance(profile_info, str):
            return '-'

        parts = profile_info.split('|')
        substrates = set()

        # Handle empty parts
        if not parts:
            return '-'

        # Get the family name (before .hmm)
        try:
            key1 = parts[0].split('.hmm')[0].split("_")[0]
        except (IndexError, AttributeError):
            return '-'

        # Process CBM families
        if key1.startswith('CBM'):
            key2 = '-'
            if (key1, key2) in subs_dict:
                substrates.add(subs_dict[(key1, key2)])
        else:
            # Process non-CBM families
            for p in parts:
                if ':' in p and '.' in p.split(':')[0] and len(p.split(':')[0].split('.')) == 4:
                    key2 = p.split(':')[0]
                    if (key1, key2) in subs_dict:
                        substrates.add(subs_dict[(key1, key2)])

            # If no substrates found, try with default EC
            if not substrates and not key1.startswith('CBM'):
                if (key1, '-') in subs_dict:
                    substrates.add(subs_dict[(key1, '-')])

        return ';'.join(sorted(substrates)) if substrates else '-'
