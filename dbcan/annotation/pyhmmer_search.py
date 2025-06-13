import pyhmmer
import logging
import psutil
import pandas as pd
import os
from dbcan.parameter import PyHMMERConfig, DBCANSUBProcessorConfig, PyHMMERSTPConfig, PyHMMERTFConfig
from dbcan.process.process_utils import process_results
from dbcan.process.process_dbcan_sub import DBCANSUBProcessor
from dbcan.constants import ( DBCAN_HMM_FILE, DBCAN_SUB_HMM_FILE,
                             TF_HMM_FILE, STP_HMM_FILE, INPUT_PROTEIN_NAME,
                             NON_CAZYME_PROTEIN_FILE, DBCAN_HMM_RESULT_FILE,
                             DBCAN_SUB_HMM_RESULT_FILE, TF_HMM_RESULT_FILE,
                             STP_HMM_RESULT_FILE, SUBSTRATE_MAPPING_FILE,
                             GT2_FAMILY_NAME, GT2_PREFIX)


logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class PyHMMERProcessor:
    """Base PyHMMER processor class using template method pattern"""

    def __init__(self, config):
        """Initialize with configuration and set up processor attributes"""
        self.config = config

        # Use template method pattern to set up all required attributes
        self._setup_processor()

    def _setup_processor(self):
        """Set up processor attributes - template method"""
        self.hmm_file = self._derive_hmm_file()
        self.input_faa = self._derive_input_faa()
        self.output_file = self._derive_output_file()
        self.e_value_threshold = self._derive_e_value_threshold()
        self.coverage_threshold = self._derive_coverage_threshold()
        self.hmmer_cpu = self._derive_threads()

        # Validate required attributes
        self._validate_attributes()

    def _validate_attributes(self):
        """Validate that all required attributes are properly set"""
        required_attrs = ['hmm_file', 'input_faa', 'output_file',
                        'e_value_threshold', 'coverage_threshold', 'hmmer_cpu']

        for attr in required_attrs:
            if getattr(self, attr, None) is None:
                raise ValueError(f"Required attribute '{attr}' was not properly set")

    def _derive_hmm_file(self):
        """Derive HMM file path - to be implemented by subclasses"""
        # Abstract method - force subclasses to implement
        raise NotImplementedError("Subclasses must implement _derive_hmm_file()")

    def _derive_input_faa(self):
        """Derive input protein sequence file path - to be implemented by subclasses"""
        # Abstract method - force subclasses to implement
        raise NotImplementedError("Subclasses must implement _derive_input_faa()")

    def _derive_output_file(self):
        """Derive output file path - to be implemented by subclasses"""
        # Abstract method - force subclasses to implement
        raise NotImplementedError("Subclasses must implement _derive_output_file()")

    def _derive_e_value_threshold(self):
        """Derive E-value threshold - to be implemented by subclasses"""
        # Abstract method - force subclasses to implement
        raise NotImplementedError("Subclasses must implement _derive_e_value_threshold()")

    def _derive_coverage_threshold(self):
        """Derive coverage threshold - to be implemented by subclasses"""
        # Abstract method - force subclasses to implement
        raise NotImplementedError("Subclasses must implement _derive_coverage_threshold()")

    def _derive_threads(self):
        """Derive number of threads to use"""
        # Default implementation uses the threads from config
        return self.config.threads

    def hmmsearch(self):
        """Execute HMMER search"""
        available_memory = psutil.virtual_memory().available
        target_size = os.stat(self.input_faa).st_size
        hmm_files = pyhmmer.plan7.HMMFile(self.hmm_file)
        results = []

        with pyhmmer.easel.SequenceFile(self.input_faa, digital=True) as seqs:
            targets = seqs.read_block() if target_size < available_memory * 0.1 else seqs
            for hits in pyhmmer.hmmsearch(hmm_files, targets, cpus=self.hmmer_cpu, domE=self.e_value_threshold):
                for hit in hits:
                    for domain in hit.domains.included:
                        coverage = (domain.alignment.hmm_to - domain.alignment.hmm_from + 1) / domain.alignment.hmm_length
                        hmm_name = domain.alignment.hmm_name.decode('utf-8')
                        if GT2_PREFIX in hmm_name:
                            hmm_name = GT2_FAMILY_NAME
                        hmm_length = domain.alignment.hmm_length
                        target_name = domain.alignment.target_name.decode('utf-8')
                        target_length = domain.alignment.target_length
                        i_evalue = domain.i_evalue
                        hmm_from = domain.alignment.hmm_from
                        hmm_to = domain.alignment.hmm_to
                        target_from = domain.alignment.target_from
                        target_to = domain.alignment.target_to
                        hmm_file_name = hmm_files.name.split("/")[-1].split(".")[0]
                        if i_evalue < self.e_value_threshold and coverage > self.coverage_threshold:
                            results.append([hmm_name, hmm_length, target_name, target_length, i_evalue, hmm_from, hmm_to, target_from, target_to, coverage, hmm_file_name])

        logging.info(f"{self.hmm_file} PyHMMER search completed. Found {len(results)} hits.")
        process_results(results, self.output_file)


class PyHMMERDBCANProcessor(PyHMMERProcessor):
    """dbCAN HMM processor"""

    def _derive_hmm_file(self):
        """Get dbCAN HMM file path"""
        return os.path.join(self.config.db_dir, DBCAN_HMM_FILE)

    def _derive_input_faa(self):
        """Get input protein sequence file path"""
        return os.path.join(self.config.output_dir, INPUT_PROTEIN_NAME)

    def _derive_output_file(self):
        """Get output file path"""
        return os.path.join(self.config.output_dir, DBCAN_HMM_RESULT_FILE)

    def _derive_e_value_threshold(self):
        """Get E-value threshold for dbCAN searches"""
        return self.config.e_value_threshold_dbcan

    def _derive_coverage_threshold(self):
        """Get coverage threshold for dbCAN searches"""
        return self.config.coverage_threshold_dbcan

    def _derive_threads(self):
        """Derive number of threads to use"""
        # Default implementation uses the threads from config
        return self.config.threads

    def run(self):
        """Run dbCAN HMM search"""
        self.hmmsearch()


class PyHMMERDBCANSUBProcessor(PyHMMERProcessor):
    """dbCAN-sub HMM processor"""

    def __init__(self, config):
        """Initialize with configuration"""
        super().__init__(config)
        # Set additional attributes specific to dbCAN-sub
        self.mapping_file = self._derive_mapping_file()

    def _derive_hmm_file(self):
        """Get dbCAN-sub HMM file path"""
        return os.path.join(self.config.db_dir, DBCAN_SUB_HMM_FILE)

    def _derive_input_faa(self):
        """Get input protein sequence file path"""
        return os.path.join(self.config.output_dir, INPUT_PROTEIN_NAME)

    def _derive_output_file(self):
        """Get output file path"""
        return os.path.join(self.config.output_dir, DBCAN_SUB_HMM_RESULT_FILE)

    def _derive_e_value_threshold(self):
        """Get E-value threshold for dbCAN-sub searches"""
        return self.config.e_value_threshold_dbsub

    def _derive_coverage_threshold(self):
        """Get coverage threshold for dbCAN-sub searches"""
        return self.config.coverage_threshold_dbsub

    def _derive_mapping_file(self):
        """Get substrate mapping file path"""
        return os.path.join(self.config.db_dir, SUBSTRATE_MAPPING_FILE)

    def _derive_threads(self):
        """Derive number of threads to use"""
        # Default implementation uses the threads from config
        return self.config.threads

    def run(self):
        """Run dbCAN-sub HMM search and process results"""
        self.hmmsearch()
        dbcansub_processor = DBCANSUBProcessor(self.config)
        dbcansub_processor.process_dbcan_sub()


class PyHMMERTFProcessor(PyHMMERProcessor):
    """Transcription Factor HMM processor"""

    def _derive_hmm_file(self):
        """Get TF HMM file path"""
        return os.path.join(self.config.db_dir, TF_HMM_FILE)

    def _derive_input_faa(self):
        """Get input protein sequence file path"""
        return os.path.join(self.config.output_dir, INPUT_PROTEIN_NAME)

    def _derive_output_file(self):
        """Get output file path"""
        return os.path.join(self.config.output_dir, TF_HMM_RESULT_FILE)

    def _derive_e_value_threshold(self):
        """Get E-value threshold for TF searches"""
        return self.config.e_value_threshold_tf

    def _derive_coverage_threshold(self):
        """Get coverage threshold for TF searches"""
        return self.config.coverage_threshold_tf

    def _derive_threads(self):
        """Derive number of threads to use"""
        # Default implementation uses the threads from config
        return self.config.threads

    def run(self):
        """Run TF HMM search"""
        self.hmmsearch()


class PyHMMERSTPProcessor(PyHMMERProcessor):
    """Signal Transduction Protein HMM processor"""

    def _derive_hmm_file(self):
        """Get STP HMM file path"""
        return os.path.join(self.config.db_dir, STP_HMM_FILE)

    def _derive_input_faa(self):
        """Get input protein sequence file path"""
        return os.path.join(self.config.output_dir, INPUT_PROTEIN_NAME)

    def _derive_output_file(self):
        """Get output file path"""
        return os.path.join(self.config.output_dir, STP_HMM_RESULT_FILE)

    def _derive_e_value_threshold(self):
        """Get E-value threshold for STP searches"""
        return self.config.e_value_threshold_stp

    def _derive_coverage_threshold(self):
        """Get coverage threshold for STP searches"""
        return self.config.coverage_threshold_stp

    def _derive_threads(self):
        """Derive number of threads to use"""
        # Default implementation uses the threads from config
        return self.config.threads

    def run(self):
        """Run STP HMM search"""
        self.hmmsearch()
