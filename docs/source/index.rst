User Guide
==========

Update: What's New in run_dbCAN
================================

The new version of **run_dbCAN** introduces multiple new features and significant performance improvements, making the pipeline more user-friendly and efficient. We highly recommend users to upgrade to this version.
If you have any questions or suggestions, please feel free to contact us:

* Dr. Yanbin Yin, Professor (yyin@unl.edu)
* Xinpeng Zhang, PhD Student (xzhang55@huskers.unl.edu)
* Dr. Haidong Yi, Software Engineer  (hyi@stjude.org)

All conda environments dependencies can be found at the following link:
`run_dbCAN Conda Environments <https://github.com/bcb-unl/run_dbcan_new/tree/master/envs>`_


Key Features and Improvements
-----------------------------

1. **Simplified Database Downloading**

   - Added a new function for downloading database files, making the process simpler than before.

2. **Enhanced Input Processing**

   - Replaced `prodigal` with `pyrodigal <https://pyrodigal.readthedocs.io/en/stable/>`_ for input processing.

3. **Improved HMMER Performance**

   - Replaced `HMMER` with `pyHMMER <https://pyhmmer.readthedocs.io/en/stable/>`_, which is faster and more efficient.

   - Redesigned memory usage to support both low-memory and high-efficiency modes (https://pyhmmer.readthedocs.io/en/stable/examples/performance_tips.html).

4. **Modular Code Structure**

   - Reorganized the logic and structure of `run_dbCAN` by splitting functions into modules and following Object Oriented Programming.

   - Rewrote non-Python code in Python for improved readability.

   - Centralized parameter management using configuration files.

   - Leveraged the power of `pandas` for efficient data processing.

   - Added extensive logging and time reporting to make the pipeline more user-friendly.

5. **Enhanced dbCAN-sub and overview Features**

   - Added coverage justifications and location information for dbCAN-sub.

   - Included CAZyme justification in the final results with an extra column called "Best Results."

   - Now follow the rule: `CAZy-sub > dbCAN-sub > dbCAN-fam` for the final results.

6. **Redesigned CGCFinder**

   - Now supports JGI, NCBI, and Prodigal gff formats.

   - Directly searches eukaryotic genomes, including fungi (**beta function**).

   - Added a new function to visualize the CGCs on the genome (**beta function**).

7.  **Faster Substrate Prediction**

   - Replaced `blastp` with `DIAMOND` for substrate prediction, significantly improving speed and efficiency.

8.  **Updated Metagenomic Protocols**

   - Improved steps for metagenomic data processing (https://www.biorxiv.org/content/10.1101/2024.01.10.575125v1).

.. hint::

   If you want to run the pipeline from raw metagenomic reads, please refer to the following part:
   **metagenomics_pipeline**


Otherwise, refer to the instructions below. Please note that some precomputed results may have different names compared to the previous version.

.. note::
   For detailed instructions, refer to the respective sections in the documentation.

.. toctree::
   :maxdepth: 1
   :caption: Getting start

   getting_started/installation
   getting_started/database_description
   getting_started/quick_start



.. toctree::
   :maxdepth: 1
   :caption: User guide

   user_guide/prepare_the_database
   user_guide/CAZyme_annotation
   user_guide/CGC_information_generation
   user_guide/CGC_annotation
   user_guide/predict_CGC_substrate
   user_guide/CGC_plots

.. toctree::
   :maxdepth: 1
   :caption: Comparison

   comparison/CAZyme_annotation_compare


.. toctree::
   :maxdepth: 1
   :caption: API

   api/index

.. toctree::
   :maxdepth: 1
   :caption: Metagenomics pipeline

   metagenomics_pipeline/run_from_raw_reads
   metagenomics_pipeline/run_from_raw_reads_am
   metagenomics_pipeline/run_from_raw_reads_pr
   metagenomics_pipeline/run_from_raw_reads_wk
   metagenomics_pipeline/run_from_raw_reads_em
   metagenomics_pipeline/supplement/run_from_raw_reads_sp_co_assem
   metagenomics_pipeline/supplement/run_from_raw_reads_sp_subsample
   metagenomics_pipeline/supplement/run_from_raw_reads_sp_assem_free

.. toctree::
   :maxdepth: 1
   :caption: Change logs

   change-logs/index



.. toctree::
   :maxdepth: 1
   :caption: References

   references/index


.. toctree::
   :maxdepth: 1
   :caption: Contributors

   contributors/index
