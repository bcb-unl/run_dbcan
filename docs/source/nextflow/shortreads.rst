.. _shortreads-mode:

Short Reads Analysis Mode
==========================

The short reads analysis mode (``--type shortreads``) is designed for processing Illumina short-read sequencing data. This mode performs assembly-based CAZyme annotation using MEGAHIT for metagenomic assembly, followed by gene prediction, CAZyme annotation, and abundance calculation.

Overview
--------

The short reads mode is the default analysis mode and is optimized for Illumina paired-end or single-end sequencing data. It provides comprehensive CAZyme and CGC (CAZyme Gene Cluster) analysis with optional RNA-seq integration for expression analysis.

Workflow
--------

The short reads workflow consists of the following main steps:

1. **Quality Control** (FastQC + TrimGalore)
   - FastQC quality assessment of raw sequencing reads
   - TrimGalore adapter trimming and quality filtering

2. **Taxonomic Filtering** (Kraken2, optional)
   - Taxonomic classification using Kraken2
   - Extraction of reads matching specified taxonomy (default: human reads, tax ID 9606)
   - Can be skipped with ``--skip_kraken_extraction``

3. **Read Processing** (optional)
   - **Subsampling**: Downsample reads before assembly (see :ref:`shortreads-subsample`)
   - **Co-assembly**: Combine all samples for joint assembly (see :ref:`shortreads-coassembly`)

4. **Assembly** (MEGAHIT)
   - Metagenomic assembly using MEGAHIT
   - Default minimum contig length: 1000 bp
   - Memory usage limited to 50% of available memory

5. **Gene Prediction** (Pyrodigal)
   - Prodigal-based gene finding optimized for metagenomic data
   - Generates protein sequences (FAA) and gene annotations (GFF)

6. **CAZyme Annotation** (run_dbCAN)
   - CAZyme identification using dbCAN database
   - CGC (CAZyme Gene Cluster) detection
   - Substrate prediction

7. **Read Mapping** (BWA-MEM)
   - Mapping of DNA reads back to assembled contigs
   - Mapping of RNA reads (if provided) for expression analysis
   - Coverage calculation for genes and CGCs

8. **Abundance Calculation**
   - Gene-level abundance calculation based on read coverage
   - CGC abundance and visualization
   - Generation of bar plots and heatmaps

9. **Report Generation** (MultiQC)
   - Aggregated quality control and analysis reports

Usage
-----

Basic Usage
~~~~~~~~~~~

The simplest command to run short reads analysis:

.. code-block:: bash

   nextflow run nf-core/dbcanmicrobiome \
     --input samplesheet.csv \
     --outdir results \
     --type shortreads \
     -profile docker

With RNA-seq Data
~~~~~~~~~~~~~~~~~

When RNA-seq transcriptome data is provided in the samplesheet, the pipeline will automatically:

- Process RNA reads through quality control
- Map RNA reads to assembled contigs
- Calculate RNA-based abundance for expression analysis
- Generate separate DNA and RNA abundance plots

.. note::
   RNA-seq processing is automatically disabled when using ``--subsample`` or ``--coassembly`` modes.

Example with RNA-seq:

.. code-block:: bash

   nextflow run nf-core/dbcanmicrobiome \
     --input samplesheet_with_rna.csv \
     --outdir results \
     --type shortreads \
     -profile docker

Skipping Steps
~~~~~~~~~~~~~~

You can skip certain steps if needed:

.. code-block:: bash

   nextflow run nf-core/dbcanmicrobiome \
     --input samplesheet.csv \
     --outdir results \
     --type shortreads \
     --skip_fastqc \
     --skip_trimming \
     --skip_kraken_extraction \
     -profile docker

Advanced Options
~~~~~~~~~~~~~~~~

Using custom databases:

.. code-block:: bash

   nextflow run nf-core/dbcanmicrobiome \
     --input samplesheet.csv \
     --outdir results \
     --type shortreads \
     --dbcan_db /path/to/dbcan_database \
     --kraken_db /path/to/kraken_database \
     -profile docker

Subsampling and Co-assembly
---------------------------

The short reads mode supports two special processing options:

- :ref:`Subsampling Mode <shortreads-subsample>`: Downsample reads before assembly to reduce computational requirements
- :ref:`Co-assembly Mode <shortreads-coassembly>`: Combine all samples for joint assembly

These modes are mutually exclusive and cannot be used together.

Output Files
------------

The pipeline generates output files organized in the following directory structure:

Assembly Results
~~~~~~~~~~~~~~~~

- ``megahit/``
  - ``*_contigs.fa.gz``: Assembled contigs in FASTA format (gzipped)
  - Assembly statistics and reports

Gene Prediction
~~~~~~~~~~~~~~~

- ``pyrodigal/``
  - ``*.faa.gz``: Predicted protein sequences (gzipped)
  - ``*.gff.gz``: Gene annotations in GFF format (gzipped)

CAZyme Annotation
~~~~~~~~~~~~~~~~~

- ``rundbcan/``
  - ``*_dbcan/``: Directory containing all run_dbCAN results
  
    Files in this directory:
    
    - ``overview.tsv``: CAZyme annotation overview
    - ``dbCAN_hmm_results.tsv``: HMM-based CAZyme predictions
    - ``dbCANsub_hmm_results.tsv``: Subfamily predictions
    - ``diamond.out``: DIAMOND search results
    - ``*_cgc.gff``: CGC annotations
    - ``*_cgc_standard_out.tsv``: CGC standard output
    - ``*_substrate_prediction.tsv``: Substrate predictions
    - ``*_synteny_pdf/``: Synteny plots for CGCs

Read Mapping
~~~~~~~~~~~~

- ``bwa/``: BWA index files
- ``bwa_index_mem/``: BAM files from read mapping
  - ``*.bam``: Aligned reads
  - ``*.bam.bai``: BAM index files

Coverage and Abundance
~~~~~~~~~~~~~~~~~~~~~~~

- ``dbcan_utils_cal_coverage/``
  - ``*_depth.txt``: Gene coverage depth files

- ``dbcan_utils_cal_abund/``
  - ``*_abund/``: Abundance calculation results
  
    Files in this directory:
    
    - ``*_abund.txt``: Gene abundance values
    - ``*_cgc_abund.txt``: CGC abundance values

- ``dbcan_plot/``
  - ``*_pdf/``: Visualization plots
  
    Files in this directory:
    
    - ``heatmap.pdf``: Abundance heatmap
    - ``ec.pdf``: EC number distribution
    - ``family.pdf``: CAZyme family distribution
    - ``subfamily.pdf``: Subfamily distribution

- ``cgc_depth_plot/``
  - ``*_cgc_depth.tsv``: CGC depth coverage data
  - ``*_cgc_depth.pdf``: CGC depth plots

Quality Control
~~~~~~~~~~~~~~~

- ``fastqc/``: FastQC reports (if not skipped)
- ``trimgalore/``: TrimGalore reports (if not skipped)
- ``multiqc/``: MultiQC aggregated report

Pipeline Information
~~~~~~~~~~~~~~~~~~~~

- ``pipeline_info/``: Execution reports, parameters, and software versions

Key Features
------------

- **Dual Analysis**: Supports both DNA and RNA-seq data for comprehensive analysis
- **Flexible Processing**: Optional subsampling and co-assembly modes
- **Comprehensive Annotation**: CAZyme identification, CGC detection, and substrate prediction
- **Abundance Calculation**: Gene-level and CGC-level abundance with visualization
- **Quality Control**: Integrated QC pipeline with MultiQC reporting

See Also
--------

- :ref:`shortreads-subsample` - Subsampling mode details
- :ref:`shortreads-coassembly` - Co-assembly mode details
- :ref:`nextflow-parameters` - Complete parameter reference
- :ref:`nextflow-output` - General output documentation
