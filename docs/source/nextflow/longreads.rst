.. _longreads-mode:

Long Reads Analysis Mode
=========================

The long reads analysis mode (``--type longreads``) is designed for processing PacBio or Nanopore long-read sequencing data. This mode uses Flye for metagenomic assembly, which is optimized for long-read technologies.

Overview
--------

The long reads mode is optimized for third-generation sequencing technologies (PacBio and Nanopore) that produce longer reads compared to Illumina short-read sequencing. Longer reads enable better assembly of complex metagenomic communities and improved detection of complete gene clusters.

Workflow
--------

The long reads workflow consists of the following main steps:

1. **Quality Control** (FastQC, optional for DNA)
   - FastQC quality assessment (primarily for RNA-seq data)
   - DNA long reads typically skip QC/trimming steps as they are less affected by adapters

2. **Taxonomic Filtering** (Kraken2, optional)
   - Taxonomic classification using Kraken2 (applied to RNA-seq data)
   - Extraction of reads matching specified taxonomy
   - Can be skipped with ``--skip_kraken_extraction``

3. **Assembly** (Flye)
   - Long-read metagenomic assembly using Flye
   - Supports multiple Flye modes for different sequencing technologies
   - Configurable via ``--flye_mode`` parameter

4. **Gene Prediction** (Pyrodigal)
   - Prodigal-based gene finding optimized for metagenomic data
   - Generates protein sequences (FAA) and gene annotations (GFF)

5. **CAZyme Annotation** (run_dbCAN)
   - CAZyme identification using dbCAN database
   - CGC (CAZyme Gene Cluster) detection
   - Substrate prediction

6. **Read Mapping** (Minimap2 for DNA, BWA-MEM for RNA)
   - Mapping of long DNA reads back to assembled contigs using Minimap2
   - Mapping of RNA reads (if provided) using BWA-MEM for expression analysis
   - Coverage calculation for genes and CGCs

7. **Abundance Calculation**
   - Gene-level abundance calculation based on read coverage
   - CGC abundance and visualization
   - Generation of bar plots and heatmaps

8. **Report Generation** (MultiQC)
   - Aggregated quality control and analysis reports

Usage
-----

.. note::
   Before running the pipeline, make sure you have cloned the repository. See :ref:`nextflow-usage` for installation instructions.

Basic Usage
~~~~~~~~~~~

The simplest command to run long reads analysis (assuming you are in the ``dbcan-nf`` directory):

.. code-block:: bash

   nextflow run main.nf \
     --input samplesheet.csv \
     --outdir results \
     --type longreads \
     -profile docker

Flye Mode Selection
~~~~~~~~~~~~~~~~~~~

The ``--flye_mode`` parameter allows you to specify the appropriate Flye mode for your sequencing technology:

.. code-block:: bash

   # PacBio HiFi reads (default)
   nextflow run main.nf \
     --input samplesheet.csv \
     --outdir results \
     --type longreads \
     --flye_mode --pacbio-hifi \
     -profile docker

   # PacBio raw reads
   nextflow run main.nf \
     --input samplesheet.csv \
     --outdir results \
     --type longreads \
     --flye_mode --pacbio-raw \
     -profile docker

   # Nanopore raw reads
   nextflow run main.nf \
     --input samplesheet.csv \
     --outdir results \
     --type longreads \
     --flye_mode --nano-raw \
     -profile docker

   # Nanopore high-quality reads
   nextflow run main.nf \
     --input samplesheet.csv \
     --outdir results \
     --type longreads \
     --flye_mode --nano-hq \
     -profile docker

With RNA-seq Data
~~~~~~~~~~~~~~~~~

When RNA-seq transcriptome data is provided in the samplesheet, the pipeline will:

- Process RNA reads through quality control (FastQC + TrimGalore)
- Apply Kraken2 taxonomic filtering (if enabled)
- Map RNA reads to assembled contigs using BWA-MEM
- Calculate RNA-based abundance for expression analysis
- Generate separate DNA and RNA abundance plots

Example with RNA-seq:

.. code-block:: bash

   nextflow run main.nf \
     --input samplesheet_with_rna.csv \
     --outdir results \
     --type longreads \
     --flye_mode --pacbio-hifi \
     -profile docker

Skipping Steps
~~~~~~~~~~~~~~

You can skip certain steps if needed:

.. code-block:: bash

   nextflow run main.nf \
     --input samplesheet.csv \
     --outdir results \
     --type longreads \
     --skip_kraken_extraction \
     -profile docker

Flye Mode Options
-----------------

The ``--flye_mode`` parameter accepts the following Flye assembly modes:

.. list-table:: Flye Mode Options
   :widths: 20 80
   :header-rows: 1

   * - Mode
     - Description
   * - ``--pacbio-hifi``
     - PacBio HiFi (high-fidelity) reads. Default mode. Best for high-quality PacBio data.
   * - ``--pacbio-raw``
     - PacBio raw (CLR) reads. Use for standard PacBio sequencing data.
   * - ``--nano-raw``
     - Nanopore raw reads. Use for standard Nanopore sequencing data.
   * - ``--nano-hq``
     - Nanopore high-quality reads. Use for Q20+ or similar high-quality Nanopore data.

Output Files
------------

The pipeline generates output files organized in the following directory structure:

Assembly Results
~~~~~~~~~~~~~~~~

- ``flye/``
  - ``*_assembly.fasta.gz``: Assembled contigs in FASTA format (gzipped)
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

- ``minimap2/``: Minimap2 alignment files for long DNA reads
  - ``*.bam``: Aligned reads
  - ``*.bam.bai``: BAM index files

- ``bwa_index_mem/``: BWA-MEM alignment files for RNA reads (if provided)
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

- ``fastqc/``: FastQC reports (for RNA-seq data)
- ``trimgalore/``: TrimGalore reports (for RNA-seq data)
- ``multiqc/``: MultiQC aggregated report

Pipeline Information
~~~~~~~~~~~~~~~~~~~~

- ``pipeline_info/``: Execution reports, parameters, and software versions

Key Features
------------

- **Long-Read Optimized**: Uses Flye assembler optimized for PacBio and Nanopore reads
- **Flexible Flye Modes**: Supports multiple Flye modes for different sequencing technologies
- **Dual Mapping**: Minimap2 for long DNA reads, BWA-MEM for RNA reads
- **Complete Gene Clusters**: Longer reads enable better assembly of complete CGCs
- **RNA-seq Integration**: Optional RNA-seq support for expression analysis

Advantages of Long Reads
------------------------

- **Better Contiguity**: Longer contigs due to long-read sequencing
- **Complete Genes**: More complete gene predictions, especially for large genes
- **Reduced Fragmentation**: Fewer fragmented gene clusters
- **Repeat Resolution**: Better resolution of repetitive regions
- **Structural Variants**: Improved detection of structural variations

Considerations
--------------

- **Computational Requirements**: Long-read assembly typically requires more memory and time
- **Error Rates**: Long reads may have higher error rates than short reads
- **Coverage**: Lower coverage requirements compared to short reads
- **Cost**: Long-read sequencing may be more expensive per base

Best Practices
--------------

1. **Choose Correct Mode**: Select the appropriate ``--flye_mode`` for your sequencing technology
2. **Quality Assessment**: Review FastQC reports for RNA-seq data
3. **Coverage Planning**: Long reads require less coverage than short reads (typically 20-30x)
4. **Resource Allocation**: Ensure sufficient memory for Flye assembly

When to Use Long Reads Mode
----------------------------

**Recommended**:
- PacBio or Nanopore sequencing data
- When complete gene clusters are important
- When dealing with complex or repetitive regions
- When structural variation detection is needed

**Not Recommended**:
- Illumina short-read data (use ``--type shortreads`` instead)
- When computational resources are very limited
- When only basic CAZyme annotation is needed

Example Results
---------------

For example visualizations from long reads mode, see :ref:`nextflow-results-examples`.

See Also
--------

- :ref:`shortreads-mode` - Short reads mode documentation
- :ref:`assemfree-mode` - Assembly-free mode documentation
- :ref:`nextflow-parameters` - Complete parameter reference
- :ref:`nextflow-results-examples` - Example results and visualizations
- `Flye Documentation <https://github.com/fenderglass/Flye>`_ - Flye assembler documentation
