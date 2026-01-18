.. _assemfree-mode:

Assembly-Free Analysis Mode
============================

The assembly-free analysis mode (``--type assemfree``) performs direct CAZyme annotation on sequencing reads without assembly. This mode is particularly useful for large datasets where assembly is computationally expensive or when you want to avoid potential assembly artifacts.

Overview
--------

The assembly-free mode bypasses the assembly step and directly searches reads against the CAZyme database using DIAMOND blastx. This approach:

- **Reduces Computational Cost**: No assembly step means lower memory and CPU requirements
- **Faster Processing**: Direct read annotation is faster than assembly-based workflows
- **Avoids Assembly Artifacts**: Eliminates potential biases introduced during assembly
- **Scalable**: Better suited for very large datasets

This mode is ideal when you have large metagenomic datasets and want quick CAZyme annotation without the overhead of assembly.

Workflow
--------

The assembly-free workflow consists of the following main steps:

1. **Quality Control** (FastQC + TrimGalore)
   - FastQC quality assessment of raw sequencing reads
   - TrimGalore adapter trimming and quality filtering

2. **Taxonomic Filtering** (Kraken2, optional)
   - Taxonomic classification using Kraken2
   - Extraction of reads matching specified taxonomy (default: human reads, tax ID 9606)
   - Can be skipped with ``--skip_kraken_extraction``

3. **Read Conversion** (seqtk seq)
   - Convert FASTQ reads to FASTA format
   - Handles both paired-end and single-end reads

4. **CAZyme Search** (DIAMOND blastx)
   - Direct alignment of reads against CAZy DIAMOND database
   - Fast protein search using DIAMOND's optimized algorithm
   - Automatic download of CAZy database if not provided

5. **Abundance Calculation** (dbcan_asmfree)
   - Read-level abundance calculation
   - CAZyme family and subfamily abundance
   - Gene-level abundance estimation

6. **Visualization** (dbcan_plot)
   - Generation of abundance plots and heatmaps
   - CAZyme family distribution plots

7. **Report Generation** (MultiQC)
   - Aggregated quality control and analysis reports

Usage
-----

Basic Usage
~~~~~~~~~~~

The simplest command to run assembly-free analysis:

.. code-block:: bash

   nextflow run nf-core/dbcanmicrobiome \
     --input samplesheet.csv \
     --outdir results \
     --type assemfree \
     -profile docker

Skipping Steps
~~~~~~~~~~~~~~

You can skip certain steps if needed:

.. code-block:: bash

   nextflow run nf-core/dbcanmicrobiome \
     --input samplesheet.csv \
     --outdir results \
     --type assemfree \
     --skip_fastqc \
     --skip_trimming \
     --skip_kraken_extraction \
     -profile docker

With Custom Options
~~~~~~~~~~~~~~~~~~~

Combine assembly-free mode with other parameters:

.. code-block:: bash

   nextflow run nf-core/dbcanmicrobiome \
     --input samplesheet.csv \
     --outdir results \
     --type assemfree \
     --skip_kraken_extraction \
     -profile docker

Key Differences from Assembly-Based Modes
-----------------------------------------

1. **No Assembly**: Reads are directly annotated without assembly
2. **DIAMOND Database**: Uses CAZy DIAMOND database instead of full dbCAN database
3. **Read-Level Analysis**: Abundance is calculated at read level, not gene level
4. **Faster Processing**: Significantly faster than assembly-based workflows
5. **No Gene Prediction**: Does not perform gene prediction (Pyrodigal)
6. **No CGC Detection**: Does not detect CAZyme Gene Clusters (CGCs)

Output Files
------------

The pipeline generates output files organized in the following directory structure:

Read Conversion
~~~~~~~~~~~~~~~

- ``seqtk/``
  - ``*.fasta``: Converted FASTA files from FASTQ reads

CAZyme Search
~~~~~~~~~~~~~

- ``diamond/``
  - ``*.out``: DIAMOND blastx output files
  - ``*.out.tc``: DIAMOND output with taxonomy classification

Abundance Calculation
~~~~~~~~~~~~~~~~~~~~~

- ``dbcan_asmfree_cal_abund/``
  - ``*_abund/``: Abundance calculation results
  
    Files in this directory:
    
    - ``*_abund.txt``: CAZyme abundance values
    - ``*_family_abund.txt``: Family-level abundance
    - ``*_subfamily_abund.txt``: Subfamily-level abundance

Visualization
~~~~~~~~~~~~~

- ``dbcan_plot/``
  - ``*_pdf/``: Visualization plots
  
    Files in this directory:
    
    - ``heatmap.pdf``: Abundance heatmap
    - ``family.pdf``: CAZyme family distribution
    - ``subfamily.pdf``: Subfamily distribution

Quality Control
~~~~~~~~~~~~~~~

- ``fastqc/``: FastQC reports (if not skipped)
- ``trimgalore/``: TrimGalore reports (if not skipped)
- ``multiqc/``: MultiQC aggregated report

Pipeline Information
~~~~~~~~~~~~~~~~~~~~

- ``pipeline_info/``: Execution reports, parameters, and software versions

Advantages
----------

- **Computational Efficiency**: No assembly step reduces memory and CPU requirements
- **Speed**: Faster processing compared to assembly-based workflows
- **Scalability**: Better suited for very large datasets
- **Simplicity**: Simpler workflow with fewer steps
- **No Assembly Bias**: Avoids potential artifacts from assembly process

Limitations
-----------

- **No Gene-Level Analysis**: Cannot identify specific genes or gene clusters
- **No CGC Detection**: Does not detect CAZyme Gene Clusters
- **Read-Level Resolution**: Analysis is at read level, not contig level
- **Less Complete Information**: May miss some CAZymes that require assembly for detection
- **No RNA-seq Support**: Does not support RNA-seq data integration

When to Use Assembly-Free Mode
--------------------------------

**Recommended**:
- Very large datasets (>100M reads per sample)
- Quick CAZyme screening or preliminary analysis
- When computational resources are limited
- When assembly is not feasible or desired
- When read-level abundance is sufficient

**Not Recommended**:
- When gene-level or CGC-level analysis is required
- When complete gene sequences are needed
- When RNA-seq integration is required
- When detailed functional analysis is needed
- When assembly quality is important

Best Practices
--------------

1. **Quality Control**: Ensure good read quality before annotation
2. **Read Length**: Longer reads generally provide better DIAMOND matches
3. **Coverage**: Ensure sufficient sequencing depth for reliable abundance estimates
4. **Database**: The CAZy database is automatically downloaded; ensure sufficient disk space
5. **Validation**: Compare results with assembly-based mode when possible

Comparison with Other Modes
----------------------------

.. list-table:: Mode Comparison
   :widths: 25 25 25 25
   :header-rows: 1

   * - Feature
     - Assembly-Free
     - Short Reads
     - Long Reads
   * - Assembly
     - No
     - Yes (MEGAHIT)
     - Yes (Flye)
   * - Gene Prediction
     - No
     - Yes (Pyrodigal)
     - Yes (Pyrodigal)
   * - CGC Detection
     - No
     - Yes
     - Yes
   * - Speed
     - Fastest
     - Medium
     - Slowest
   * - Memory Usage
     - Lowest
     - Medium
     - Highest
   * - Resolution
     - Read-level
     - Gene-level
     - Gene-level
   * - RNA-seq Support
     - No
     - Yes
     - Yes

See Also
--------

- :ref:`shortreads-mode` - Short reads assembly-based mode
- :ref:`longreads-mode` - Long reads assembly-based mode
- :ref:`nextflow-parameters` - Complete parameter reference
- `DIAMOND Documentation <https://github.com/bbuchfink/diamond>`_ - DIAMOND BLASTX documentation
