Nextflow Pipeline Documentation
=================================

The Nextflow pipeline for CAZyme annotation in microbiome data provides an automated, scalable workflow for processing metagenomic sequencing data. This documentation covers how to use the pipeline and understand its outputs.

The pipeline is built using `Nextflow <https://www.nextflow.io/>`_ and follows best practices for reproducible bioinformatics workflows. It supports multiple analysis modes including short-read assembly, long-read assembly, and assembly-free approaches for CAZyme annotation in microbiome datasets.

Quick Start
-----------

The pipeline supports three main analysis modes:

- **Short Reads** (``--type shortreads``): Assembly-based analysis for Illumina short-read data using MEGAHIT. Supports subsampling and co-assembly options.
- **Long Reads** (``--type longreads``): Assembly-based analysis for PacBio/Nanopore long-read data using Flye.
- **Assembly Free** (``--type assemfree``): Direct annotation without assembly using DIAMOND blastx, ideal for large datasets.

Installation and basic usage:

.. code-block:: bash

   # Clone the repository
   git clone https://github.com/bcb-unl/dbcan-nf.git
   cd dbcan-nf

   # Run the pipeline
   nextflow run main.nf \
     --input samplesheet.csv \
     --outdir results \
     --type shortreads \
     -profile docker \
     --skip_kraken_extraction # based on the database size of kraken2, you can skip this step if the database is too large.

For detailed information about each mode, see the corresponding documentation pages below. Example results and visualizations are available in the :ref:`results examples <nextflow-results-examples>` section.

Documentation Sections
----------------------

Getting Started
~~~~~~~~~~~~~~~

.. toctree::
   :maxdepth: 1

   usage

Analysis Modes
~~~~~~~~~~~~~~

.. toctree::
   :maxdepth: 2

   shortreads
   shortreads_subsample
   shortreads_coassembly
   longreads
   assemfree

Reference
~~~~~~~~~

.. toctree::
   :maxdepth: 1

   parameters
   output
   results_examples
   hpc

