.. _cazyme-annotation:

CAZyme Annotation
==================

Introduction
-------------

CAZyme annotation is a critical step in identifying and classifying Carbohydrate-Active Enzymes (CAZymes) in biological sequences. The ``run_dbcan`` tool enables comprehensive annotation of CAZymes from various input types:

* Prokaryotic genomes (nucleotide sequences)
* Metagenomic contigs (nucleotide sequences)
* Protein sequences (prokaryotic or eukaryotic)

The annotation process integrates multiple analytical tools to ensure high sensitivity and specificity in CAZyme identification.

Command Syntax
----------------

.. code-block:: shell

   run_dbcan CAZyme_annotation --input_raw_data <INPUT_FILE> --output_dir <OUTPUT_DIRECTORY> --db_dir <DATABASE_DIRECTORY> --mode <MODE>

Key Parameters
~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - Parameter
     - Description
   * - ``--input_raw_data``
     - Path to input sequence file (FASTA format)
   * - ``--output_dir``
     - Directory for output files
   * - ``--db_dir``
     - Directory containing database files
   * - ``--mode``
     - Analysis mode: ``prok`` (prokaryote), ``meta`` (metagenome), or ``protein`` (protein sequences)
   * - ``--methods``
     - Optional: Specify tools to use (``diamond``, ``hmm``, and/or ``dbCANsub``, default is ``all``)
       Usage: ``--methods diamond --methods hmm --methods dbCANsub`` for multiple methods or just choose one/two.

Usage Examples
---------------

Analyzing Prokaryotic Genomes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When working with bacterial or archaeal genomes, use the ``prok`` mode:

.. code-block:: shell

   run_dbcan CAZyme_annotation --input_raw_data EscheriaColiK12MG1655.fna --output_dir output_EscheriaColiK12MG1655_fna --db_dir db --mode prok

Analyzing Protein Sequences
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For pre-translated protein sequences, use the ``protein`` mode:

.. code-block:: shell

   run_dbcan CAZyme_annotation --input_raw_data EscheriaColiK12MG1655.faa --output_dir output_EscheriaColiK12MG1655_faa --db_dir db --mode protein

Analyzing Eukaryotic Proteins
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Eukaryotic proteins are processed the same way using ``protein`` mode:

.. code-block:: shell

   run_dbcan CAZyme_annotation --input_raw_data Xylona_heveae_TC161.faa --output_dir output_Xylona_heveae_TC161_faa --db_dir db --mode protein

.. code-block:: shell

   run_dbcan CAZyme_annotation --input_raw_data Xylhe1_GeneCatalog_proteins_20130827.aa.fasta --output_dir output_Xylhe1_faa --db_dir db --mode protein

.. tip::

   For large eukaryotic datasets, consider change the computational resources with ``--threads`` to specify the number of CPU cores.
   The default is ``all cores`` of your machine.

Output Files
--------------

The annotation process generates several key output files in your specified output directory:

* ``uniInput.faa`` - Unified input file for all tools
* ``overview.txt`` - Summary of identified CAZymes
* ``dbCAN_hmm_results.tsv`` - Detailed HMMER results
* ``diamond.out`` - DIAMOND search results
* ``dbCANsub_hmm_results.tsv`` - dbCAN sub-HMM results including substrate specificity

Customizing the Analysis
----------------------------

To customize which analytical methods are used:

.. code-block:: shell
   :caption: Using specific tools

   run_dbcan CAZyme_annotation --input_raw_data input.fna --output_dir output --db_dir db --mode prok --methods hmm --methods diamond

Available method combinations: ``hmm``, ``diamond``, ``dbCANsub``, or any combination.

.. admonition:: Next Steps

   After completing CAZyme annotation, you may want to proceed to :doc:`CGC Information Generation <CGC_information_generation>` to identify CAZyme gene clusters.

