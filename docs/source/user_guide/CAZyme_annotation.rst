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
     - Path to input sequence file (FASTA, optionally ``.gz``)
   * - ``--output_dir``
     - Directory for output files (created if needed)
   * - ``--db_dir``
     - Directory containing database files
   * - ``--mode``
     - Analysis mode: ``prok`` (prokaryotic genome, default), ``meta`` (metagenomic contigs), or ``protein`` (protein sequences)
   * - ``--methods``
     - Optional, comma-separated list of tools: ``diamond``, ``hmm``, ``dbCANsub``. Default: ``diamond,hmm,dbCANsub`` (all three).
       Example: ``--methods diamond,hmm``. Do not repeat the option; only the last ``--methods`` value is used.
   * - ``--threads``
     - Number of CPU threads for DIAMOND and pyHMMER. Default: all logical CPUs of the machine.
       (Gene prediction with pyrodigal currently always uses all CPUs; this is a known limitation.)
   * - ``--e_value_threshold``
     - Maximum E-value for DIAMOND hits against CAZy (default: ``1e-102``). The best hit per protein is kept.
   * - ``--e_value_threshold_dbcan`` / ``--coverage_threshold_dbcan``
     - Per-domain thresholds for dbCAN HMM hits (defaults: ``1e-15`` / ``0.35``).
   * - ``--e_value_threshold_dbsub`` / ``--coverage_threshold_dbsub``
     - Per-domain thresholds for dbCAN-sub HMM hits (defaults: ``1e-15`` / ``0.35``).

Run ``run_dbcan CAZyme_annotation --help`` for the full list, including memory-related options
(see :doc:`../FAQ/troubleshooting`).

.. note::

   **How the HMM thresholds are applied.** For every domain reported by pyHMMER, run_dbcan keeps the hit when its
   independent E-value (``i-Evalue``, per domain) is below the E-value threshold **and** its coverage is above the
   coverage threshold. Coverage is the fraction of the **HMM profile** covered by the alignment,
   ``(HMM To - HMM From + 1) / HMM Length``, not the fraction of the protein. Overlapping hits on the same protein
   (more than 50% of the shorter hit) are then reduced to the one with the better i-Evalue.

Usage Examples
---------------

Analyzing Prokaryotic Genomes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When working with bacterial or archaeal genomes, use the ``prok`` mode. Genes are predicted with pyrodigal and
written to ``uniInput.faa`` (proteins) and ``uniInput.gff`` (gene coordinates):

.. code-block:: shell

   run_dbcan CAZyme_annotation --input_raw_data EscheriaColiK12MG1655.fna --output_dir output_EscheriaColiK12MG1655_fna --db_dir db --mode prok

For metagenomic contigs or MAGs, use ``--mode meta`` (pyrodigal metagenomic mode, no training step).

Analyzing Protein Sequences
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For pre-translated protein sequences, use the ``protein`` mode:

.. code-block:: shell

   run_dbcan CAZyme_annotation --input_raw_data EscheriaColiK12MG1655.faa --output_dir output_EscheriaColiK12MG1655_faa --db_dir db --mode protein

In ``protein`` mode the FASTA IDs are copied into ``uniInput.faa`` with these changes:

* only the first word of each header is kept;
* IDs containing ``JGI`` (case-insensitive) have every ``|`` replaced by ``-``;
* IDs longer than 80 characters are truncated to 80 characters;
* duplicated IDs are renamed by appending ``_1``, ``_2``, ... to the second and later copies.

These changes are not applied to your GFF file, so for CGC prediction keep protein IDs unique and at most 80 characters.

.. warning::

   Use protein IDs that contain at least one letter. Purely numeric IDs (``>0``, ``>1``, ...) and IDs such as
   ``NA``, ``NULL``, ``nan`` or ``None`` are read as numbers or missing values when the result tables are merged. This is a
   known limitation: ``overview.tsv`` then lists the proteins without any annotation (``#ofTools`` = 0) or is not
   written at all. Prefix such IDs (for example ``>prot_0``) before running run_dbcan.

Analyzing Eukaryotic Proteins
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There is no eukaryotic genome mode: pyrodigal predicts prokaryotic genes only. For fungi and other eukaryotes,
predict genes with a eukaryotic gene finder (or use a published annotation) and run the predicted proteins in
``protein`` mode:

.. code-block:: shell

   run_dbcan CAZyme_annotation --input_raw_data Xylona_heveae_TC161.faa --output_dir output_Xylona_heveae_TC161_faa --db_dir db --mode protein

.. code-block:: shell

   run_dbcan CAZyme_annotation --input_raw_data Xylhe1_GeneCatalog_proteins_20130827.aa.fasta --output_dir output_Xylhe1_faa --db_dir db --mode protein

To predict CGCs for eukaryotes you also need a matching GFF file; see :doc:`CGC_information_generation`.

.. tip::

   For large datasets, set the number of CPU cores with ``--threads``.
   The default is all logical CPUs of your machine.

Output Files
--------------

The annotation process generates several key output files in your specified output directory:

* ``uniInput.faa`` - Unified protein input used by all tools (pyrodigal translations, or your cleaned protein FASTA)
* ``uniInput.gff`` - pyrodigal gene coordinates (``prok`` and ``meta`` modes only)
* ``overview.tsv`` - Summary of identified CAZymes (see below)
* ``dbCAN_hmm_results.tsv`` - dbCAN HMM hits after filtering
* ``diamond.out`` - DIAMOND search results against CAZy
* ``dbCANsub_hmm_results.tsv`` - dbCAN-sub HMM hits, including subfamily EC numbers and substrates
* ``dbCANsub_hmm_raw.tsv`` - unfiltered dbCAN-sub hits used to build the previous file

.. _overview-columns:

Columns of ``overview.tsv``
~~~~~~~~~~~~~~~~~~~~~~~~~~~

``overview.tsv`` has one row per protein that has at least one hit from any of the selected tools.

.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - Column
     - Content
   * - ``Gene ID``
     - Protein ID from ``uniInput.faa``.
   * - ``EC#``
     - EC numbers attached to the matched dbCAN-sub subfamilies, as ``EC:count`` tokens (``count`` = number of
       characterized proteins with that EC in the subfamily). Tokens of one subfamily are separated by ``;``,
       subfamilies by ``|``. ``-`` when there is no dbCAN-sub hit or no EC.
   * - ``dbCAN_hmm``
     - dbCAN HMM families with protein coordinates, e.g. ``GH23(101-244)``; several domains joined with ``+``.
   * - ``dbCAN_sub``
     - dbCAN-sub subfamilies with coordinates, e.g. ``GH23_e819(102-244)``; joined with ``+``.
   * - ``DIAMOND``
     - CAZy families of the best DIAMOND hit in CAZy, joined with ``+`` (e.g. ``CBM50+GH23``).
   * - ``#ofTools``
     - Number of tools (0-3) that reported at least one hit for the protein.
   * - ``Recommend Results``
     - Consensus call, only filled when ``#ofTools`` >= 2, otherwise ``-``. Domains are joined with ``|`` in
       protein order. If both dbCAN HMM and dbCAN-sub hit, overlapping domains are merged and one name is kept per
       domain, preferring a subfamily from dbCAN HMM (name containing ``_``, e.g. ``GH5_4``), then the dbCAN-sub
       subfamily, then the dbCAN HMM family. Otherwise the dbCAN HMM families are used, then dbCAN-sub, then DIAMOND.
   * - ``Substrate``
     - Substrates predicted from the dbCAN-sub hits: each subfamily's family and EC numbers are looked up in
       ``fam-substrate-mapping.tsv``. Multiple substrates are joined with ``;``. ``-`` when unknown.
       Note: in the current release, when a subfamily has no EC number (or its EC is not in the mapping file), the
       substrate falls back to the family-level entry of the mapping file, or to the family's only substrate if the
       family has exactly one. Such family-level substrates can differ from the dbCAN-sub web server, which reports
       none in this case (GitHub issue #85).
   * - ``SignalP`` / ``DeepTMHMM``
     - Only present with ``--run_signalp`` / ``--run_deeptmhmm``; see :doc:`signalp6_deeptmhmm`.

We recommend using proteins supported by at least two tools (``#ofTools`` >= 2). The CGC steps follow the same
rule: only proteins with ``#ofTools`` >= 2 are marked as CAZymes in ``cgc.gff``, and their ``Recommend Results``
value is used as the CAZyme annotation.

Customizing the Analysis
----------------------------

To customize which analytical methods are used:

.. code-block:: shell
   :caption: Using specific tools

   run_dbcan CAZyme_annotation --input_raw_data input.fna --output_dir output --db_dir db --mode prok --methods hmm,diamond

Available method combinations: ``hmm``, ``diamond``, ``dbCANsub``, or any comma-separated combination. Skipped tools
still get an empty result file, so ``overview.tsv`` is always generated. With a single method, ``#ofTools`` is at most 1
and ``Recommend Results`` stays ``-``, so no protein is used as a CAZyme in CGC prediction.

.. admonition:: Next Steps

   After completing CAZyme annotation, you may want to proceed to :doc:`CGC Information Generation <CGC_information_generation>` to identify CAZyme gene clusters.

.. tip::

   Optional **signal peptide** and **transmembrane topology** columns (SignalP 6.0 / DeepTMHMM) require separate installation and testing. See :doc:`signalp6_deeptmhmm`.
