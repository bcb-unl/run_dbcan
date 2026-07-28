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
* ``dbCAN_hmm_alignment_details.tsv`` - per-domain alignment diagnostics
  (alignment span, match-state intervals, i-Evalue, coverage); written for each
  HMM library while match-state overlap resolution is active, see below

Customizing the Analysis
----------------------------

To customize which analytical methods are used:

.. code-block:: shell
   :caption: Using specific tools

   run_dbcan CAZyme_annotation --input_raw_data input.fna --output_dir output --db_dir db --mode prok --methods hmm --methods diamond

Available method combinations: ``hmm``, ``diamond``, ``dbCANsub``, or any combination.

Experimental dbCAN HMM policies
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Family-specific bit-score/coverage cutoffs and pair-specific overlap
coexistence are opt-in research interfaces. They do not change the default
pipeline. A policy is validated against the exact ``dbCAN.hmm`` SHA-256,
profile count, schema, and HMM/family names before the search starts.

.. code-block:: shell

   run_dbcan CAZyme_annotation \
     --input_raw_data proteins.faa \
     --output_dir output \
     --db_dir db \
     --mode protein \
     --methods hmm \
     --dbcan-family-cutoff-file family_cutoffs.json \
     --dbcan-overlap-coexistence-file overlap_pairs.json \
     --allow-experimental-dbcan-cutoffs

An experimental policy cannot run without
``--allow-experimental-dbcan-cutoffs``. A family cutoff is applied before the
normal overlap filter; families absent from the policy retain the legacy global
E-value and coverage rule. Policy files calibrated for one taxon or one HMM
release must not be reused automatically for another.

Match-state overlap resolution
----------------------------------

Domain conflicts are resolved using target residues assigned to profile-HMM
match states rather than the continuous sequence alignment span. Residues
assigned to HMM insert states do not contribute to match-state overlap.

HMMER reports each domain as a continuous ``Target From``-``Target To`` span,
but a profile alignment may contain match (M), insert (I) and delete (D) states.
A target residue aligned to an insert state lies between two consensus match
positions and still falls inside the continuous span, so a large catalytic
profile can geometrically span an inserted accessory domain without having any
match-state evidence over it:

.. code-block:: text

   GH18 alignment span     1 ------------------------------------- 400
   GH18 match states       1 -------------- 250      340 --------- 400
   GH18 insert states                    251 ------------ 339
   CBM18                                        300 --- 334

Under the historical rule the CBM18 span lies almost entirely inside the GH18
span, so the CBM was removed. Under match-state semantics the two share no
match-state residues, there is no conflict, and both domains are kept.

.. note::

   An insert state is an HMM alignment concept, not a biological one. A
   lineage-specific loop or an inserted structural region may still belong to
   the domain biologically. Match-state occupancy is used here because it is
   the evidence the profile actually provides for that stretch of sequence, not
   because inserted residues are "outside" the domain.

Geometric overlap alone does not suppress a hit. The decision is:

* no substantial match-state overlap → keep both;
* substantial match-state overlap covered by an explicit coexistence policy
  (``--dbcan-overlap-coexistence-file``) → keep both;
* substantial match-state overlap between competing models → resolve by
  evidence (best i-Evalue first).

"Substantial" means strictly more than 50% of the shorter domain's match-state
length; an overlap of exactly 50% is not a conflict. Resolution runs over a
per-protein conflict graph, so the outcome does not depend on the order hits
appear in the input.

This behaviour is on by default for every HMM library in a run (dbCAN,
dbCAN-sub, TF, STP, Pfam). ``--dbcan-legacy-overlap`` restores the previous
alignment-span *geometry* and exists for reproducibility only. Note that the
order-independent conflict-graph resolution applies in both modes: the old
implementation compared each hit only against the most recently retained one,
which could keep or drop a hit depending on row order, and that is a bug rather
than behaviour worth reproducing.

The public result tables keep their historical columns in both modes:
``Target From`` and ``Target To`` remain the reported alignment location and are
not replaced by match-state intervals. While match-state resolution is active,
the intervals used for the decision are additionally written to a companion
``*_alignment_details.tsv`` file (for example
``dbCAN_hmm_alignment_details.tsv``) with the HMM name, target, alignment span,
match-state intervals, i-Evalue and coverage. Downstream steps read that file so
they consume already-resolved domains instead of re-deriving the architecture
from the continuous span.

.. note::

   HMM coverage is unchanged: it is still ``(hmm_to - hmm_from + 1) /
   hmm_length``, i.e. HMM coordinate-span coverage. That is not identical to
   match-state occupancy when delete states occur, but the ``0.35`` default is
   tied to historical calibration, so coverage semantics are deliberately left
   as they are.

.. admonition:: Next Steps

   After completing CAZyme annotation, you may want to proceed to :doc:`CGC Information Generation <CGC_information_generation>` to identify CAZyme gene clusters.

.. tip::

   Optional **signal peptide** and **transmembrane topology** columns (SignalP 6.0 / DeepTMHMM) require separate installation and testing. See :doc:`signalp6_deeptmhmm`.
