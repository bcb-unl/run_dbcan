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
     - Input type (case-insensitive; any other value is rejected): ``prok`` (bacterial/archaeal genome, DNA;
       default), ``meta`` (metagenome assembly or MAGs, DNA) or ``protein`` (any proteome, including fungi and other
       eukaryotes annotated with funannotate, BRAKER, NCBI or JGI). ``prok`` and ``meta`` predict genes with
       pyrodigal, which is prokaryote-only.
   * - ``--methods``
     - Optional, comma-separated list of tools: ``diamond``, ``hmm``, ``dbCANsub``. Default: ``diamond,hmm,dbCANsub`` (all three).
       Example: ``--methods diamond,hmm``. Do not repeat the option; only the last ``--methods`` value is used.
   * - ``--threads``
     - Number of CPU threads for gene prediction (pyrodigal, ``prok``/``meta`` modes), DIAMOND and pyHMMER.
       Default: all logical CPUs of the machine.
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
   ``(HMM To - HMM From + 1) / HMM Length``, not the fraction of the protein. Hits on the same protein that share
   more than 50% of the shorter hit's match-state residues are then reduced to the one with the better i-Evalue
   (see `Match-state overlap resolution`_).

Usage Examples
---------------

Analyzing Prokaryotic Genomes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When working with bacterial or archaeal genomes, use the ``prok`` mode. Genes are predicted with pyrodigal and
written to ``uniInput.faa`` (proteins) and ``uniInput.gff`` (gene coordinates):

.. code-block:: shell

   run_dbcan CAZyme_annotation --input_raw_data EscheriaColiK12MG1655.fna --output_dir output_EscheriaColiK12MG1655_fna --db_dir db --mode prok

For metagenomic contigs or MAGs, use ``--mode meta`` (pyrodigal metagenomic mode, no training step).

.. warning::

   pyrodigal models prokaryotic genes (no introns). On a fungal or other eukaryotic genome it splits most genes
   into short ORFs: on *Xylona heveae* (24 Mb) it predicts 32,471 ORFs instead of 8,201 genes, with a median length
   of 104 aa, and most CAZymes are split across several ORFs. After gene prediction run_dbcan logs a ``WARNING`` when
   the proteins look like this (median below 150 aa and more than 30% shorter than 100 aa); the run is not
   stopped. Use ``--mode protein`` with an external annotation instead (see `Analyzing Eukaryotic Proteins`_).

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

The protein sequences are also cleaned for DIAMOND: a trailing stop symbol (``*`` or ``.``) is removed and internal
``*`` or ``.`` are replaced by ``X``. ``gffread -y`` writes ``.`` for stop codons, which DIAMOND rejects. The number
of changed sequences is logged.

The GFF file is not rewritten. For CGC prediction the identifiers of the GFF are passed through the same ID
cleaning before they are compared with ``uniInput.faa``, so a JGI or long ID still matches; the ``_1``, ``_2``
suffixes of duplicated IDs cannot be matched, so keep protein IDs unique.

.. warning::

   Use protein IDs that contain at least one letter. Purely numeric IDs (``>0``, ``>1``, ...) and IDs such as
   ``NA``, ``NULL``, ``nan`` or ``None`` are read as numbers or missing values when the result tables are merged. This is a
   known limitation: ``overview.tsv`` then lists the proteins without any annotation (``#ofTools`` = 0) or is not
   written at all. Prefix such IDs (for example ``>prot_0``) before running run_dbcan.

Analyzing Eukaryotic Proteins
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There is no eukaryotic genome mode, and run_dbcan does not predict eukaryotic genes. For fungi and other
eukaryotes, annotate the genome with a eukaryotic pipeline such as `funannotate <https://github.com/nextgenusfs/funannotate>`_
/ funannotate2 or `BRAKER <https://github.com/Gaius-Augustus/BRAKER>`_ (or use the NCBI or JGI MycoCosm
annotation), and run its protein FASTA in ``protein`` mode:

.. code-block:: shell

   run_dbcan CAZyme_annotation --input_raw_data Xylona_heveae_TC161.faa --output_dir output_Xylona_heveae_TC161_faa --db_dir db --mode protein

.. code-block:: shell

   run_dbcan CAZyme_annotation --input_raw_data Xylhe1_GeneCatalog_proteins_20130827.aa.fasta --output_dir output_Xylhe1_faa --db_dir db --mode protein

To predict CGCs for eukaryotes, also pass the GFF3/GTF of the same annotation run with ``--input_gff`` and
``--gff_type`` (``funannotate``, ``BRAKER``, ``GFF3``, ``GTF``, ``NCBI_euk`` or ``JGI``); see
:ref:`cgc-information-generation` for the fungal workflow.

.. tip::

   If you translate a GFF yourself, use ``gffread annotation.gff3 -g genome.fna -y proteins.faa -S``: ``-S`` writes
   ``*`` instead of ``.`` for stop codons. run_dbcan cleans both, but other tools may not.

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
* ``dbCAN_hmm_alignment_details.tsv`` / ``dbCANsub_hmm_alignment_details.tsv`` - per-domain alignment diagnostics
  (alignment span, match-state intervals, i-Evalue, coverage); written while match-state overlap resolution is
  active (the default; also for the TF, STP and Pfam HMM tables of the CGC steps), see
  `Match-state overlap resolution`_ below

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
       characterized proteins with that EC in the subfamily; the eHMM name omits a count of 1, which is written as ``:1``). Tokens of one subfamily are separated by ``;``,
       subfamilies by ``|``. ``-`` when there is no dbCAN-sub hit or no EC.
   * - ``dbCAN_hmm``
     - dbCAN HMM families with protein coordinates, e.g. ``GH23(101-244)``; several domains joined with ``+``.
   * - ``dbCAN_sub``
     - dbCAN-sub subfamilies with coordinates, e.g. ``GH23_e819(102-244)``; joined with ``+``.
   * - ``DIAMOND``
     - CAZy families of the best DIAMOND hit in CAZy, joined with ``+`` (e.g. ``CBM50+GH23``). Only whole family
       or subfamily names (``GH13``, ``GH13_6``, ``CBM50``) are taken from the CAZy sequence header; accessions
       (also those starting with ``AA``/``CE``/``GT``, e.g. ``AAK49401.1``), EC numbers, JGI protein numbers and
       file names are dropped, and a family repeated in the header is listed once.
   * - ``#ofTools``
     - Number of tools (0-3) that reported at least one hit for the protein.
   * - ``Recommend Results``
     - Consensus call, only filled when ``#ofTools`` >= 2, otherwise ``-``. Domains are joined with ``|`` in
       protein order. If both dbCAN HMM and dbCAN-sub hit, overlapping domains are merged and one name is kept per
       domain, preferring a subfamily from dbCAN HMM (name containing ``_``, e.g. ``GH5_4``), then the dbCAN-sub
       subfamily, then the dbCAN HMM family. Otherwise the dbCAN HMM families are used, then dbCAN-sub, then DIAMOND.
   * - ``Substrate``
     - Substrates predicted from the dbCAN-sub hits. For every EC number of a matched subfamily:

       1. if the database directory contains ``dbCAN-sub-substrate.tsv`` (the per-subfamily substrate table of the
          dbCAN-sub web server), its entry for (subfamily, EC) is used, so the result matches the web server;
       2. otherwise (family, EC) is looked up in ``fam-substrate-mapping.tsv`` for the subfamily's own family and
          then for the other families of its composition (the EC of a CBM or multi-domain subfamily usually belongs
          to its catalytic family, e.g. ``CBM1_e3`` with GH45 and 3.2.1.4);
       3. a complete EC number without an entry falls back to the family-level entry of the mapping file, or to
          the family's only substrate when the family has exactly one; an incomplete EC (``a.b.c.-``) has no
          fallback.

       Subfamilies without any EC number get the family-level substrate (step 3), for example GH71 ->
       ``alpha-glucan``; the web server shows no substrate for them. Unique substrates of all hits are joined with
       ``;``; ``-`` when unknown.
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

``overview.tsv`` only merges the results of the methods selected in the current run. When the output directory
still holds result files of other methods from an earlier run (for example ``diamond.out`` after re-running with
``--methods hmm``), those files are left in place but ignored, and a warning says so (``Ignoring diamond.out ...``).
The result files of the selected methods are deleted when the run starts, so a selected method that fails in this
run contributes nothing (an empty table) to ``overview.tsv`` instead of its result from an earlier run; the command
then exits with a non-zero status.

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

``--dbcan-taxon-profile fungi`` loads the bundled fungal policy
(``dbcan/data/taxon_profiles/fungi.cutoff.json``) instead of a policy file; it
is experimental, locked to the ``dbCAN.hmm`` release it was calibrated on (the
run stops on another release), and explicit ``--dbcan-*-file`` options override
it. In ``easy_CGC`` / ``easy_substrate`` it also turns on ``--fungi`` (TF HMM
instead of TF DIAMOND) unless ``--fungi`` or ``--no-fungi`` is given.

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
