.. _cazyme-annotation:

CAZyme Annotation
=================

Step 1 of the pipeline. Three tools search every protein: DIAMOND against CAZy, pyHMMER against the dbCAN family
HMMs and against the dbCAN-sub subfamily HMMs. ``overview.tsv`` merges them, one row per protein.

.. code-block:: shell

   run_dbcan CAZyme_annotation --mode prok    --input_raw_data genome.fna   --output_dir out --db_dir db
   run_dbcan CAZyme_annotation --mode meta    --input_raw_data contigs.fna  --output_dir out --db_dir db
   run_dbcan CAZyme_annotation --mode protein --input_raw_data proteins.faa --output_dir out --db_dir db

A protein FASTA is enough here; a GFF is needed only for the CGC steps.

Input modes
-----------

.. list-table::
   :widths: 14 86
   :header-rows: 1

   * - ``--mode``
     - Input
   * - ``prok``
     - Bacterial or archaeal genome (DNA). Genes are predicted with pyrodigal, trained on the input.
   * - ``meta``
     - Metagenome assembly, MAGs or many genomes in one file (DNA). pyrodigal metagenome mode, no training, streamed.
   * - ``protein``
     - Any proteome, including fungi and other eukaryotes annotated with
       `funannotate <https://github.com/nextgenusfs/funannotate>`_, `BRAKER <https://github.com/Gaius-Augustus/BRAKER>`_,
       NCBI or JGI.

Always give ``--mode`` (``--help`` marks it required). The input may be gzipped; any other value is rejected
(case-insensitive).

.. warning::

   pyrodigal predicts **prokaryotic** genes only. On a fungal genome it splits most genes at their introns
   (*Xylona heveae*: 32,471 ORFs of median 104 aa instead of 8,201 genes). run_dbcan logs a ``WARNING`` when the
   median protein is below 150 aa and more than 30% are shorter than 100 aa, but does not stop. Use
   ``--mode protein`` with a eukaryotic annotation instead.

Main options
------------

.. list-table::
   :widths: 40 60
   :header-rows: 1

   * - Option
     - Description (default)
   * - ``--methods``
     - Comma-separated tools: ``diamond``, ``hmm``, ``dbCANsub`` (all three). Give the option once.
   * - ``--threads``
     - Threads for pyrodigal, DIAMOND and pyHMMER (all logical CPUs).
   * - ``--e_value_threshold``
     - DIAMOND E-value vs CAZy; the best hit per protein is kept (``1e-102``).
   * - ``--e_value_threshold_dbcan`` / ``--coverage_threshold_dbcan``
     - dbCAN HMM per-domain i-Evalue / HMM coverage (``1e-15`` / ``0.35``).
   * - ``--e_value_threshold_dbsub`` / ``--coverage_threshold_dbsub``
     - Same for dbCAN-sub (``1e-15`` / ``0.35``).
   * - ``--run_signalp``, ``--run_deeptmhmm``
     - Optional signal peptide / topology columns: :doc:`signalp6_deeptmhmm`.

Memory options (``--large``, ``--batch_size``, ...) are in ``run_dbcan CAZyme_annotation --help`` and in
:doc:`../FAQ/troubleshooting`.

.. note::

   **HMM thresholds.** A domain is kept when its independent E-value (``i-Evalue``) is below the threshold and its
   coverage of the **HMM profile**, ``(HMM To - HMM From + 1) / HMM Length``, is above it (not the fraction of the
   protein). Hits that overlap by more than 50% of the shorter one are reduced to the one with the better i-Evalue.

Output files
------------

.. list-table::
   :widths: 32 68
   :header-rows: 1

   * - File
     - Content
   * - ``overview.tsv``
     - Merged result, one row per protein with at least one hit (columns below)
   * - ``uniInput.faa``
     - Proteins searched by all tools (pyrodigal translations, or your cleaned FASTA)
   * - ``uniInput.gff``
     - pyrodigal genes (``prok``/``meta`` only); used automatically by the CGC steps
   * - ``diamond.out``
     - DIAMOND hits against CAZy
   * - ``dbCAN_hmm_results.tsv``
     - dbCAN HMM domains after filtering
   * - ``dbCANsub_hmm_results.tsv``
     - dbCAN-sub domains with subfamily EC numbers and substrates (``dbCANsub_hmm_raw.tsv``: unfiltered)

.. _overview-columns:

Columns of ``overview.tsv``
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - Column
     - Content
   * - ``Gene ID``
     - Protein ID from ``uniInput.faa``.
   * - ``EC#``
     - EC numbers of the matched dbCAN-sub subfamilies as ``EC:count`` (count = characterized proteins with that EC
       in the subfamily; written ``:1`` when the eHMM name has no count). ``;`` between the ECs of one subfamily,
       ``|`` between subfamilies, ``-`` if none.
   * - ``dbCAN_hmm``
     - dbCAN families with protein coordinates, e.g. ``GH23(101-244)``; domains joined with ``+``.
   * - ``dbCAN_sub``
     - dbCAN-sub subfamilies with coordinates, e.g. ``GH23_e819(102-244)``; joined with ``+``.
   * - ``DIAMOND``
     - CAZy families of the best CAZy hit, joined with ``+`` (e.g. ``CBM50+GH23``). Only family and subfamily names
       are kept from the CAZy header (no accessions, EC numbers or file names), each once.
   * - ``#ofTools``
     - Number of tools (0-3) with a hit.
   * - ``Recommend Results``
     - Consensus, only when ``#ofTools`` >= 2 (else ``-``); one name per domain, ``|``-joined in protein order.
       Preference: a dbCAN HMM subfamily (name with ``_``, e.g. ``GH5_4``) > the dbCAN-sub subfamily > the dbCAN HMM
       family; without both HMM tools: dbCAN HMM, then dbCAN-sub, then DIAMOND.
   * - ``Substrate``
     - Substrates of the dbCAN-sub hits, unique and sorted, joined with ``;`` (``-`` if unknown); see below.
   * - ``SignalP`` / ``DeepTMHMM``
     - Only with ``--run_signalp`` / ``--run_deeptmhmm``.

**Use proteins with** ``#ofTools`` **>= 2.** The CGC steps do the same: only these proteins are CAZymes in
``cgc.gff``, named by their ``Recommend Results``.

**How the** ``Substrate`` **column is filled.** For every EC number of a matched subfamily:

#. ``dbCAN-sub-substrate.tsv`` in ``--db_dir`` (optional; the dbCAN-sub web server table) is used when it has the
   (subfamily, EC) pair, so the result matches the web server;
#. otherwise (family, EC) is looked up in ``fam-substrate-mapping.tsv``, first for the subfamily's own family, then
   for the other families of its composition (the EC of ``CBM1_e3`` belongs to GH45);
#. a complete EC without an entry, and a subfamily without EC, fall back to the family-level substrate (or to the
   family's only substrate); an incomplete EC (``a.b.c.-``) has no fallback.

.. note::

   Subfamilies without any EC therefore get the family-level substrate, e.g. GH71 -> ``alpha-glucan``, while the
   dbCAN-sub web server shows **no** substrate for them. This is the expected difference between the two.

In ``dbCANsub_hmm_results.tsv`` the ``Substrate`` column has one entry per EC, in the order of ``Subfam EC``,
separated by ``;``. Several substrates of one (family, EC) are joined with ``,`` inside that entry; ``-`` = no
substrate for that EC.

Protein input: IDs and sequences
--------------------------------

In ``--mode protein`` the FASTA is copied into ``uniInput.faa`` with these changes:

* only the first word of each header is kept; IDs containing ``JGI`` (any case) get ``-`` instead of ``|``;
* IDs longer than 80 characters are truncated; duplicated IDs get ``_1``, ``_2``, ... (keep IDs unique: these
  suffixes cannot be matched to a GFF);
* a trailing stop (``*`` or ``.``) is removed and internal ones become ``X`` (``gffread -y`` writes ``.``, which
  DIAMOND rejects; ``gffread -y ... -S`` writes ``*``).

The GFF is not rewritten; its IDs are cleaned the same way when they are matched for CGC prediction.

IDs are always read as text, so numeric (``>0``, ``>007``) and ``NA``/``None``-like IDs are kept (5.2.9 and earlier
lost them in ``overview.tsv``; there, prefix them, e.g. ``>prot_0``).

Choosing tools with ``--methods``
---------------------------------

.. code-block:: shell

   run_dbcan CAZyme_annotation --input_raw_data genome.fna --mode prok --output_dir out --db_dir db --methods hmm,diamond

* Tools that are not selected get an empty result file, so ``overview.tsv`` is always written.
* With **one** tool ``#ofTools`` is at most 1: ``Recommend Results`` stays ``-`` and no protein is a CAZyme for
  CGC prediction.
* Only the selected tools are merged. Result files of other tools left in the folder by an earlier run are ignored
  (with a warning), and the files of the selected tools are deleted at the start, so a tool that fails contributes
  an empty table and the command exits with a non-zero status.

Next: :doc:`CGC_information_generation` (step 2).
