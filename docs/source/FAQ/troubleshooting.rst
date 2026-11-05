.. _troubleshooting:

Known Limitations and Troubleshooting
=====================================

This page collects questions that come up frequently on the issue tracker, together with known limitations of the
current release. Run any command with ``-v`` (or ``--log-level INFO``) and ``--log-file run.log`` to get detailed logs
when reporting a problem.

Input data
----------

**Can run_dbcan annotate a fungal or other eukaryotic genome?**
   Not from DNA. There is no eukaryotic genome mode: ``--mode prok`` and ``--mode meta`` use pyrodigal, which predicts
   prokaryotic genes only. Predict genes with a eukaryotic gene finder (or use an NCBI/JGI annotation), run the
   proteins with ``--mode protein``, and for CGCs provide the GFF with ``--input_gff`` and ``--gff_type NCBI_euk``
   (or ``JGI``). GFF3 files from funannotate or BRAKER need a ``protein_id`` attribute on their ``CDS`` rows first;
   see :ref:`cgc-information-generation`.

**Which protein IDs are safe?**
   Use unique IDs of at most 80 characters that contain at least one letter. In ``--mode protein``:

   * only the first word of each FASTA header is kept;
   * IDs longer than 80 characters are truncated;
   * duplicated IDs get ``_1``, ``_2``, ... appended;
   * IDs containing ``JGI`` have ``|`` replaced by ``-``.

   These changes are applied to ``uniInput.faa`` only, not to your GFF, so changed IDs no longer match the GFF.

**overview.tsv is empty or has no annotations, although the HMM/DIAMOND result files contain hits.**
   Known limitation: purely numeric protein IDs (``>0``, ``>1``, ...) are read as numbers when the result tables are
   merged, so the rows of ``overview.tsv`` show ``-`` everywhere and ``#ofTools`` = 0. IDs such as ``NA``, ``NULL``, ``nan``
   or ``None`` are read as missing values, and ``overview.tsv`` may not be written at all (the error is only logged).
   Rename such sequences with a prefix, for example:

   .. code-block:: shell

      awk '/^>/ { sub(/^>/, ">prot_") } { print }' proteins.faa > proteins.prefixed.faa

**Only one or two tools found a CAZyme. Which result should I use?**
   ``Recommend Results`` is only filled when at least two tools report a hit for the protein (``#ofTools`` >= 2; the tools do not have to report the same family), and only these proteins
   are used as CAZymes for CGC prediction. With ``--methods`` set to a single tool, no protein reaches this threshold.
   See :ref:`overview-columns` for how ``Recommend Results`` is chosen.

Thresholds
----------

**What does the HMM coverage threshold mean?**
   ``--coverage_threshold_dbcan`` / ``--coverage_threshold_dbsub`` (and the TF, STP and Pfam equivalents) are compared
   with the fraction of the **HMM profile** covered by a domain, ``(HMM To - HMM From + 1) / HMM Length``, not with
   the fraction of the protein. (The ``--help`` text of ``--coverage_threshold_dbcan`` says "fraction of query length";
   the code uses the HMM length.) The E-value thresholds are compared with the per-domain independent E-value
   (``i-Evalue``). A domain is kept when ``i-Evalue < threshold`` and ``coverage > threshold``.

**DIAMOND coverage options use percent, HMM options use fractions.**
   ``--coverage_threshold_tc``, ``--coverage_threshold_tf_diamond``, ``--coverage_threshold_sulfatase`` and
   ``--coverage_threshold_peptidase`` are passed to DIAMOND ``--query-cover`` and use 0-100 (default 35).
   The pyHMMER coverage options use 0-1 (default 0.35).

**Which unit does substrate_prediction -iden use?**
   Percent identity on a **0-100 scale** (DIAMOND ``pident``), although ``--help`` says 0.0-1.0. ``-iden 0.3`` keeps
   almost everything; use ``-iden 30`` for 30% identity. ``-cov`` is a fraction (0-1) of the CGC protein length.

CGC prediction
--------------

**No CGC is found (cgc_standard_out.tsv is empty).**
   Check, in this order:

   1. ``cgc.gff`` exists. ``gff_process`` only logs an error when the input GFF is missing. It has no ``--mode``
      option and defaults to ``<output_dir>/uniInput.gff``, which exists only after a ``prok``/``meta`` run; with
      protein input always pass ``--input_gff`` and ``--gff_type``.
   2. Genes are matched to proteins. ``grep -c "protein_id=unknown" cgc.gff`` should be small and
      ``grep -c "CGC_annotation=CAZyme" cgc.gff`` should be greater than zero. Unmatched genes are written silently as
      ``protein_id=unknown``/``CGC_annotation=null``. The attribute used per ``--gff_type`` is listed in
      :ref:`cgc-information-generation`.
   3. Enough CAZymes are supported by two tools (``#ofTools`` >= 2 in ``overview.tsv``).
   4. The default rule requires at least one CAZyme **and** one transporter (``--additional_genes TC``) within
      ``--num_null_gene 2`` genes of each other. Relax it, for example with
      ``--additional_genes TC --additional_genes TF --additional_genes STP --additional_logic any``.

**--feature_type has no effect.**
   Known limitation: ``cgc_finder`` currently uses all lines of ``cgc.gff`` regardless of ``--feature_type``.

Substrate prediction
--------------------

**substrate_prediction asks for --mode and --input_raw_data.**
   Fixed after 5.2.9 (GitHub issue #81): ``substrate_prediction`` only needs ``--output_dir`` and ``--db_dir``.
   With 5.2.9 and earlier, pass the same ``--mode``/``--input_raw_data`` values as in the annotation step.

   .. code-block:: shell

      run_dbcan substrate_prediction --output_dir output_dir --db_dir db

   See :doc:`../user_guide/predict_CGC_substrate` for the files that must exist in ``--output_dir``.

Databases
---------

**HTTP or --aws_s3?**
   The default HTTP source is the moving ``db_current`` snapshot; ``--aws_s3`` downloads the release pinned in the
   installed package (currently ``db_v5-2-9_11-3-2026``). Use ``--aws_s3`` for reproducible results and record the
   release name, the ``run_dbcan version`` output and the file checksums (see :doc:`../user_guide/prepare_the_database`).

**How much disk space do I need?**
   About 8.5 GB for all files (``dbCAN-sub.hmm`` ~5.1 GB, ``CAZy.dmnd`` ~2.6 GB, ``dbCAN.hmm`` ~120 MB, CGC databases
   ~0.5 GB). ``Pfam-A.hmm`` (~1.6 GB) is only needed for ``Pfam_null_cgc`` and must be downloaded manually.

**A download was interrupted.**
   Rerun the same command. Partial files (``*.part``) are resumed, complete files are skipped when their size matches,
   and ``--no-overwrite`` skips every file that already exists.

Memory and run time
-------------------

**How many threads are used?**
   ``--threads`` defaults to the number of logical CPUs of the machine and is used for DIAMOND and pyHMMER. On shared
   servers or cluster jobs, set it to the number of cores you were allocated. Gene prediction with pyrodigal
   currently uses all CPUs regardless of ``--threads`` (known limitation). ``substrate_prediction`` run on its own
   always uses all CPUs for DIAMOND.

**The HMM search runs out of memory on a very large input.**
   * ``dbCAN-sub.hmm`` is 5.1 GB. HMM profiles are loaded into memory only when about five times the HMM file size
     fits into 70% of the available RAM (roughly 37 GB free for dbCAN-sub); otherwise they are streamed from disk,
     which needs less memory but is slower.
   * Inputs larger than ``--large_input_threshold_mb`` (default 5000 MB) switch to a streaming-safe "large" mode
     automatically. Force it with ``--large`` (dbCAN HMM) and ``--large_dbsub`` (dbCAN-sub).
   * Set a fixed batch size with ``--batch_size`` / ``--batch_size_dbsub``, or lower
     ``--memory_safety_factor`` / ``--memory_safety_factor_dbsub`` for smaller automatic batches. After a
     ``MemoryError`` the search is retried up to ``--max_retries`` times with smaller batches.
   * Fewer ``--threads`` also reduce peak memory.

**Can I split a large protein file into chunks?**
   Yes, for CAZyme annotation. The HMM searches use a fixed ``Z`` (the number of profiles in the database), so E-values
   do not depend on how many proteins are searched together, and ``overview.tsv`` has one independent row per
   protein. Run each chunk into its own ``--output_dir`` and concatenate the ``overview.tsv`` files (keep one header).
   Do not split genomes for CGC prediction, because CGCs need all genes of a contig.

**Prokaryotic genome mode is slow or uses much memory on a metagenome.**
   ``--mode prok`` reads all sequences into memory to train pyrodigal on them. For metagenomic assemblies or many
   genomes concatenated together, use ``--mode meta``, which processes sequences in a stream without training.

Optional tools
--------------

**SignalP or DeepTMHMM columns are missing.**
   Both tools must be installed separately. DeepTMHMM currently runs only when ``--run_signalp`` is also given, and
   failures are logged without stopping ``CAZyme_annotation``. See :doc:`../user_guide/signalp6_deeptmhmm`.
