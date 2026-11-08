.. _troubleshooting:

Known Limitations and Troubleshooting
=====================================

Frequent questions from the issue tracker and known limitations of the current release. When you report a
problem, add ``-v`` (or ``--log-level INFO``) and ``--log-file run.log`` to the command and attach the log.

Input data
----------

**Can run_dbcan annotate a fungal or other eukaryotic genome?**
   Not from DNA. ``--mode prok`` and ``--mode meta`` use pyrodigal, which predicts prokaryotic genes only (bacteria,
   archaea, metagenomes); on a fungal genome it splits most genes at their introns, and run_dbcan logs a warning
   when the predicted proteins look like that. run_dbcan does not predict eukaryotic genes. Annotate the genome
   with funannotate/funannotate2 or BRAKER (or use the NCBI/JGI annotation), run its proteins with
   ``--mode protein``, and for CGCs pass its GFF3/GTF with ``--input_gff`` and ``--gff_type funannotate``,
   ``BRAKER``, ``GFF3``, ``GTF``, ``NCBI_euk`` or ``JGI``; see :ref:`cgc-information-generation`.

**Do I need a GFF file for protein input?**
   Only for CGCs. ``CAZyme_annotation --mode protein`` needs only the protein FASTA. For CGC prediction
   (``gff_process``, ``easy_CGC``, ``easy_substrate``) from proteins, give the GFF of the same annotation with
   ``--input_gff`` and its ``--gff_type``; DNA input (``prok``/``meta``) uses the pyrodigal GFF automatically.

**cgc_standard_out.tsv is empty / gff_process stops with "protein-coding genes ... were matched".**
   The genes of ``--input_gff`` could not be linked to the protein IDs of the FASTA. ``gff_process`` (and
   ``easy_CGC``/``easy_substrate`` in protein mode, before step 1) stops when fewer than ``--gff_min_match``
   (default 0.9) of the protein-coding genes match. The message shows unmatched genes, the identifiers looked up
   and the ``--gff_type`` that matches the file best; ``gff_resolution.tsv`` lists every gene. Use the protein
   FASTA and the GFF of the same annotation run, and the suggested ``--gff_type``.

**DIAMOND fails with "Invalid character in sequence: '.'"**
   Proteins translated with ``gffread -y`` contain ``.`` for stop codons. In ``--mode protein`` run_dbcan now
   removes trailing stops and replaces internal ``*``/``.`` with ``X``; with other tools, translate with
   ``gffread -y ... -S``.

**Which protein IDs are safe?**
   Use unique IDs of at most 80 characters. In ``--mode protein``:

   * only the first word of each FASTA header is kept;
   * IDs longer than 80 characters are truncated;
   * duplicated IDs get ``_1``, ``_2``, ... appended;
   * IDs containing ``JGI`` have ``|`` replaced by ``-``.

   These changes are applied to ``uniInput.faa`` only, not to your GFF. For CGC prediction the GFF identifiers are
   cleaned the same way before matching, except for the ``_1``, ``_2`` suffixes of duplicated IDs.

**overview.tsv shows** ``-`` **everywhere for numeric IDs (5.2.9 and earlier).**
   Older releases read purely numeric IDs (``>0``, ``>1``) as numbers and ``NA``/``NULL``/``nan``/``None`` as
   missing values when merging the tables. This is fixed after 5.2.9 (GitHub issue #82). With an older release,
   prefix the IDs:

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
   ``--coverage_threshold_tc``, ``--coverage_threshold_tf_diamond`` and ``--coverage_threshold_sulfatase`` are
   passed to DIAMOND ``--query-cover`` (full-length protein databases) and use 0-100 (default 35).
   Sulfatase hits must also pass ``--subject_coverage_threshold_sulfatase`` (DIAMOND ``--subject-cover``, 0-100,
   default 70), because SulfAtlas entries are often sulfatases fused to other domains.
   ``--coverage_threshold_peptidase`` is passed to ``--subject-cover``: the peptidase database holds MEROPS
   peptidase units (domains), so the coverage is measured on the unit (0-100, default 35).
   The pyHMMER coverage options use 0-1 (default 0.35).

**Which unit does substrate_prediction -iden use?**
   Percent identity on a **0-100 scale** (DIAMOND ``pident``; older releases said 0.0-1.0 in ``--help``). ``-iden 0.3`` keeps
   almost everything; use ``-iden 30`` for 30% identity. ``-cov`` is a fraction (0-1) of the CGC protein length.

CGC prediction
--------------

**No CGC is found (cgc_standard_out.tsv is empty).**
   Check, in this order:

   1. ``cgc.gff`` exists. ``gff_process`` has no ``--mode`` option and defaults to ``<output_dir>/uniInput.gff``,
      which exists only after a ``prok``/``meta`` run; with protein input pass ``--input_gff`` and ``--gff_type``.
   2. Genes are matched to proteins. ``grep -c "protein_id=unknown" cgc.gff`` should be small and
      ``grep -c "CGC_annotation=CAZyme" cgc.gff`` should be greater than zero. Unmatched genes are written silently as
      ``protein_id=unknown``/``CGC_annotation=null``. The attribute used per ``--gff_type`` is listed in
      :ref:`cgc-information-generation`.
   3. Enough CAZymes are supported by two tools (``#ofTools`` >= 2 in ``overview.tsv``).
   4. The default rule requires at least one CAZyme **and** one transporter (``--additional_genes TC``) within
      ``--num_null_gene 2`` genes of each other. Relax it, for example with
      ``--additional_genes TC --additional_genes TF --additional_genes STP --additional_logic any``.

**--feature_type has no effect / "is deprecated and ignored".**
   ``cgc_finder`` uses every gene of ``cgc.gff``; the GFF feature type is chosen per ``--gff_type`` by
   ``gff_process``. The option is only kept so that old command lines still run. The same applies to
   ``--run_pfam``, ``--force_topology`` and the ``substrate_prediction`` options ``--pul``, ``-o``, ``-w``,
   ``-rerun`` and ``-env``.

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
   installed package (currently ``db_v5-2-9_5-5-2026``). Use ``--aws_s3`` for reproducible results and record the
   release name, the ``run_dbcan version`` output and the file checksums (see :doc:`../user_guide/prepare_the_database`).

**How much disk space do I need?**
   About 8 GB for all files (``dbCAN-sub.hmm`` ~5.1 GB, ``CAZy.dmnd`` ~2.2 GB, ``dbCAN.hmm`` ~130 MB, CGC databases
   ~0.5 GB). ``Pfam-A.hmm`` (~1.6 GB) is only needed for ``Pfam_null_cgc`` and must be downloaded manually.

**A download was interrupted.**
   Rerun the same command. Partial files (``*.part``) are resumed, complete files are skipped when their size matches,
   and ``--no-overwrite`` skips every file that already exists.

Memory and run time
-------------------

**How many threads are used?**
   ``--threads`` defaults to the number of logical CPUs of the machine and is used for gene prediction (pyrodigal),
   DIAMOND and pyHMMER. On shared servers or cluster jobs, set it to the number of cores you were allocated.
   ``substrate_prediction`` run on its own always uses all CPUs for DIAMOND.

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
   Both tools must be installed separately and are only run with ``--run_signalp`` / ``--run_deeptmhmm`` (each works
   on its own). A requested tool that cannot run makes ``CAZyme_annotation`` exit with a non-zero status; the log
   says which one. See :doc:`../user_guide/signalp6_deeptmhmm`.
