.. _cgc-information-generation:

CAZyme Gene Cluster (CGC) Annotation
====================================

Introduction
-------------

A CAZyme Gene Cluster (CGC) refers to a group of genes co-located on the genome that collectively participate in glycan metabolism. These clusters encode various enzymes and regulatory proteins that work together to process specific carbohydrate substrates. CGCs are particularly important in microbial genomes, where they enable efficient utilization of diverse carbohydrate sources in the environment.

Generating CGC Information
---------------------------

To identify and analyze CGCs, users must first prepare annotation information by converting their GFF file into a specialized CGC-ready format. The ``gff_process`` command must be run in the output directory of a finished ``CAZyme_annotation`` run (it reads ``uniInput.faa`` and ``overview.tsv``) and performs these steps:

1. Uses DIAMOND to annotate transporters (TC, against TCDB), sulfatases (SulfAtlas), peptidases (MEROPS) and, for
   prokaryotes, transcription factors (TF, against ``TF.dmnd``)
2. Uses pyHMMER to identify signal transduction proteins (STP) and, with ``--fungi``, fungal transcription factors
   (``TF.hmm``; the prokaryotic TF DIAMOND search is then skipped)
3. Combines all hits that pass their E-value and coverage thresholds into ``total_cgc_info.tsv``

.. note::

   The DIAMOND databases contain full-length proteins, so a TC, TF (DIAMOND), sulfatase or peptidase hit is a
   protein-level label; its query interval is an alignment span, not a domain, so it cannot be compared with a
   CAZyme domain. On a protein that is a CGC CAZyme (``#ofTools`` >= 2) the CAZyme annotation therefore wins: its
   TC, TF (DIAMOND) and sulfatase labels are not written to ``cgc.gff`` (they stay in ``total_cgc_info.tsv``). For
   example a GH13 protein that also hits TC 8.A.9 (the rBAT/SLC3 family of TCDB, itself a GH13 homologue) is a
   CAZyme only, and a cluster of glycosyltransferases cannot become a CGC because one of them is also in TCDB.
   Domain-level labels (STP, fungal TF HMM, MEROPS peptidase units) are kept, and CGCFinder ignores them only when
   they lie on the CAZyme domain. On proteins without a CAZyme annotation, domain-level labels are kept side by side
   with one full-length label (one protein can be, for example, ``TC`` and ``STP``); such a gene still fills only
   one role when CGCs are called. Each DIAMOND search keeps one best hit per protein (highest bit score); STP and
   TF HMM domains are resolved within their own HMM library. The peptidase database is built from MEROPS ``pepunit.lib`` (peptidase units); hits to MEROPS
   peptidase *inhibitor* units (families ``I*``) are not reported as ``Peptidase`` (they stay in the raw
   ``diamond.out.peptidase`` file). In ``cgc.gff`` all labels of a protein are joined with ``+``; the ``Gene Type`` column of
   ``cgc_standard_out.tsv`` shows the first of CAZyme > TC > TF > STP > Sulfatase > Peptidase.

   **One full-length label per protein.** A protein that is not a CGC CAZyme can hit more than one full-length
   database, for example TCDB and SulfAtlas, or TCDB and the prokaryotic TF database. These hits are competing
   names for the whole protein, not domains at different positions, so only the best one is written to
   ``cgc.gff`` (all of them stay in ``total_cgc_info.tsv``). The best hit is the one with the highest DIAMOND
   **bit score** (``Bitscore``, the last column of ``total_cgc_info.tsv``):

   * the E-value cannot be used: it grows with the size of the database searched (E = K·m·n·e\ :sup:`-λS`, with
     *n* the total length of the database), and TCDB, ``TF.dmnd`` and SulfAtlas differ in size, so the same
     alignment gets a different E-value in each of them;
   * the bit score is the raw alignment score normalised by the scoring system (BLOSUM62 and the same gap
     penalties for all three searches) and does not depend on database size, so the scores of one query protein
     against the three databases are on the same scale;
   * the three searches share the query, so ranking by bit score is the same as ranking by the blast score ratio
     (bit score / bit score of the query aligned with itself): the self-score is one number per query protein.

   Ties (equal bit scores) go to the higher min(query coverage, subject coverage) of the alignment, then to the
   lower E-value, then to TC > TF > Sulfatase. For example *E. coli* LacI hits ``LACI_ECOLI`` in the TF database
   (bit score 677) and TCDB ``3.A.1.2.1`` (55; the ribose ABC transporter, whose binding protein shares the LacI
   ligand-binding fold): it is labelled ``TF`` only. Domain-level labels (STP, fungal TF HMM, MEROPS peptidase units)
   are not part of this choice. A ``total_cgc_info.tsv`` written before the ``Bitscore`` column was added keeps
   all labels.

   SulfAtlas entries are full-length proteins, and many of them are sulfatases fused to other domains (for example
   GH2 or GH32). A protein that aligns only to the non-sulfatase part of such an entry would take its sulfatase
   family (a GH32 protein labelled ``S1_7``). A sulfatase hit must therefore also cover most of the SulfAtlas
   protein: ``--subject_coverage_threshold_sulfatase`` (DIAMOND ``--subject-cover``, percent, default 70; 0 turns
   it off), in addition to ``--coverage_threshold_sulfatase`` (``--query-cover``, default 35). On *E. coli* K-12
   this keeps 6 of 113 sulfatase hits and on *Xylona heveae* 10 of 267; the CGCs and their substrate predictions
   do not change.

4. Marks CAZymes using ``overview.tsv``: only proteins with ``#ofTools`` >= 2 count as CAZymes, annotated with their ``Recommend Results``
5. Integrates these annotations into the user-submitted GFF file to generate a ``cgc.gff`` file

``easy_CGC`` and ``easy_substrate`` run this step automatically.

Input GFF File Types
-----------------------

``cgc.gff`` is built by matching each gene in your GFF to a protein ID in ``uniInput.faa``. How a gene is linked
to its protein depends on ``--gff_type`` (case-insensitive; any other value is rejected). With ``--mode protein``
``--gff_type`` must be given; ``prodigal`` is only the default for the ``uniInput.gff`` written by pyrodigal in
``prok``/``meta`` mode.

.. list-table::
   :widths: 16 30 54
   :header-rows: 1

   * - ``--gff_type``
     - Use for
     - How the gene is linked to its protein
   * - ``prodigal``
     - ``uniInput.gff`` of ``prok``/``meta`` mode; Prodigal/pyrodigal GFF
     - Every feature line is a gene; the protein ID is its ``ID``. Features are kept in file order.
   * - ``NCBI_prok``
     - NCBI prokaryote GFF3 (``*_genomic.gff``)
     - ``gene`` features; ``protein_id`` of the gene's ``CDS`` child. Genes with a non-CDS child (tRNA, rRNA, ...)
       are written as ``protein_id=NA`` with ``CGC_annotation=Other|<type>``.
   * - ``NCBI_euk``
     - NCBI eukaryote GFF3
     - ``gene`` features with ``mRNA`` children; ``protein_id`` of the first ``CDS`` under the first ``mRNA`` that
       has one. Genes with a child that is not ``mRNA`` (tRNA, ncRNA, ...) are written as ``protein_id=NA``.
   * - ``JGI``
     - JGI MycoCosm GFF3
     - ``proteinId`` attribute of the gene, compared with the third ``-``-separated field of the FASTA IDs
       (``jgi|Xylhe1|123456|name`` becomes ``jgi-Xylhe1-123456-name`` in protein mode).
   * - ``funannotate``
     - funannotate and funannotate2 GFF3
     - ``gene`` → ``mRNA`` → ``CDS``; proteins are named by the ``mRNA`` ID (``FUN_000001-T1``, or the locus tag
       after ``funannotate annotate``); the ``CDS`` ID (``FUN_000001-T1.cds``) and gene ID are also tried.
   * - ``BRAKER``
     - BRAKER / AUGUSTUS GFF3 (``braker.gff3``, ``augustus --gff3=on``)
     - ``gene`` → ``mRNA`` or ``transcript`` → ``CDS``; proteins are named by the transcript ID (``g1.t1``, as in
       ``braker.aa``); ``CDS`` IDs such as ``g1.t1.CDS1`` or ``g1.t1.cds`` are also tried.
   * - ``GFF3``
     - any other ``gene`` → ``mRNA``/``transcript`` → ``CDS`` GFF3 (Helixer, Ensembl, MAKER, gffread, ...)
     - In this order, the first identifier found in the FASTA: ``CDS`` ``protein_id``/``proteinId``/``Name``/``ID``,
       transcript ``ID``/``transcript_id``, gene ``ID``/``Name``/``locus_tag``.
   * - ``GTF``
     - GTF (``braker.gtf``, AUGUSTUS GTF, NCBI/Ensembl ``.gtf``)
     - ``CDS`` rows grouped by ``transcript_id``; ``protein_id``, then ``transcript_id``, then ``gene_id``.

For ``funannotate``, ``BRAKER``, ``GFF3`` and ``GTF`` (the resolver):

* every identifier is also tried without the prefixes ``cds-``, ``rna-``, ``gene-``, ``CDS:``, ``transcript:``,
  ``gene:`` and without the suffixes ``.cds``, ``-CDS``, ``.CDS.<n>``, ``-P<n>``; a version-less accession
  (``XP_1`` for ``XP_1.1``) or a JGI ``proteinId`` is used only if it points to a single FASTA protein;
* identifiers go through the same ID cleaning as the FASTA headers in ``--mode protein`` (first word, ``|`` → ``-``
  for JGI IDs, at most 80 characters);
* exact identifiers win over weak ones whatever the order of the genes in the file: proteins are first assigned
  through ``ID``-like attributes (``CDS`` ``protein_id``/``proteinId``/``ID``, transcript ``ID``/``transcript_id``)
  that equal a FASTA ID, and only then through ``Name``/``locus_tag``, gene-level IDs and version-less/JGI aliases;
* transcripts without a gene parent (gffread output without gene rows) are grouped into genes by their ``geneID``,
  ``gene_id`` or ``gene`` attribute, so that their isoforms give one gene. Only transcripts on the same contig and
  strand that overlap are grouped: ``gene=`` often holds a gene symbol, and two copies of ``tnpA`` or two ``trnL``
  tRNAs stay two genes. A repeated value names each locus after its first transcript; a gene's coordinates span all
  its transcripts;
* one protein is kept per gene. When a gene has several transcripts (isoforms), the one annotated as a CAZyme
  (else as another CGC signature) is kept, otherwise the one with the longest CDS. All isoform proteins of the gene
  are listed in ``gff_resolution.tsv``;
* genes without a CDS (tRNA, rRNA, ncRNA, pseudogenes) are kept as ``protein_id=NA`` rows, as with ``NCBI_euk``,
  so that they count as null genes between signature genes. Pseudogenes (a ``pseudogene`` row, or a gene whose
  ``CDS`` rows carry ``pseudo=true`` and no ``protein_id``) are such rows too (``CGC_annotation=Other|pseudogene``)
  and do not count as unmatched protein-coding genes: on the NCBI GFF3 of *E. coli* K-12, ``--gff_type GFF3``
  matches 4140/4140 protein-coding genes, as ``NCBI_prok`` does;
* ``cgc.gff`` has one ``gene`` row per gene, contigs in file order, genes sorted by position. On the NCBI GFF3 of
  *Xylona heveae*, ``--gff_type GFF3`` writes a ``cgc.gff`` identical to ``--gff_type NCBI_euk``.

Gene/protein match check
~~~~~~~~~~~~~~~~~~~~~~~~

For every ``--gff_type``, ``gff_process`` first matches the genes of ``--input_gff`` to the protein IDs (this
takes a few seconds) and writes ``gff_resolution.tsv`` (one row per gene: ``gene_id``, ``contig``, ``start``,
``end``, ``strand``, ``protein_id``, ``matched_by`` (the attribute that matched, ``~`` = after removing a prefix or
suffix), ``status`` (``matched``, ``unmatched``, ``noncoding``), ``n_transcripts`` and ``isoform_proteins``). A
summary line is logged, for example::

   GFF resolution (--gff_type funannotate): 8201/8201 protein-coding genes matched to a protein (100.0%);
   145 non-coding genes kept as null rows; ...

If fewer than ``--gff_min_match`` (default ``0.9``) of the protein-coding genes are matched, the command stops
with a non-zero exit status before the TC/TF/STP/sulfatase/peptidase searches. The error lists unmatched genes
with the identifiers that were looked up, the first FASTA IDs, and the ``--gff_type`` that matches most genes of
the file. In ``easy_CGC`` and ``easy_substrate`` with ``--mode protein`` the check runs before step 1, so a wrong
``--gff_type`` or GFF fails in seconds. ``--gff_min_match 0`` turns the stop off (the table is still written).
The check is skipped for the ``uniInput.gff`` that pyrodigal writes in ``prok``/``meta`` mode.

Eukaryotic genomes (fungi, etc.)
---------------------------------

run_dbcan does not predict eukaryotic genes, and ``--mode prok``/``meta`` must not be used for them: pyrodigal is a
prokaryotic gene finder and splits intron-containing genes into short ORFs (a ``WARNING`` is logged when the
predicted proteins look eukaryotic). For fungi and other eukaryotes:

1. Annotate the genome with a eukaryotic pipeline, for example funannotate / funannotate2 or BRAKER, or download
   the NCBI or JGI MycoCosm annotation. Keep the protein FASTA and the GFF3/GTF of the **same** run.
2. Run ``easy_CGC`` (or ``easy_substrate``) with ``--mode protein`` on the proteins, ``--input_gff`` with the GFF and
   the matching ``--gff_type``. Add ``--fungi`` to search fungal transcription factors with ``TF.hmm`` instead of
   the prokaryotic TF DIAMOND search.

.. code-block:: shell

   # funannotate / funannotate2: <name>.proteins.fa + <name>.gff3
   run_dbcan easy_CGC --mode protein --input_raw_data Xylona.proteins.fa \
     --input_gff Xylona.gff3 --gff_type funannotate --fungi \
     --output_dir output_fungus --db_dir db

   # BRAKER: braker.aa + braker.gff3 (or braker.gtf with --gff_type GTF)
   run_dbcan easy_CGC --mode protein --input_raw_data braker/braker.aa \
     --input_gff braker/braker.gff3 --gff_type BRAKER --fungi \
     --output_dir output_fungus --db_dir db

   # NCBI or JGI annotation
   run_dbcan easy_CGC --mode protein --input_raw_data Xylona_heveae_TC161.faa \
     --input_gff Xylona_heveae_TC161.gff --gff_type NCBI_euk --fungi \
     --output_dir output_fungus --db_dir db

If your pipeline only writes a GFF3 (for example Helixer), translate it with gffread and use ``--gff_type GFF3``;
``-S`` writes ``*`` instead of ``.`` for stop codons (run_dbcan replaces both, see :doc:`CAZyme_annotation`):

.. code-block:: shell

   gffread annotation.gff3 -g genome.fna -y proteins.faa -S
   run_dbcan easy_CGC --mode protein --input_raw_data proteins.faa \
     --input_gff annotation.gff3 --gff_type GFF3 --fungi --output_dir output_fungus --db_dir db

.. note::

   Earlier versions of this guide suggested adding ``protein_id`` to the ``CDS`` rows with ``awk`` and using
   ``--gff_type NCBI_euk``. This is no longer needed. If a file is not matched by any ``--gff_type`` (the error
   message says so), a GFF with one line per gene and ``ID=<protein ID of the FASTA>`` in column 9 can be used
   with ``--gff_type prodigal`` as a last resort; genes are then taken in file order and tRNA genes are not kept.

.. note::

   The CGC rule (CAZyme plus signature genes separated by at most ``--num_null_gene`` genes) was developed for
   bacterial genomes. Fungal genes are about three times less dense than bacterial genes and rarely form operons,
   so fungal CGCs are less specific; ``--use_distance`` with 15 kb spans about 5 fungal genes. dbCAN-PUL contains
   no fungal PULs, so PUL-based substrate prediction for fungal CGCs relies on bacterial references.

Command Examples
------------------

.. code-block:: shell

   run_dbcan gff_process --output_dir output_dir --db_dir db --input_gff EscheriaColiK12MG1655.gff --gff_type NCBI_prok

.. code-block:: shell

   run_dbcan gff_process --output_dir output_dir --db_dir db --input_gff output_dir/uniInput.gff --gff_type prodigal

.. code-block:: shell

   run_dbcan gff_process --output_dir output_dir --db_dir db --input_gff Xylhe1_GeneCatalog_proteins.gff --gff_type JGI

.. code-block:: shell

   run_dbcan gff_process --output_dir output_dir --db_dir db --input_gff eukaryotic_genome.gff --gff_type NCBI_euk --fungi

.. code-block:: shell

   run_dbcan gff_process --output_dir output_dir --db_dir db --input_gff braker.gff3 --gff_type BRAKER --fungi

.. note::

   ``gff_process`` has no ``--mode`` option. If ``--input_gff`` is omitted, it uses ``<output_dir>/uniInput.gff`` with
   ``--gff_type prodigal``, which only exists after a ``prok``/``meta`` run. Any other ``--input_gff`` needs
   ``--gff_type``. A missing GFF file or too few matched genes stop the command with an error.

Key Parameters
~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - Parameter
     - Description
   * - ``--output_dir``
     - Directory for both input and output files, it is generated from the ``CAZyme_annotation`` step
   * - ``--db_dir``
     - Directory containing database files, it is generated from the ``database`` step
   * - ``--input_gff``
     - Path to your GFF annotation file (default: ``<output_dir>/uniInput.gff``)
   * - ``--gff_type``
     - Format of the GFF file: ``prodigal``, ``NCBI_prok``, ``NCBI_euk``, ``JGI``, ``funannotate``, ``BRAKER``,
       ``GFF3`` or ``GTF`` (see `Input GFF File Types`_). Required for any ``--input_gff`` other than
       ``<output_dir>/uniInput.gff`` (default there: ``prodigal``)
   * - ``--gff_min_match``
     - Stop when fewer than this fraction of the protein-coding genes are matched to a protein (default ``0.9``;
       ``0`` disables the stop)
   * - ``--fungi`` / ``--no-fungi``
     - Fungal input: annotate TFs with the TF HMM search (default: off)
   * - ``--prokaryotic`` / ``--no-prokaryotic``
     - Run the prokaryotic TF DIAMOND search (default: on without ``--fungi``, off with ``--fungi``;
       pass ``--prokaryotic --fungi`` to run both)

E-value and coverage thresholds for each search (``--e_value_threshold_tc``, ``--coverage_threshold_tc``, ...)
are listed by ``run_dbcan gff_process --help``. DIAMOND coverage options are percentages (0-100) of the query
length, except for peptidases, where the coverage is measured on the matched MEROPS peptidase unit (subject);
pyHMMER coverage options are fractions of the HMM length (0-1).

Output Files
-------------

This step generates several output files:

* ``diamond.out.tc`` - DIAMOND output for Transporters (TCs) via TCDB database (https://www.tcdb.org/)
* ``diamond.out.tf`` - DIAMOND output for prokaryotic Transcription Factors (TFs)
* ``TF_hmm_results.tsv`` - pyHMMER output for fungal Transcription Factors (only with ``--fungi``)
* ``STP_hmm_results.tsv`` - pyHMMER output for Signal Transduction Proteins (STPs) using the STP-HMM database
* ``diamond.out.sulfatase`` and ``diamond.out.peptidase`` - DIAMOND output for sulfatases and peptidases
* ``total_cgc_info.tsv`` - Combined table of all non-CAZyme signature genes (TC, TF, STP, Sulfatase, Peptidase).
  ``Source`` names the search (``TC``, ``TF_prok``, ``TF_fungi``, ``STP``, ``Sulfatase``, ``Peptidase``) and
  ``Bitscore`` (last column) the DIAMOND bit score (empty for the HMM searches). The ``diamond.out.*`` files have a
  ``BITSCORE`` column after ``Database``.
* ``cgc.gff`` - Simplified GFF with ``protein_id`` and ``CGC_annotation`` attributes, the input of ``cgc_finder``
* ``gff_resolution.tsv`` - Gene → protein match of every gene of ``--input_gff`` (not written for the pyrodigal
  ``uniInput.gff``); see `Gene/protein match check`_

.. admonition:: Next Steps

   After generating the CGC annotation information, proceed to :doc:`CGC identification <CGC_annotation>` to identify and analyze CAZyme gene clusters in your genome.
