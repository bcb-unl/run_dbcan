.. _cgc-information-generation:

CGC Signature Annotation (``gff_process``)
==========================================

Step 2 of the pipeline. A CAZyme gene cluster (CGC) is a group of neighbouring genes that work on one glycan:
CAZymes plus "signature" genes such as transporters. ``gff_process`` finds the signature genes and writes
``cgc.gff``, your gene annotation with one CAZyme/signature label per gene, which ``cgc_finder`` (step 3) reads.

Run it in the output directory of a finished ``CAZyme_annotation`` (it reads ``uniInput.faa`` and ``overview.tsv``):

.. code-block:: shell

   # DNA input (prok/meta): the pyrodigal GFF uniInput.gff is used automatically
   run_dbcan gff_process --output_dir output_dir --db_dir db

   # protein input: the GFF of the same annotation and its type
   run_dbcan gff_process --output_dir output_dir --db_dir db --input_gff EscheriaColiK12MG1655.gff --gff_type NCBI_prok
   run_dbcan gff_process --output_dir output_dir --db_dir db --input_gff Xylhe1_GeneCatalog_proteins_20130827.gff --gff_type JGI --fungi

``easy_CGC`` and ``easy_substrate`` run this step for you.

What it searches
----------------

.. list-table::
   :widths: 16 34 50
   :header-rows: 1

   * - Label
     - Search
     - Kind of hit
   * - ``TC``
     - DIAMOND vs TCDB
     - full-length protein
   * - ``TF``
     - DIAMOND vs ``TF.dmnd`` (default); pyHMMER ``TF.hmm`` with ``--fungi``
     - full-length (DIAMOND) / domain (HMM)
   * - ``STP``
     - pyHMMER ``STP.hmm``
     - domain
   * - ``Sulfatase``
     - DIAMOND vs SulfAtlas
     - full-length protein
   * - ``Peptidase``
     - DIAMOND vs MEROPS peptidase units
     - domain (unit)
   * - ``CAZyme``
     - ``overview.tsv``: proteins with ``#ofTools`` >= 2, named by ``Recommend Results``
     - domains

Each DIAMOND search keeps the best hit per protein (highest bit score). The labels are then combined per gene
(``+``-joined in ``cgc.gff``), as drawn in :ref:`the CGC rule figure <cgc-rules-figure>`:

.. list-table::
   :widths: 45 55
   :header-rows: 1

   * - Protein has
     - Labels in ``cgc.gff``
   * - CAZyme (``#ofTools`` >= 2) + TC / TF (DIAMOND) / sulfatase hit
     - CAZyme only (a full-length hit cannot be placed next to the CAZyme domain)
   * - no CAZyme, several full-length hits (TC, TF, sulfatase)
     - the one with the highest **bit score** (ties: min(query, subject) coverage, E-value, TC > TF > Sulfatase)
   * - domain labels (STP, fungal TF HMM, peptidase)
     - always kept, also next to a CAZyme or a full-length label (``TC|...+STP|...``); the gene still fills only
       **one** role in ``cgc_finder``, which also ignores a domain lying on the CAZyme domain

All hits stay in ``total_cgc_info.tsv``.

.. note::

   * **Why the bit score**: E-values grow with database size, and TCDB, ``TF.dmnd`` and SulfAtlas differ in size.
     *E. coli* LacI: ``LACI_ECOLI`` (TF, 677 bits) beats TCDB ``3.A.1.2.1`` (55 bits), so it is ``TF``. A
     ``total_cgc_info.tsv`` without the ``Bitscore`` column (older runs) keeps all labels.
   * **Why CAZyme wins**: a GH13 protein that also hits TC 8.A.9 (rBAT/SLC3, itself a GH13 homologue) would
     otherwise be CAZyme *and* transporter, and any glycosyltransferase cluster with a TCDB homologue a CGC.
   * **Sulfatase**: a hit must also cover >= 70% of the SulfAtlas protein
     (``--subject_coverage_threshold_sulfatase``; ``0`` = off), because many entries are sulfatases fused to GH
     domains (*E. coli*: 6 of 113 hits remain, CGCs and substrates unchanged).
   * **Peptidase**: coverage is measured on the MEROPS unit; hits to inhibitor units (``I*``) are never reported
     as ``Peptidase`` (see :ref:`rebuild-peptidase-db`).

Choosing ``--gff_type``
-----------------------

``cgc.gff`` is built by matching each gene of the GFF to a protein ID of ``uniInput.faa``. ``gff_process`` needs
``--gff_type`` for every ``--input_gff`` you give (case-insensitive); only the pyrodigal ``uniInput.gff`` defaults
to ``prodigal``. ``easy_CGC``/``easy_substrate`` with ``--mode protein`` need ``--input_gff`` and ``--gff_type``.

.. list-table::
   :widths: 15 30 55
   :header-rows: 1

   * - ``--gff_type``
     - Use for
     - Protein ID taken from
   * - ``prodigal``
     - ``uniInput.gff``, Prodigal/pyrodigal GFF
     - ``ID`` of each line (one line = one gene, file order)
   * - ``NCBI_prok``
     - NCBI prokaryote GFF3 (``*_genomic.gff``)
     - ``protein_id`` of the gene's ``CDS``; genes with an RNA child become ``protein_id=NA``
   * - ``NCBI_euk``
     - NCBI eukaryote GFF3
     - ``protein_id`` of the first ``CDS`` of the first ``mRNA``; non-mRNA genes become ``protein_id=NA``
   * - ``JGI``
     - JGI MycoCosm GFF3
     - ``proteinId``, compared with the 3rd field of ``jgi-Xylhe1-123456-name``
   * - ``funannotate``
     - funannotate / funannotate2 GFF3
     - ``mRNA`` ID (``FUN_000001-T1``); ``CDS`` and gene IDs also tried
   * - ``BRAKER``
     - BRAKER / AUGUSTUS GFF3
     - transcript ID (``g1.t1``, as in ``braker.aa``); ``g1.t1.CDS1``-style IDs also tried
   * - ``GFF3``
     - other gene → mRNA/transcript → CDS GFF3 (Helixer, Ensembl, MAKER, gffread)
     - first ID found in the FASTA: ``CDS`` ``protein_id``/``proteinId``/``Name``/``ID``, transcript ID, gene ID
   * - ``GTF``
     - GTF (``braker.gtf``, NCBI/Ensembl ``.gtf``)
     - ``protein_id``, then ``transcript_id``, then ``gene_id`` of the ``CDS`` rows

.. note::

   **How** ``funannotate``, ``BRAKER``, ``GFF3`` **and** ``GTF`` **are matched.**

   * IDs are also tried without common prefixes/suffixes (``cds-``, ``rna-``, ``gene-``, ``CDS:``, ``transcript:``,
     ``gene:``; ``.cds``, ``-CDS``, ``.CDS.<n>``, ``-P<n>``), and as version-less accession or JGI ``proteinId`` if
     that points to a single protein; IDs are cleaned like FASTA headers (first word, ``|`` → ``-`` for JGI, 80 chars).
   * Exact ``ID``-like matches are assigned before ``Name``/``locus_tag``, gene IDs and aliases.
   * One protein per gene: of several isoforms, the CAZyme (else another signature, else the longest CDS) is kept;
     all isoforms are listed in ``gff_resolution.tsv``.
   * Transcripts without gene rows (gffread) are grouped by ``geneID``/``gene_id``/``gene`` only when they share
     contig and strand and overlap, so two ``tnpA`` copies stay two genes.
   * Genes without CDS (tRNA, rRNA, pseudogenes) are kept as ``protein_id=NA`` null genes; a pseudogene (a
     ``pseudogene`` row, or ``CDS`` rows with ``pseudo=true`` and no ``protein_id``) gets
     ``CGC_annotation=Other|pseudogene``.

.. _gff-match-check:

Gene/protein match check
~~~~~~~~~~~~~~~~~~~~~~~~

Before any search, the genes of ``--input_gff`` are matched to the proteins (a few seconds) and a summary is
logged::

   GFF resolution (--gff_type funannotate): 8201/8201 protein-coding genes matched to a protein (100.0%); ...

If fewer than ``--gff_min_match`` (default ``0.9``) of the protein-coding genes match, the command **stops** and
names unmatched genes, the IDs it tried and the ``--gff_type`` that fits the file best. ``easy_CGC`` and
``easy_substrate`` in protein mode run this check before step 1, so a wrong GFF fails in seconds.
``--gff_min_match 0`` disables the stop. The check is skipped for the pyrodigal ``uniInput.gff``.

``gff_resolution.tsv`` lists every gene: ``gene_id``, ``contig``, ``start``, ``end``, ``strand``, ``protein_id``,
``matched_by`` (``~`` = after removing a prefix/suffix), ``status`` (``matched``, ``unmatched``, ``noncoding``),
``n_transcripts``, ``isoform_proteins``.

Fungi and other eukaryotes
--------------------------

run_dbcan does not predict eukaryotic genes. Annotate the genome with funannotate/funannotate2 or BRAKER (or use
the NCBI or JGI annotation), keep the protein FASTA and the GFF of the **same** run, and add ``--fungi``:

.. code-block:: shell

   # funannotate: <name>.proteins.fa + <name>.gff3
   run_dbcan easy_CGC --mode protein --input_raw_data Xylona.proteins.fa \
     --input_gff Xylona.gff3 --gff_type funannotate --fungi --output_dir output_fungus --db_dir db

   # BRAKER: braker.aa + braker.gff3 (or braker.gtf with --gff_type GTF)
   run_dbcan easy_CGC --mode protein --input_raw_data braker/braker.aa \
     --input_gff braker/braker.gff3 --gff_type BRAKER --fungi --output_dir output_fungus --db_dir db

   # a pipeline that writes only a GFF3 (e.g. Helixer): translate it first
   gffread annotation.gff3 -g genome.fna -y proteins.faa -S
   run_dbcan easy_CGC --mode protein --input_raw_data proteins.faa \
     --input_gff annotation.gff3 --gff_type GFF3 --fungi --output_dir output_fungus --db_dir db

.. note::

   The CGC rule was developed for bacteria. Fungal genes are about three times less dense and rarely form operons,
   so fungal CGCs are less specific (``--use_distance`` with 15 kb spans about 5 fungal genes), and dbCAN-PUL has
   no fungal PULs. Last resort for a GFF that no ``--gff_type`` matches: one line per gene with
   ``ID=<protein ID>`` and ``--gff_type prodigal`` (file order, no tRNA genes).

Options
-------

.. list-table::
   :widths: 38 62
   :header-rows: 1

   * - Option
     - Description (default)
   * - ``--output_dir`` / ``--db_dir``
     - Output folder of ``CAZyme_annotation``; database folder
   * - ``--input_gff``
     - Annotation GFF3/GTF (``<output_dir>/uniInput.gff``)
   * - ``--gff_type``
     - See the table above; required with any ``--input_gff`` (``prodigal`` for ``uniInput.gff``)
   * - ``--gff_min_match``
     - Minimum fraction of matched protein-coding genes (``0.9``; ``0`` = no stop)
   * - ``--fungi`` / ``--no-fungi``
     - TF HMM search instead of TF DIAMOND (off)
   * - ``--prokaryotic`` / ``--no-prokaryotic``
     - TF DIAMOND search (on without ``--fungi``; ``--prokaryotic --fungi`` runs both)
   * - ``--coverage_threshold_tc``, ``_tf_diamond``, ``_sulfatase``
     - DIAMOND ``--query-cover``, percent 0-100 (35)
   * - ``--subject_coverage_threshold_sulfatase``
     - DIAMOND ``--subject-cover`` for SulfAtlas, percent (70)
   * - ``--coverage_threshold_peptidase``
     - Coverage of the MEROPS unit (``--subject-cover``), percent (35)
   * - ``--coverage_threshold_stp``, ``_tf``
     - HMM coverage, fraction 0-1 (0.35)
   * - ``--e_value_threshold_*``
     - E-value of each search (``1e-4``)

``gff_process`` has no ``--mode`` option. Without ``--input_gff`` it needs the ``uniInput.gff`` of a
``prok``/``meta`` run.

Output files
------------

.. list-table::
   :widths: 35 65
   :header-rows: 1

   * - File
     - Content
   * - ``cgc.gff``
     - One ``gene`` row per gene with ``protein_id`` and ``CGC_annotation`` (e.g. ``CAZyme|GH13_31``,
       ``TC|2.A.1.5.1``, ``null``): the input of ``cgc_finder``
   * - ``total_cgc_info.tsv``
     - All signature hits; ``Source`` = search (``TC``, ``TF_prok``, ``TF_fungi``, ``STP``, ``Sulfatase``,
       ``Peptidase``), ``Bitscore`` = DIAMOND bit score (last column, empty for HMMs)
   * - ``diamond.out.tc``, ``.tf``, ``.sulfatase``, ``.peptidase``
     - Raw DIAMOND hits (``BITSCORE`` after ``Database``)
   * - ``STP_hmm_results.tsv``, ``TF_hmm_results.tsv``
     - pyHMMER hits (TF only with ``--fungi``)
   * - ``gff_resolution.tsv``
     - Gene → protein match of ``--input_gff`` (see :ref:`gff-match-check`)

Next: :doc:`CGC_annotation` (step 3).
