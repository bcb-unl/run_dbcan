.. _dbcan-utils-abundance:

Read Abundance with ``dbcan_utils``
===================================

``dbcan_utils`` turns reads mapped to the annotated contigs (metagenome or metatranscriptome) into abundances of
CAZyme families, subfamilies, ECs, substrates and CGCs, in two steps:

.. code-block:: shell

   samtools sort -o sample.bam sample.sam && samtools index sample.bam
   # 1. reads per gene
   dbcan_utils cal_coverage -g sample.dbCAN/cgc.gff -i sample.bam -o sample.depth.txt -t 6
   # 2. abundance tables (-i = run_dbcan output directory)
   dbcan_utils fam_abund           -bt sample.depth.txt -i sample.dbCAN -a TPM -o sample_abund
   dbcan_utils fam_substrate_abund -bt sample.depth.txt -i sample.dbCAN -a TPM -o sample_abund
   dbcan_utils CGC_abund           -bt sample.depth.txt -i sample.dbCAN -a TPM -o sample_abund
   dbcan_utils CGC_substrate_abund -bt sample.depth.txt -i sample.dbCAN -a TPM -o sample_abund

Every command exits with a non-zero status when it cannot write its output (missing input, BAM without index, no
GFF contig in the BAM, no gene ID shared by the count table and the run_dbcan tables).

Outputs
-------

.. list-table::
   :header-rows: 1
   :widths: 24 30 46

   * - Command
     - Output files
     - Content
   * - ``cal_coverage``
     - ``-o`` file
     - No header: gene ID, gene length (bp), read count
   * - ``fam_abund``
     - ``fam_abund.out``, ``subfam_abund.out``, ``EC_abund.out``
     - Summed abundance and number of proteins per CAZy family (dbCAN HMM + dbCAN-sub), subfamily and EC (every EC
       of ``EC#``)
   * - ``fam_substrate_abund``
     - ``fam_substrate_abund.out``
     - Summed abundance of the proteins whose dbCAN-sub hits name each substrate, with their IDs
   * - ``CGC_abund``
     - ``CGC_abund.out``
     - Mean abundance of the protein-coding genes of each CGC, with their IDs, abundances and ``Fams``
       (CAZy family, ``TC|3.A.1.1.27``-style label of the Gene Type role, or ``null``); tRNAs etc. are left out
   * - ``CGC_substrate_abund``
     - ``CGC_substrate_PUL_homology.out``, ``CGC_substrate_majority_voting.out``
     - Summed CGC abundance per substrate of the dbCAN-PUL and the dbCAN-sub prediction. A CGC with several
       substrates (``arabinan,galactan``) counts for each; a ``,`` between digits (``alpha-1,4-glucan``) is part of
       the name

``-o`` is created if needed (default: current directory). ``-a`` sets the normalisation, with *count* the read count
of a gene and *N* the sum over all genes:

* **RPM** (default) = count / (N / 10^6)
* **FPKM** (= **RPKM**) = RPM / (gene length / 1000)
* **TPM** = (count / length) / sum over genes (count / length) x 10^6

Choosing the GFF for ``cal_coverage``
-------------------------------------

The gene IDs of the count table must be the protein IDs of the run_dbcan tables.

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - run_dbcan input
     - GFF for ``cal_coverage``
   * - prokaryote / metagenome (``prok``, ``meta``, ``NCBI_prok``, ``prodigal``)
     - ``cgc.gff`` of the output directory (one row per gene = its CDS)
   * - eukaryote with NCBI, funannotate, BRAKER, GFF3 or GTF annotation
     - The **original annotation GFF**: NCBI with the default ``--id-attr auto``; funannotate/BRAKER with
       ``--id-attr Parent --feature-type CDS``. Not ``cgc.gff``: its gene rows include introns and UTRs
       (*Xylona*: median 1.37 x the CDS), which lowers FPKM/TPM and counts intronic reads
   * - eukaryote with JGI annotation
     - ``cgc.gff`` (the JGI ``proteinId`` values do not match the run_dbcan IDs); lengths then include introns

.. note::

   **How** ``--id-attr auto`` **names genes.** Rows with ``protein_id`` (NCBI, ``cgc.gff``) or ``proteinId`` (JGI) are
   named by it; otherwise by ``ID``, then ``Name``, then ``locus_tag`` (Prodigal, Prokka, Bakta). When ``CDS`` rows
   are present only they are counted, and rows with the same ID (exons) form one gene whose length is their union.
   ``protein_id=NA``/``unknown`` rows (tRNA, rRNA) are skipped. Prodigal GFFs (``ID=1_1``) need ``dbcan_utils
   gff_fix`` first.

Which reads are counted
-----------------------

An alignment counts for a gene when it is

* primary (secondary, supplementary, QC-fail and duplicate alignments are skipped; ``--include-duplicates`` keeps
  duplicates);
* MAPQ >= ``--mapping-quality`` (30);
* inside the gene for >= ``--overlap-base-ratio`` (0.2) of the read (``--hifi``: of the shorter of read and gene);
* identity >= ``--identity`` (0.98), with identity = 1 - edits / alignment columns, edits from ``NM``, else ``MD`` +
  insertions, else ``=``/``X`` CIGAR. Alignments without these tags are counted (with a warning; add them with
  ``samtools calmd``). ``--identity 0`` turns the filter off.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Option
     - Effect
   * - ``--count-mode fragments`` (default) / ``reads``
     - A read pair counts once per gene (``reads``: each mate, as in 5.2.9 and earlier). Changes raw counts only,
       not RPM/TPM. Use ``fragments`` for DESeq2/edgeR. A pair spanning two genes counts for both.
   * - ``--strandedness none`` (default) / ``reverse`` / ``forward``
     - Strand of read 1 relative to the gene: ``reverse`` for dUTP/TruSeq stranded RNA-seq
   * - ``--id-attr``, ``--feature-type``
     - GFF attribute that names genes, and feature type to count (see above)

.. warning::

   Genes present in several identical copies (transposases, rRNA operons, recent paralogs) get MAPQ 0-1 with most
   aligners, so **0 reads** at ``--mapping-quality 30``. Use ``--mapping-quality 0`` if they matter (reads then
   count at their primary alignment only).

The log reports skipped alignments per reason (``secondary_or_supplementary``, ``qc_fail``, ``duplicate``,
``low_mapq``, ``low_identity``) and rejected alignment/gene pairs (``wrong_strand``, ``low_overlap``). Identity now
counts inserted bases as edits, slightly stricter than 5.2.9. The old options ``-1/--r1``, ``-2/--r2`` and
``--db-dir`` are accepted with a warning and ignored.

.. _abundance-charts:

Charts: ``dbcan_plot bar_plot`` and ``heatmap_plot``
----------------------------------------------------

Both charts take one ``*_abund`` file per sample (comma-separated, named with ``--samples``) and write a PDF
(``--pdf``). They keep ``--top`` rows (20): with two or more files, those with the largest absolute difference
between the **first two** files; ``--col``/``--value`` select rows by name instead. All options:
:ref:`dbcan-plot-options`.

.. code-block:: shell

   dbcan_plot bar_plot -i ctrl1_abund/fam_abund.out,treat1_abund/fam_abund.out --samples ctrl1,treat1 \
       --top 15 --pdf fam_bar.pdf -o figures
   dbcan_plot heatmap_plot -i ctrl1_abund/fam_abund.out,ctrl2_abund/fam_abund.out,treat1_abund/fam_abund.out \
       --samples ctrl1,ctrl2,treat1 --top 15 --pdf fam_heatmap.pdf -o figures

.. figure:: ../_static/img/figures/abund_fam_bar.png
   :width: 55%
   :alt: Bar chart of CAZyme family abundance (TPM) in two samples

   ``bar_plot`` of two samples (largest difference on top). ``--vertical-bar`` turns the bars.

.. figure:: ../_static/img/figures/abund_fam_heatmap.png
   :width: 45%
   :alt: Heatmap of CAZyme family abundance (TPM) in six samples

   ``heatmap_plot`` of six samples (white below 0.1). ``--show-abund`` prints values, ``--cluster-map`` clusters
   rows. Full commands: :ref:`figure_gallery`.
