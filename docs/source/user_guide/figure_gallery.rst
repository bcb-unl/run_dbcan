.. _figure_gallery:

Figure Gallery
==============

Every figure that run_dbcan, ``dbcan_plot`` and ``dbcan_utils`` can draw, with the command that made it. How to
read each figure and all options: :ref:`cgc_visualization` and :ref:`abundance-charts`.

.. list-table:: Example data (real runs of the current code)
   :widths: 30 70
   :header-rows: 1

   * - Output directory
     - Input
   * - ``output_EscheriaColiK12MG1655_fna``
     - *E. coli* K-12 MG1655, one chromosome ``NC_000913.3`` with its NCBI GFF (``--gff_type NCBI_prok``):
       52 CGCs, 16 with a substrate
   * - ``output_Xylona_heveae_TC161_faa``
     - *Xylona heveae* TC161 (fungus, 23 scaffolds), NCBI proteins + GFF (``--gff_type NCBI_euk --fungi``):
       110 CGCs on 15 scaffolds, 6 with a substrate
   * - ``output_PUL0008``
     - Genomic region of dbCAN-PUL ``PUL0008`` (fructan) run as a genome (``--gff_type prodigal``)
   * - ``ctrl1``-``ctrl3``, ``treat1``-``treat3``
     - Simulated paired-end RNA-seq of *E. coli* mapped to ``NC_000913.3``; some CGCs, among them CGC7 (*lac*),
       are expressed four times lower in the treatment

The images are PNG at 150 dpi; the commands write PDF (and PNG at 200 dpi with ``--format png``).

Genome view: circos (one complete genome)
-----------------------------------------

.. figure:: ../_static/img/figures/cgc_circos_ecoli.png
   :width: 80%
   :alt: Circos plot of the CGCs of E. coli K-12 MG1655

   The 52 CGCs of *E. coli*. Rings explained in :ref:`cgc-genome-view`.

.. code-block:: shell

   run_dbcan cgc_circle_plot --output_dir output_EscheriaColiK12MG1655_fna
   # -> output_EscheriaColiK12MG1655_fna/cgc_circos/cgc_circos_plot.png and .pdf

Genome view: circos with differential expression
------------------------------------------------

.. figure:: ../_static/img/figures/cgc_circos_ecoli_deg.png
   :width: 80%
   :alt: Circos plot of E. coli with up- and down-regulated genes

   The same plot with ``DEG.tsv``: up- (red) and down-regulated (blue) genes and the log2 fold change. Here
   ``DEG.tsv`` lists the 224 genes with \|log2FC\| >= 1.5 (treatment vs control, mean TPM, genes with TPM >= 20 in
   one condition); usually it is filtered DESeq2 or edgeR output.

.. code-block:: shell

   cp DEG.tsv output_EscheriaColiK12MG1655_fna/
   run_dbcan cgc_circle_plot --output_dir output_EscheriaColiK12MG1655_fna

Genome view: overview (many contigs)
------------------------------------

.. figure:: ../_static/img/figures/cgc_overview_xylona.png
   :width: 80%
   :alt: Overview of the CGCs on the scaffolds of Xylona heveae

   One bar per scaffold with a CGC; ``--layout auto`` picks this layout for several contigs.

.. code-block:: shell

   run_dbcan cgc_circle_plot --output_dir output_Xylona_heveae_TC161_faa
   # -> output_Xylona_heveae_TC161_faa/cgc_circos/cgc_overview.png and .pdf

CGC map: ``CGC_plot``
---------------------

.. figure:: ../_static/img/figures/cgc_map_ecoli_CGC7.png
   :width: 100%
   :alt: Gene map of E. coli CGC7

   CGC7 of *E. coli* (lactose). See :ref:`cgc-map`.

.. code-block:: shell

   dbcan_plot CGC_plot -i output_EscheriaColiK12MG1655_fna --cgcid 'NC_000913.3|CGC7' -o figures --format pdf,png
   # -> figures/NC_000913.3_CGC7.cgc.pdf and .png

CGC map with read depth: ``CGC_coverage_plot``
----------------------------------------------

.. figure:: ../_static/img/figures/cgc_coverage_ecoli_CGC7.png
   :width: 100%
   :alt: Read depth of two samples above the gene map of E. coli CGC7

   One panel per sample on a shared y axis: CGC7 is lower in ``treat1``. See :ref:`cgc-coverage`.

.. code-block:: shell

   samtools depth -a -r 'NC_000913.3:352706-366305' ctrl1.bam  > ctrl1.CGC7.depth.txt
   samtools depth -a -r 'NC_000913.3:352706-366305' treat1.bam > treat1.CGC7.depth.txt
   dbcan_plot CGC_coverage_plot -i output_EscheriaColiK12MG1655_fna --cgcid 'NC_000913.3|CGC7' \
       --reads-count ctrl1.CGC7.depth.txt --reads-count treat1.CGC7.depth.txt --samples ctrl1,treat1 \
       -o figures --format pdf,png
   # -> figures/NC_000913.3_CGC7.cgc-coverage.pdf and .png

Synteny with dbCAN-PUL: ``CGC_synteny_plot``
--------------------------------------------

.. figure:: ../_static/img/figures/cgc_synteny_ecoli_CGC7_PUL0751.png
   :width: 100%
   :alt: E. coli CGC7 compared with dbCAN-PUL PUL0751

   CGC7 above its best dbCAN-PUL hit, PUL0751 (lactose). See :ref:`cgc-synteny`.

.. code-block:: shell

   dbcan_plot CGC_synteny_plot -i output_EscheriaColiK12MG1655_fna --cgcid 'NC_000913.3|CGC7' --db-dir db \
       -o figures --format pdf,png
   # -> figures/NC_000913.3_CGC7_PUL0751-syntenic.pdf and .png

Synteny with read depth: ``CGC_synteny_coverage_plot``
------------------------------------------------------

.. figure:: ../_static/img/figures/cgc_synteny_coverage_ecoli_CGC7.png
   :width: 100%
   :alt: Read depth panels above the synteny plot of E. coli CGC7 and PUL0751

   The synteny figure with the depth panels of ``CGC_coverage_plot`` above the CGC.

.. code-block:: shell

   dbcan_plot CGC_synteny_coverage_plot -i output_EscheriaColiK12MG1655_fna --cgcid 'NC_000913.3|CGC7' --db-dir db \
       --reads-count ctrl1.CGC7.depth.txt --reads-count treat1.CGC7.depth.txt --samples ctrl1,treat1 \
       -o figures --format pdf,png
   # -> figures/NC_000913.3_CGC7_PUL0751-syntenic-cov.pdf and .png

Synteny PDFs of ``substrate_prediction``
----------------------------------------

.. figure:: ../_static/img/figures/synteny_PUL0008_CGC1.png
   :width: 100%
   :alt: CGC1 of the PUL0008 region compared with dbCAN-PUL PUL0008

   One of the ``synteny_pdf/`` files, as PNG: the PUL0008 region finds itself; the paralogous ``GH32`` genes also
   hit each other's counterpart (multiple matches).

.. code-block:: shell

   run_dbcan substrate_prediction --output_dir output_PUL0008 --db_dir db
   # -> output_PUL0008/synteny_pdf/PUL0008_GCMeta_00079453_C0_CGC1_PUL0008-syntenic.pdf

Abundance bar chart: ``bar_plot``
---------------------------------

.. figure:: ../_static/img/figures/abund_fam_bar.png
   :width: 60%
   :alt: Bar chart of CAZyme family abundance (TPM) in two samples

   ``fam_abund.out`` (TPM) of two samples: the families that differ most. See :ref:`abundance-charts`.

.. code-block:: shell

   for s in ctrl1 treat1; do
     dbcan_utils cal_coverage -g output_EscheriaColiK12MG1655_fna/cgc.gff -i $s.bam -o $s.depth.txt -t 2
     dbcan_utils fam_abund -bt $s.depth.txt -i output_EscheriaColiK12MG1655_fna -a TPM -o ${s}_abund
   done
   dbcan_plot bar_plot -i ctrl1_abund/fam_abund.out,treat1_abund/fam_abund.out --samples ctrl1,treat1 \
       --top 15 --pdf fam_bar.pdf -o figures
   # -> figures/fam_bar.pdf

Abundance heatmap: ``heatmap_plot``
-----------------------------------

.. figure:: ../_static/img/figures/abund_fam_heatmap.png
   :width: 50%
   :alt: Heatmap of CAZyme family abundance (TPM) in six samples

   ``fam_abund.out`` (TPM) of six samples (white below 0.1); rows chosen by the first two files.

.. code-block:: shell

   dbcan_plot heatmap_plot \
       -i ctrl1_abund/fam_abund.out,ctrl2_abund/fam_abund.out,ctrl3_abund/fam_abund.out,treat1_abund/fam_abund.out,treat2_abund/fam_abund.out,treat3_abund/fam_abund.out \
       --samples ctrl1,ctrl2,ctrl3,treat1,treat2,treat3 --top 15 --pdf fam_heatmap.pdf -o figures
   # -> figures/fam_heatmap.pdf
