Metagenomics protocol: Wastyk 2021 (human gut)
==============================================

`Wastyk et al., Cell 2021 <https://www.cell.com/cell/fulltext/S0092-8674(21)00754-6>`_ compared the gut microbiome
of healthy adults before and after a high-fiber diet; CAZyme abundance increased with fiber. Here: two stool
metagenomes (Illumina 2x146 bp) of one participant, baseline (``fefifo_8022_1``) and after 10 weeks of high-fiber
diet (``fefifo_8022_7``), individual assembly route.

The steps are those of :ref:`metagenomics-protocol`; this page lists only what differs. Files of every step:
https://bcb.unl.edu/dbCAN_tutorial/dataset2-Wastyk2021/.

Data
----

Download and rename the reads to ``<sample>_1/2.fastq.gz``, the names used by the main protocol:

.. code-block:: shell

   SAMPLES="fefifo_8022_1 fefifo_8022_7"
   for s in $SAMPLES; do
     for r in 1 2; do
       wget -O ${s}_${r}.fastq.gz https://bcb.unl.edu/dbCAN_tutorial/dataset2-Wastyk2021/${s}__shotgun_${r}.fastq.gz
     done
   done

Steps
-----

With ``SAMPLES`` set as above, run :ref:`P1-P12 <mg-module1>` of the main protocol unchanged. Kraken2 (P1)
finds very little contamination, so no host reads are removed.

Plots (P13-P15):

.. code-block:: shell

   dbcan_plot heatmap_plot --samples fefifo_8022_1,fefifo_8022_7 --top 20 --show-abund \
       -i fefifo_8022_1_abund/fam_substrate_abund.out,fefifo_8022_7_abund/fam_substrate_abund.out --pdf substrate_heatmap.pdf
   for t in fam subfam EC; do
     dbcan_plot bar_plot --samples fefifo_8022_1,fefifo_8022_7 --top 20 --vertical-bar \
         -i fefifo_8022_1_abund/${t}_abund.out,fefifo_8022_7_abund/${t}_abund.out --pdf ${t}_abund.pdf
   done
   # CGC of the example: k141_2168|CGC1 of fefifo_8022_1 (P11 region k141_2168:4235-19858)
   dbcan_plot CGC_synteny_coverage_plot -i fefifo_8022_1.dbCAN --cgcid 'k141_2168|CGC1' \
       --reads-count fefifo_8022_1.cgc.depth.txt --db-dir db

Timing
------

.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - Steps
     - Time
   * - P1 / P2
     - ~10 min / ~20 min
   * - P3 assembly
     - ~84 min
   * - P4 gene prediction
     - ~24 min
   * - P5 / P6 / P7
     - ~10 min / ~15 min / ~5 h
   * - P8-P9 mapping and sorting
     - ~16 min
   * - P10-P12, P13-P15
     - ~1 min each
