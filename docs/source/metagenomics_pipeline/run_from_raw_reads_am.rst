Metagenomics protocol: Amelia 2024 (soil after wildfire)
========================================================

`Amelia et al., ISME Communications 2024 <https://academic.oup.com/ismecommun/article/4/1/ycae108/7742758>`_
studied the effect of wildfires on the soil microbiome. Here: two soil metagenomes (Illumina paired-end),
``HTOOS`` and ``HCZGU`` (originally `SRR24887509 <https://www.ncbi.nlm.nih.gov/sra/?term=SRR24887509>`_ and
`SRR24887498 <https://www.ncbi.nlm.nih.gov/sra/?term=SRR24887498>`_), individual assembly route.

The steps are those of :ref:`metagenomics-protocol`; this page lists only what differs. Files of every step:
https://bcb.unl.edu/dbCAN_tutorial/dataset5-Amelia2024/.

Data
----

.. code-block:: shell

   SAMPLES="HTOOS HCZGU"
   for s in $SAMPLES; do
     for r in 1 2; do
       wget https://bcb.unl.edu/dbCAN_tutorial/dataset5-Amelia2024/${s}_${r}.fastq.gz
     done
   done

Steps
-----

With ``SAMPLES`` set as above, run :ref:`P1-P12 <mg-module1>` of the main protocol unchanged. Kraken2 (P1)
finds very little contamination, so no host reads are removed.

Plots (P13-P15):

.. code-block:: shell

   dbcan_plot heatmap_plot --samples HCZGU,HTOOS --top 20 --show-abund \
       -i HCZGU_abund/fam_substrate_abund.out,HTOOS_abund/fam_substrate_abund.out --pdf substrate_heatmap.pdf
   for t in fam subfam EC; do
     dbcan_plot bar_plot --samples HCZGU,HTOOS --top 20 --vertical-bar \
         -i HCZGU_abund/${t}_abund.out,HTOOS_abund/${t}_abund.out --pdf ${t}_abund.pdf
   done
   # CGC of the example: k141_1169008|CGC1 of HCZGU (P11 region k141_1169008:2-5684)
   dbcan_plot CGC_synteny_coverage_plot -i HCZGU.dbCAN --cgcid 'k141_1169008|CGC1' \
       --reads-count HCZGU.cgc.depth.txt --db-dir db

Timing
------

.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - Steps
     - Time
   * - Download
     - ~20 min
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
