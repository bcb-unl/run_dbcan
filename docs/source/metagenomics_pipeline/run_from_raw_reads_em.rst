Metagenomics protocol: Emilson 2024 (soil depths, metagenome + metatranscriptome)
=================================================================================

`Emilson et al., Nature Communications 2024 <https://www.nature.com/articles/s41467-023-44431-4>`_ studied
carbohydrate-related genes of soil microbiomes at different depths and hillslope positions. Here: one site at two
depths (``HS-05``, ``HS-60``), each with a metagenome (``MAG_*``) and a metatranscriptome (``RNA_*``), Illumina
paired-end. The metagenome of each depth is assembled and annotated; both read sets are mapped to it, which gives
CAZyme abundance (DNA) and expression (RNA).

The steps are those of :ref:`metagenomics-protocol`; this page lists only what differs. Files of every step:
https://bcb.unl.edu/dbCAN_tutorial/dataset6-Emilson2024/.

Data
----

DNA: `SRR20140034 <https://www.ncbi.nlm.nih.gov/sra/?term=SRR20140034>`_,
`SRR20140038 <https://www.ncbi.nlm.nih.gov/sra/?term=SRR20140038>`_; RNA:
`SRR20140125 <https://www.ncbi.nlm.nih.gov/sra/?term=SRR20140125>`_,
`SRR20140128 <https://www.ncbi.nlm.nih.gov/sra/?term=SRR20140128>`_ (renamed after the sample).

.. code-block:: shell

   for s in MAG_HS-05 MAG_HS-60 RNA_HS-05 RNA_HS-60; do
     for r in 1 2; do
       wget https://bcb.unl.edu/dbCAN_tutorial/dataset6-Emilson2024/${s}_${r}.fastq.gz
     done
   done

Steps
-----

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Steps
     - Run with
   * - :ref:`P1-P2 <mg-p1>`
     - ``SAMPLES="MAG_HS-05 MAG_HS-60 RNA_HS-05 RNA_HS-60"`` (Kraken2 finds very little contamination)
   * - :ref:`P3-P7 <mg-p3>`
     - ``SAMPLES="MAG_HS-05 MAG_HS-60"``: only the metagenomes are assembled and annotated
   * - P8-P12
     - every read set against the contigs and ``cgc.gff`` of the metagenome of the same depth (below)

.. code-block:: shell

   for d in HS-05 HS-60; do
     bwa index megahit_MAG_$d/MAG_$d.contigs.fa
     for t in MAG RNA; do
       s=${t}_$d
       bwa mem -t 32 megahit_MAG_$d/MAG_$d.contigs.fa ${s}_1_val_1.fq.gz ${s}_2_val_2.fq.gz | samtools sort -@ 32 -o $s.bam
       samtools index $s.bam
       dbcan_utils cal_coverage -g MAG_$d.dbCAN/cgc.gff -i $s.bam -o $s.depth.txt -t 6
       for f in fam_abund fam_substrate_abund CGC_abund CGC_substrate_abund; do
         dbcan_utils $f -bt $s.depth.txt -i MAG_$d.dbCAN -a TPM -o ${s}_abund
       done
     done
   done

.. note::
   For a stranded RNA-seq library add ``--strandedness reverse`` (dUTP/TruSeq stranded) or ``forward`` to
   ``cal_coverage`` for the ``RNA_*`` samples; see :ref:`dbcan-utils-abundance`.

Plots (P13-P15). ``--top`` picks the rows that differ most between the **first two** files:

.. code-block:: shell

   S=MAG_HS-05,MAG_HS-60,RNA_HS-05,RNA_HS-60
   dbcan_plot heatmap_plot --samples $S --top 20 --show-abund --pdf substrate_heatmap.pdf \
       -i MAG_HS-05_abund/fam_substrate_abund.out,MAG_HS-60_abund/fam_substrate_abund.out,RNA_HS-05_abund/fam_substrate_abund.out,RNA_HS-60_abund/fam_substrate_abund.out
   for t in fam subfam EC; do
     dbcan_plot bar_plot --samples $S --top 20 --vertical-bar --pdf ${t}_abund.pdf \
         -i MAG_HS-05_abund/${t}_abund.out,MAG_HS-60_abund/${t}_abund.out,RNA_HS-05_abund/${t}_abund.out,RNA_HS-60_abund/${t}_abund.out
   done
   # CGC of the example: k141_282901|CGC1 of MAG_HS-05 (P11 region k141_282901:2644-8574);
   # one coverage panel per read set
   samtools depth -a -r k141_282901:2644-8574 MAG_HS-05.bam > MAG_HS-05.cgc.depth.txt
   samtools depth -a -r k141_282901:2644-8574 RNA_HS-05.bam > RNA_HS-05.cgc.depth.txt
   dbcan_plot CGC_synteny_coverage_plot -i MAG_HS-05.dbCAN --cgcid 'k141_282901|CGC1' --db-dir db \
       --reads-count MAG_HS-05.cgc.depth.txt --reads-count RNA_HS-05.cgc.depth.txt --samples DNA,RNA

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
     - ~1 h
   * - P5 / P6 / P7
     - ~10 min / ~15 min / ~5 h
   * - P8-P9 mapping and sorting (4 read sets)
     - ~20 min
   * - P10-P12, P13-P15
     - a few minutes
