Metagenomics protocol: Priest 2023 (HiFi metagenomes + metatranscriptomes)
==========================================================================

`Priest et al., ISME Communications 2023 <https://www.nature.com/articles/s43705-023-00324-7>`_ studied glycan
degradation in Arctic seawater with PacBio HiFi metagenomes (MG) and Illumina metatranscriptomes (MT, 2x151 bp).
Here: one site at two depths, surface water (``SRF``) and bottom of the surface mixed layer (``BML``), each with an
MG and an MT sample (``FRAM_WSC20_S25_*``). The HiFi reads are assembled with Flye; MG and MT reads are mapped to
the contigs of the same depth (abundance and expression).

The steps are those of :ref:`metagenomics-protocol`; this page gives every command that differs (long reads, host
removal, Flye, minimap2). Files of every step: https://bcb.unl.edu/dbCAN_tutorial/dataset3-Priest2023/.

Data and extra tools
--------------------

Add ``flye`` (2.9.5) and ``minimap2`` to the conda environment of :ref:`mg-setup`.

.. code-block:: shell

   URL=https://bcb.unl.edu/dbCAN_tutorial/dataset3-Priest2023
   for d in BML SRF; do
     wget $URL/${d}_MG.fastq.gz $URL/${d}_MT_1.fastq.gz $URL/${d}_MT_2.fastq.gz
   done

Module 1: reads to genes
------------------------

**P1. Contamination check** (~10 min). HiFi reads are single files (no ``--paired``):

.. code-block:: shell

   for d in BML SRF; do
     kraken2 --threads 32 --quick --db K2 --report ${d}_MG.kreport --output ${d}_MG.kraken.output ${d}_MG.fastq.gz
     kraken2 --threads 32 --quick --paired --db K2 --report ${d}_MT.kreport --output ${d}_MT.kraken.output \
         ${d}_MT_1.fastq.gz ${d}_MT_2.fastq.gz
   done

**P1b. Remove human reads** (~40 min). Kraken2 finds many human reads in this dataset. Keep the reads that do not
align to the human genome (minimap2 for HiFi, BWA for Illumina):

.. code-block:: shell

   wget https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
   bwa index -p hg38 Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
   for d in BML SRF; do
     minimap2 -a -x map-hifi -t 32 Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz ${d}_MG.fastq.gz \
         | samtools fastq -f 4 -@ 32 -0 ${d}_MG.clean.fq.gz -
     bwa mem -t 32 hg38 ${d}_MT_1.fastq.gz ${d}_MT_2.fastq.gz \
         | samtools fastq -f 12 -@ 32 -1 ${d}_MT_1.clean.fq.gz -2 ${d}_MT_2.clean.fq.gz -
   done

The KrakenTools alternative is shown in the :ref:`host-removal note <mg-host-removal>` (for HiFi reads use ``-s``
instead of ``-s1``/``-s2``).

**P2. Trim the Illumina MT reads** (~20 min). HiFi reads are not trimmed.

.. code-block:: shell

   for d in BML SRF; do
     trim_galore --paired --illumina -j 8 ${d}_MT_1.clean.fq.gz ${d}_MT_2.clean.fq.gz
   done

Output: ``BML_MT_1.clean_val_1.fq.gz`` and ``BML_MT_2.clean_val_2.fq.gz``.

**P3. Assemble the HiFi reads with Flye** (~4 h 20 min).

.. code-block:: shell

   for d in BML SRF; do
     flye --threads 32 --meta --pacbio-hifi ${d}_MG.clean.fq.gz --hifi-error 0.01 --keep-haplotypes --out-dir flye_${d}_MG
   done

Output: ``flye_BML_MG/assembly.fasta``. ``--meta`` assembles a metagenome; ``--hifi-error 0.01`` is the usual HiFi
error rate.

**P4. Predict genes** (~14 min), as in the main protocol but on the Flye contigs:

.. code-block:: shell

   for d in BML SRF; do
     pyrodigal -p meta -j 36 -i flye_${d}_MG/assembly.fasta -a ${d}_MG.faa -d ${d}_MG.cds -f gff -o ${d}_MG.gff
   done

Module 2: run_dbcan
-------------------

:ref:`P5-P7 <mg-module2>` unchanged with ``SAMPLES="BML_MG SRF_MG"`` (output ``BML_MG.dbCAN``, ``SRF_MG.dbCAN``).

Module 3: abundance
-------------------

**P8-P10.** HiFi reads are mapped with minimap2 and counted with ``--hifi`` (the overlap ratio is then relative to
the shorter of read and gene); MT reads as in the main protocol (~35 min in total):

.. code-block:: shell

   for d in BML SRF; do
     ref=flye_${d}_MG/assembly.fasta
     minimap2 -a -x map-hifi -t 32 $ref ${d}_MG.clean.fq.gz | samtools sort -@ 32 -o ${d}_MG.bam
     bwa index $ref
     bwa mem -t 32 $ref ${d}_MT_1.clean_val_1.fq.gz ${d}_MT_2.clean_val_2.fq.gz | samtools sort -@ 32 -o ${d}_MT.bam
     samtools index ${d}_MG.bam && samtools index ${d}_MT.bam
     dbcan_utils cal_coverage -g ${d}_MG.dbCAN/cgc.gff -i ${d}_MG.bam -o ${d}_MG.depth.txt -t 6 --hifi
     dbcan_utils cal_coverage -g ${d}_MG.dbCAN/cgc.gff -i ${d}_MT.bam -o ${d}_MT.depth.txt -t 6
   done

.. note::
   minimap2 writes the ``NM`` tag that ``cal_coverage`` uses for ``--identity``; ``--MD`` is not needed. With
   ``--identity 0.98`` (default) HiFi reads with more errors are not counted.

**P11-P12.** Both read sets of a depth use the annotation of that depth's MG:

.. code-block:: shell

   for d in BML SRF; do
     for t in MG MT; do
       for f in fam_abund fam_substrate_abund CGC_abund CGC_substrate_abund; do
         dbcan_utils $f -bt ${d}_$t.depth.txt -i ${d}_MG.dbCAN -a TPM -o ${d}_${t}_abund
       done
     done
   done
   # depth over the example CGC contig_12725|CGC1 of SRF_MG
   for t in MG MT; do samtools depth -a -r contig_12725:10309-25115 SRF_$t.bam > SRF_$t.cgc.depth.txt; done

Module 4: plots
---------------

.. code-block:: shell

   S=BML_MG,BML_MT,SRF_MG,SRF_MT
   dbcan_plot heatmap_plot --samples $S --top 20 --show-abund --pdf substrate_heatmap.pdf \
       -i BML_MG_abund/fam_substrate_abund.out,BML_MT_abund/fam_substrate_abund.out,SRF_MG_abund/fam_substrate_abund.out,SRF_MT_abund/fam_substrate_abund.out
   for t in fam subfam EC; do
     dbcan_plot bar_plot --samples $S --top 20 --vertical-bar --pdf ${t}_abund.pdf \
         -i BML_MG_abund/${t}_abund.out,BML_MT_abund/${t}_abund.out,SRF_MG_abund/${t}_abund.out,SRF_MT_abund/${t}_abund.out
   done
   dbcan_plot CGC_synteny_coverage_plot -i SRF_MG.dbCAN --cgcid 'contig_12725|CGC1' --db-dir db \
       --reads-count SRF_MG.cgc.depth.txt --reads-count SRF_MT.cgc.depth.txt --samples MG,MT

``--top`` picks the rows that differ most between the **first two** files (here BML MG and MT).

Timing
------

.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - Steps
     - Time
   * - P1 / P1b / P2
     - ~10 min / ~40 min / ~20 min
   * - P3 Flye assembly
     - ~4 h 20 min
   * - P4 gene prediction
     - ~14 min
   * - P5 / P6 / P7
     - ~10 min / ~15 min / ~5 h
   * - P8-P10 mapping and counting
     - ~35 min (HiFi mapping ~20 min)
   * - P11-P12, plots
     - a few minutes
