Supplementary protocol: co-assembly (Carter 2023)
=================================================

The reads of all samples are assembled together; each sample's reads are then mapped to the shared contigs. All
samples share one annotation, so contig and CGC IDs **are** comparable between samples (route 6-10 of the
:ref:`routes figure <metagenomics-protocol>`). Co-assembly needs more RAM than individual assembly.

The steps are those of :ref:`metagenomics-protocol` (same data, ``SAMPLES="Wet2014 Dry2014"``); this page lists
only what differs.

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Steps
     - Change
   * - :ref:`S1-S3, P1-P2 <mg-setup>`
     - none
   * - P3-P4
     - one assembly of the pooled trimmed reads, named ``coasm`` (below)
   * - :ref:`P5-P7 <mg-module2>`
     - once, on ``coasm.faa`` / ``coasm.gff`` → ``coasm.dbCAN``
   * - P8-P12
     - every sample against the ``coasm`` contigs and ``coasm.dbCAN`` (below)

.. code-block:: shell

   # P3-P4: co-assembly and genes
   zcat Wet2014_1_val_1.fq.gz Dry2014_1_val_1.fq.gz | gzip > coasm_1.fq.gz
   zcat Wet2014_2_val_2.fq.gz Dry2014_2_val_2.fq.gz | gzip > coasm_2.fq.gz
   megahit -m 0.5 -t 32 --min-contig-len 1000 -o megahit_coasm --out-prefix coasm -1 coasm_1.fq.gz -2 coasm_2.fq.gz
   pyrodigal -p meta -j 36 -i megahit_coasm/coasm.contigs.fa -a coasm.faa -d coasm.cds -f gff -o coasm.gff
   # P5-P7
   run_dbcan easy_substrate --mode protein --input_raw_data coasm.faa --input_gff coasm.gff --gff_type prodigal \
       --db_dir db --output_dir coasm.dbCAN --threads 32
   # P8-P12: per sample, on the shared contigs
   bwa index megahit_coasm/coasm.contigs.fa
   for s in Wet2014 Dry2014; do
     bwa mem -t 32 megahit_coasm/coasm.contigs.fa ${s}_1_val_1.fq.gz ${s}_2_val_2.fq.gz | samtools sort -@ 32 -o $s.bam
     samtools index $s.bam
     dbcan_utils cal_coverage -g coasm.dbCAN/cgc.gff -i $s.bam -o $s.depth.txt -t 6
     for f in fam_abund fam_substrate_abund CGC_abund CGC_substrate_abund; do
       dbcan_utils $f -bt $s.depth.txt -i coasm.dbCAN -a TPM -o ${s}_abund
     done
   done

Plots: :ref:`P13-P14 <mg-module4>` unchanged (``Wet2014_abund``, ``Dry2014_abund``). For one CGC (P11, P15), both
samples can be drawn on the same CGC:

.. code-block:: shell

   CONTIG=k141_597014; CGC=CGC1     # choose a CGC from coasm.dbCAN/substrate_prediction.tsv
   REGION=$(awk -F'\t' -v c="$CONTIG" -v n="$CGC" '$1==n && $2==c {print c":"$3"-"$4}' \
            coasm.dbCAN/cgc_standard_out_summary.tsv)
   for s in Wet2014 Dry2014; do samtools depth -a -r "$REGION" $s.bam > $s.cgc.depth.txt; done
   dbcan_plot CGC_synteny_coverage_plot -i coasm.dbCAN --cgcid "$CONTIG|$CGC" --db-dir db \
       --reads-count Wet2014.cgc.depth.txt --reads-count Dry2014.cgc.depth.txt --samples Wet2014,Dry2014

Timing: as the main protocol, except that the single assembly takes longer than one individual assembly and
P5-P7 run once.
