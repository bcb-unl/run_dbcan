Supplementary protocol: assembly-free (Carter 2023)
===================================================

Reads are searched directly against the CAZy sequences of the dbCAN database with DIAMOND blastx, without
assembly (route 11 of the :ref:`routes figure <metagenomics-protocol>`). This is fast and needs little RAM, but
gives family, subfamily, EC and substrate abundance only: no genes, no CGCs, no CGC substrates.

The setup and P1-P2 are those of :ref:`metagenomics-protocol` (same data, ``SAMPLES="Wet2014 Dry2014"``). Add
``seqtk`` and ``diamond`` to the conda environment.

**S2b. Download the assembly-free mapping tables** (~1 min). ``run_dbcan database`` gives ``db/CAZy.dmnd`` and
``db/fam-substrate-mapping.tsv``; the subfamily, EC and substrate steps of P5 also need two tables that it does not
download (the nf-core/dbcan pipeline takes them from the same place):

.. code-block:: shell

   for f in CAZyID_subfam_mapping.tsv subfam_EC_mapping.tsv; do
     wget -O db/$f https://dbcan.s3.us-west-2.amazonaws.com/dbcan_asmfree/$f
   done

.. list-table::
   :header-rows: 1
   :widths: 36 34 30

   * - File in ``db/``
     - Maps
     - Used by
   * - ``CAZyID_subfam_mapping.tsv``
     - CAZy sequence → eCAMI subfamily (``GH5_e12``)
     - ``diamond_subfam_abund``
   * - ``subfam_EC_mapping.tsv``
     - subfamily → EC, subfamily → substrate
     - ``diamond_EC_abund``, ``diamond_substrate_abund``
   * - ``fam-substrate-mapping.tsv``
     - EC → substrate
     - ``diamond_substrate_abund``

**P3. Convert the trimmed reads to FASTA** (~9 min).

.. code-block:: shell

   for s in $SAMPLES; do
     seqtk seq -a ${s}_1_val_1.fq.gz > ${s}_1.fna
     seqtk seq -a ${s}_2_val_2.fq.gz > ${s}_2.fna
   done

**P4. Search the reads against CAZy** (~5 h).

.. code-block:: shell

   for s in $SAMPLES; do
     for r in 1 2; do
       diamond blastx --db db/CAZy.dmnd --query ${s}_$r.fna --out ${s}_$r.blastx --threads 32 \
           --id 80 --query-cover 90 --max-target-seqs 1 --quiet \
           --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp qlen slen
     done
   done

Keep this ``--outfmt 6`` column list: ``dbcan_asmfree`` reads it. ``--id 80`` and ``--query-cover 90`` strongly
change the results; ``--max-target-seqs 1`` assigns each read to its best CAZy sequence only.

**P5. Abundance of families, subfamilies, ECs and substrates** (~11 min).

.. code-block:: shell

   for s in $SAMPLES; do
     dbcan_asmfree diamond_fam_abund -paf1 ${s}_1.blastx -paf2 ${s}_2.blastx --raw_reads ${s}_1_val_1.fq.gz -n FPKM -o ${s}_fam_abund
     dbcan_asmfree diamond_subfam_abund -paf1 ${s}_1.blastx -paf2 ${s}_2.blastx --raw_reads ${s}_1_val_1.fq.gz -n FPKM -o ${s}_subfam_abund
     dbcan_asmfree diamond_EC_abund -i ${s}_subfam_abund -o ${s}_EC_abund
     dbcan_asmfree diamond_substrate_abund -i ${s}_subfam_abund -o ${s}_substrate_abund
   done

``-o`` is an output **file** here. ``-n`` is ``TPM`` (default), ``FPKM`` or ``RPM``; ``-d`` the database folder
(default ``./db``). Reads are grouped by their best CAZy sequence, then by family and eCAMI subfamily. EC
abundance sums the subfamilies of each EC (``subfam_EC_mapping.tsv``); substrate abundance sums the subfamilies
linked to a substrate directly (``subfam_EC_mapping.tsv``) or through their EC (``fam-substrate-mapping.tsv``).

.. note::
   ``--raw_reads`` only gives the total read number (the normalisation denominator); it may be FASTQ or FASTA
   (``.fq``/``.fastq``/``.fa``/``.fasta``/``.fna``, optionally ``.gz``). The families of a hit are the family or
   subfamily names (``GH5``, ``GH5_7``, ``CBM1``) in the CAZy sequence header, counted once each. Reads are counted
   per reference sequence: the accession of NCBI-style headers (``CAD7272370.1|GH33``), ``<protein number>|<file>``
   of JGI-style headers (``CE16|390012|Krezon1_GeneCatalog_proteins_20220127.aa.fasta``). Subfamily, EC and
   substrate abundance count only the CAZy sequences whose full header is listed in ``CAZyID_subfam_mapping.tsv``.

**P6-P7. Plots** (1 min each).

.. code-block:: shell

   for t in fam subfam EC; do
     dbcan_plot bar_plot --samples Wet2014,Dry2014 --top 20 --vertical-bar \
         -i Wet2014_${t}_abund,Dry2014_${t}_abund --pdf ${t}_abund.pdf
   done
   dbcan_plot heatmap_plot --samples Wet2014,Dry2014 --top 20 --show-abund \
       -i Wet2014_substrate_abund,Dry2014_substrate_abund --pdf substrate_heatmap.pdf
