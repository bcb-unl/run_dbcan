Run from Raw Reads(Amelia 2024): Automated CAZyme and Glycan Substrate Annotation in Microbiomes: A Step-by-Step Protocol
=========================================================================================================================



.. _amelia_2024: https://academic.oup.com/ismecommun/article/4/1/ycae108/7742758#480953371

Example 2: Amelia 2024 Dataset `amelia_2024`_
---------------------------------------------

The Amelia 2024 dataset `amelia_2024`_ was published in 2024 from a soil metagenomic study. In the published paper. They analyzed the effects of wildfires on the soil microbiome.

The raw read data, intermediate data from each analysis step, and final result data and visualization files are organized in nested folders available on our website https://bcb.unl.edu/dbCAN_tutorial/dataset5-Amelia2024/).

Procedure
---------

Software and versions
`````````````````````

- **Anaconda** (`Anaconda <https://www.anaconda.com>`_, version 23.7.3)
- **MEGAHIT** (`MEGAHIT <https://github.com/voutcn/megahit>`_, version 1.2.9)
- **pyHMMER** (`pyHMMER <https://pyhmmer.readthedocs.io/en/stable/>`_, newest version)
- **DIAMOND** (`DIAMOND <https://github.com/bbuchfink/diamond>`_, version 2.1.8)
- **TrimGalore** (`TrimGalore <https://github.com/FelixKrueger/TrimGalore>`_, version 0.6.0)
- **Pyrodigal** (`Pyrodigal <https://pyrodigal.readthedocs.io/en/stable/>`_, newest version)
- **BWA** (`BWA <https://github.com/lh3/bwa>`_, version 0.7.17)
- **Samtools** (`Samtools <https://github.com/samtools/samtools>`_, version 1.7)
- **Kraken2** (`Kraken2 <https://ccb.jhu.edu/software/kraken2/>`_, version 2.1.1)
- **KrakenTools** (`KrakenTools <https://github.com/jenniferlu717/KrakenTools?tab=readme-ov-file>`_, version 1.2)
- **run_dbcan** (`run_dbcan <https://github.com/Xinpeng021001/run_dbCAN_new>`_, version 5.0.3)
- **MMseqs2** (`MMseqs2 <https://github.com/soedinglab/MMseqs2>`_, release 15-6f452)
- **Pysam** (`Pysam <https://github.com/pysam-developers/pysam>`_, version 0.22.1)

.. code-block:: shell

    conda env create -f CAZyme_annotation_am.yml
    conda activate CAZyme_annotation_am

To install the databases, execute the following commands:

.. include:: prepare_the_database.rst



Module 1: Reads processing to obtain contigs
`````````````````````````````````````````````````````

S1| Download Amelia2024 (Table 2) raw reads (~20min)
.. code-block:: shell

    wget https://bcb.unl.edu/dbCAN_tutorial/dataset5-Amelia2024/HTOOS_1.fastq.gz
    wget https://bcb.unl.edu/dbCAN_tutorial/dataset5-Amelia2024/HTOOS_2.fastq.gz
    wget https://bcb.unl.edu/dbCAN_tutorial/dataset5-Amelia2024/HCZGU_1.fastq.gz
    wget https://bcb.unl.edu/dbCAN_tutorial/dataset5-Amelia2024/HCZGU_2.fastq.gz


These raw data were originally downloaded from https://www.ncbi.nlm.nih.gov/sra/?term=SRR24887509 and https://www.ncbi.nlm.nih.gov/sra/?term=SRR24887498 and renamed to indicate their collected samples (Table 2).

.. code-block:: shell

    kraken2-build --standard --db K2

or could download database directly from the https://benlangmead.github.io/aws-indexes/k2,
and extract the tarball to the db folder, here we use standard database as an example:

.. code-block:: shell

    wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20240605.tar.gz
    mkdir K2
    tar -xvf k2_standard_20240605.tar.gz -C K2

Use `kraken2` to check for contaminated reads:

.. code-block:: shell

    kraken2 --threads 32 --quick --paired --db K2 --report HTOOS.kreport --output HTOOS.kraken.output HTOOS_1.fastq.gz HTOOS_2.fastq.gz

    kraken2 --threads 32 --quick --paired --db K2 --report HCZGU.kreport --output HCZGU.kraken.output HCZGU_1.fastq.gz HCZGU_2.fastq.gz

Kraken2 found very little contamination in the Amelia2024 data. Consequently, there was no need for the contamination removal step.

If contamination is identified, users can align the reads to the reference genomes of potential contamination source organisms to remove
the aligned reads (Box 1). The most common source in human microbiome studies is from human hosts.

Box 1: Example to remove contamination reads from human
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Kraken2 will produce the following output files:

.. code-block:: shell

    -rw-r--r-- 1 jinfang yinlab 9.2G May 31 05:59 HTOOS.kraken.output
    -rw-r--r-- 1 jinfang yinlab 1.2M May 31 05:59 HTOOS.kreport
    -rw-r--r-- 1 jinfang yinlab 4.9G May 31 03:21 HCZGU.kraken.output
    -rw-r--r-- 1 jinfang yinlab 1.1M May 31 03:21 HCZGU.kreport

Suppose from these files, we have identified humans as the contam


Suppose from these files, we have identified humans as the contamination source, we can use the following commands to remove the contamination reads by aligning reads to the human reference genome.

.. code-block:: shell

    wget https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
    bwa index -p hg38 Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
    bwa mem hg38 HTOOS_1.fastq.gz HTOOS_2.fastq.gz -t 32 -o HTOOS.hg38.sam
    bwa mem hg38 HCZGU_1.fastq.gz HCZGU_2.fastq.gz -t 32 -o HCZGU.hg38.sam
    samtools view -f 12 HCZGU.hg38.sam > HCZGU.hg38.unmap.bam
    samtools view -f 12 HTOOS.hg38.sam > HTOOS.hg38.unmap.bam
    samtools fastq -1 HCZGU_1.clean.fq.gz -2 HCZGU_2.clean.fq.gz HCZGU.hg38.unmap.bam
    samtools fastq -1 HTOOS_1.clean.fq.gz -2 HTOOS_2.clean.fq.gz HTOOS.hg38.unmap.bam

KrakenTools could also extract host reads quickly and easied which is recommended. We use tax 2759 (plant) as an example.
Please read KrakenTools README for more information (https://github.com/jenniferlu717/KrakenTools?tab=readme-ov-file).

.. code-block:: shell

    extract_kraken_reads.py -k HTOOS.kraken.output  -taxid 9606 -exclude -s1 HTOOS_1.fq.gz -s2 HTOOS_2.fastq.gz -o HTOOS_1.clean.fq -o2 HTOOS_2.clean.fq
    gzip HTOOS_1.clean.fq
    gzip HTOOS_2.clean.fq


    extract_kraken_reads.py -k HCZGU.kraken.output -taxid 9606 -exclude -s1 HCZGU_1.fastq.gz -s2 HCZGU_2.fastq.gz -o HCZGU_1.clean.fq -o2 HCZGU_2.clean.fq
    gzip HCZGU_1.clean.fq
    gzip HCZGU_2.clean.fq


P2. Trim adapter and low-quality reads (TIMING ~20min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    trim_galore --paired HTOOS_1.fastq.gz HTOOS_2.fastq.gz --illumina -j 36

    trim_galore --paired HCZGU_1.fastq.gz HCZGU_2.fastq.gz --illumina -j 36

We specified --illumina to indicate that the reads were generated using the Illumina sequencing platform.
Nonetheless, trim_galore can automatically detect adapters, providing flexibility for users who may know the specific sequencing platform.
Details of trimming are available in the trimming report file (Box 2).

Box 2: Example output of `trim_galore`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    In addition to the trimmed read files, `Trim_galore`` also generates a trimming report file.
    The trimming report contains details on read trimming, such as the number of trimmed reads.

.. code-block:: shell

    -rw-r--r-- 1 jinfang yinlab 5.1K May 31 03:27 HCZGU_1.fastq.gz_trimming_report.txt
    -rw-r--r-- 1 jinfang yinlab 8.5G May 31 03:45 HCZGU_1_val_1.fq.gz
    -rw-r--r-- 1 jinfang yinlab 5.3K May 31 03:45 HCZGU_2.fastq.gz_trimming_report.txt
    -rw-r--r-- 1 jinfang yinlab 8.5G May 31 03:45 HCZGU_2_val_2.fq.gz
    -rw-r--r-- 1 jinfang yinlab 5.2K May 31 06:10 HTOOS_1.fastq.gz_trimming_report.txt
    -rw-r--r-- 1 jinfang yinlab  16G May 31 06:45 HTOOS_1_val_1.fq.gz
    -rw-r--r-- 1 jinfang yinlab 5.4K May 31 06:45 HTOOS_2.fastq.gz_trimming_report.txt
    -rw-r--r-- 1 jinfang yinlab  16G May 31 06:45 HTOOS_2_val_2.fq.gz

.. warning::

    During the trimming process, certain reads may be entirely removed due to low quality in its entirety.
    Using the ``--retain_unpaired`` parameter in ``trim_galore`` allows for the preservation of single-end reads.
    In this protocol, this option was not selected, so that both reads of a forward-revise pair were removed.

P3. Assemble reads into contigs (TIMING ~84min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Use Megahit for assembling reads into contigs:

.. code-block:: shell

    megahit -m 0.5 -t 32 -o megahit_HTOOS -1 HTOOS_1_val_1.fq.gz -2 HTOOS_2_val_2.fq.gz --out-prefix HTOOS --min-contig-len 1000

    megahit -m 0.5 -t 32 -o megahit_HCZGU -1 HCZGU_1_val_1.fq.gz -2 HCZGU_2_val_2.fq.gz --out-prefix HCZGU --min-contig-len 1000

MEGAHIT generates two output folders: megahit_HCZGU and megahit_HTOOS.
Each contains five files and one sub-folder (Box 3). HTOOS.contigs.fa is the final contig sequence file. We set --min-contig-len 1000, a common practice to retain all contigs longer than 1,000 base pairs.

Box 3: Example output of `MEGAHIT`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: shell

    -rw-r--r-- 1 jinfang yinlab   262 Jun  3 04:58 checkpoints.txt
    -rw-r--r-- 1 jinfang yinlab     0 Jun  3 04:58 done
    -rw-r--r-- 1 jinfang yinlab  2.1G Jun  3 04:58 HTOOS.contigs.fa
    -rw-r--r-- 1 jinfang yinlab  1.4M Jun  3 04:58 HTOOS.log
    drwxr-xr-x 2 jinfang yinlab   41K Jun  3 04:57 intermediate_contigs
    -rw-r--r-- 1 jinfang yinlab  1.1K Jun  1 18:59 options.json


P4. Predict genes by `Pyrodigal` (TIMING ~24min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    pyrodigal -p meta -i megahit_HTOOS/HTOOS.contigs.fa -d fefifo_8022_1.cds -a fefifo_8022_1.faa -f gff -o fefifo_8022_1.gff -j 36
    pyrodigal -p meta -i megahit_HCZGU/HCZGU.contigs.fa -d fefifo_8022_7.cds -a fefifo_8022_7.faa -f gff -o fefifo_8022_7.gff -j 36


Box 4: Example output of `Pyrodigal`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: shell

    -rw-r--r-- 1 jinfang yinlab 2.1G Jun  3 07:15 HTOOS.cds
    -rw-r--r-- 1 jinfang yinlab 921M Jun  3 07:15 HTOOS.faa
    -rw-r--r-- 1 jinfang yinlab 857M Jun  3 07:15 HTOOS.gff

Hint: If you apply ``prodigal`` not ``pyrodigal``, use this to fix gene ID in gff files by dbcan_utils (TIMING ~1min)

.. code-block:: shell

    dbcan_utils gff_fix -i HTOOS.faa -g HTOOS.gff
    dbcan_utils gff_fix -i HCZGU.faa -g HCZGU.gff


Module 2. run_dbcan annotation (Fig. 3) to obtain CAZymes, CGCs, and substrates
```````````````````````````````````````````````````````````````````````````````

**CRITICAL STEP**

Users can skip P5 and P6, and directly run P7 (much slower though), if they want to predict not only CAZymes and CGCs, but also substrates.

P5. CAZyme annotation at the CAZyme family level (TIMING ~10min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    run_dbcan CAZyme_annotation --input_raw_data HTOOS.faa --mode protein --output_dir HTOOS.CAZyme --db_dir db
    run_dbcan CAZyme_annotation --input_raw_data HCZGU.faa --mode protein --output_dir HCZGU.CAZyme --db_dir db

Two arguments are required for ``run_dbcan``: the input sequence file (faa files) and the sequence type (protein).
By default, ``run_dbcan`` will use three methods (``pyHMMER`` vs ``dbCAN HMMdb``, ``DIAMOND`` vs ``CAZy``, ``pyHMMER`` vs ``dbCAN-sub HMMdb``) for
CAZyme annotation.

Box 5: CAZyme annotation with default setting
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The sequence type can be `protein`, `prok`, `meta`. If the input sequence file contains metagenomic contig sequences (`fna` file),
the sequence type has to be `meta`, and `prodigal` will be called to predict genes.

P6. CGC prediction (TIMING ~15 min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following commands will re-run run_dbcan to not only predict CAZymes but also CGCs with protein `faa` and gene location `gff` files.

.. code-block:: shell

    run_dbcan easy_CGC --input_raw_data HTOOS.faa --mode protein  --output_dir HTOOS.PUL --input_format NCBI --gff_type prodigal --input_gff HTOOS.gff --db_dir db
    run_dbcan easy_CGC --input_raw_data HCZGU.faa --mode protein  --output_dir HCZGU.PUL --input_format NCBI --gff_type prodigal --input_gff HCZGU.gff --db_dir db

.. warning::

    **Creating own gff file**
    If the users would like to create their own ``gff`` file (instead of using Prokka or Prodigal),
    it is important to make sure the value of ID attribute in the ``gff`` file matches the protein ID in the protein ``faa`` file.

    **[Troubleshooting]CGC not found**
    If no result is found in CGC output file, it is most likely because the sequence IDs in ``gff`` file and ``faa`` file do not match.
    Another less likely reason is that the contigs are too short and fragmented and not suitable for CGC prediction.


P7. Substrate prediction for CAZymes and CGCs (TIMING ~5h)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following commands will re-run run_dbcan to predict CAZymes, CGCs, and their substrates.

.. code-block:: shell

    run_dbcan easy_substrate --input_raw_data HTOOS.faa --mode protein  --input_format NCBI --input_gff HTOOS.gff --gff_type prodigal   --output_dir HTOOS.dbCAN --db_dir db
    run_dbcan easy_substrate --input_raw_data HCZGU.faa --mode protein  --input_format NCBI --input_gff HCZGU.gff --gff_type prodigal   --output_dir HCZGU.dbCAN --db_dir db



In the `HCZGU.dbCAN <https://bcb.unl.edu/dbCAN_tutorial/dataset5-Amelia2024/HCZGU.dbCAN/>_` directory, a total of 17 files and 1 folder are generated.



Module 3. Read mapping (Fig. 3) to calculate abundance for CAZyme families, subfamilies, CGCs, and substrates
``````````````````````````````````````````````````````````````````````````````````````````````````````````````

P8. Read mapping to all contigs of each sample (TIMING ~10 min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    mkdir samfiles
    bwa index megahit_HCZGU/HCZGU.contigs.fa
    bwa index megahit_HTOOS/HTOOS.contigs.fa
    bwa mem -t 32 -o samfiles/HCZGU.sam megahit_HCZGU/HCZGU.contigs.fa HCZGU_1_val_1.fq.gz HCZGU_2_val_2.fq.gz
    bwa mem -t 32 -o samfiles/HTOOS.sam megahit_HTOOS/HTOOS.contigs.fa HTOOS_1_val_1.fq.gz HTOOS_2_val_2.fq.gz



P9. Sort SAM files by coordinates (TIMING ~6min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    cd samfiles
    samtools sort -@ 32 -o HCZGU.bam HCZGU.sam
    samtools sort -@ 32 -o HTOOS.bam HTOOS.sam
    samtools index HCZGU.bam
    samtools index HTOOS.bam
    rm -rf *sam
    cd ..


P10. Read count calculation for all proteins of each sample using dbcan_utils (TIMING ~1min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    mkdir HCZGU_abund && cd HCZGU_abund
    dbcan_utils cal_coverage -g ../HCZGU.gff -i  ../samfiles/HCZGU.bam -o HCZGU.depth.txt -t 6 --overlap_base_ratio 0.2 --mapping_quality 30 --identity 0.98

    cd .. && mkdir HTOOS_abund && cd HTOOS_abund
    dbcan_utils cal_coverage -g ../HTOOS.fix.gff -i ../samfiles/HTOOS.bam -o HTOOS.depth.txt -t 6 --overlap_base_ratio 0.2 --mapping_quality 30 --identity 0.98
    cd ..


P11. Read count calculation for a given region of contigs using Samtools (TIMING ~1min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    cd HCZGU_abund
    samtools depth -r k141_1169008:2-5684 ../samfiles/HCZGU.bam > HCZGU.cgc.depth.txt


The parameter -r k141_1169008:2-5684 specifies a region in a contig. For any CGC, its positional range can be found in the file

.. warning::

    The contig IDs are automatically generated by MEGAHIT. There is a small chance that a same contig ID appears in both samples. However, the two contigs in the two samples do not match each other even the ID is the same. For example, the contig ID ``k141_2168`` is most likely only found in the ``fefifo_8022_1`` sample. Even if there is a ``k141_2168`` in ``fefifo_8022_7``, the actual contigs in two samples are different.

P12. dbcan_utils to calculate the abundance of CAZyme families, subfamilies, CGCs, and substrates (TIMING ~1min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    dbcan_utils fam_abund -bt HCZGU.depth.txt -i ../HCZGU.dbCAN -a TPM
    dbcan_utils fam_substrate_abund -bt HCZGU.depth.txt -i ../HCZGU.dbCAN -a TPM
    dbcan_utils CGC_abund -bt HCZGU.depth.txt -i ../HCZGU.dbCAN -a TPM
    dbcan_utils CGC_substrate_abund -bt HCZGU.depth.txt -i ../HCZGU.dbCAN -a TPM

.. code-block:: shell

    cd .. && cd HTOOS_abund
    dbcan_utils fam_abund -bt HTOOS.depth.txt -i ../HTOOS.dbCAN -a TPM
    dbcan_utils fam_substrate_abund -bt HTOOS.depth.txt -i ../HTOOS.dbCAN -a TPM
    dbcan_utils CGC_abund -bt HTOOS.depth.txt -i ../HTOOS.dbCAN -a TPM
    dbcan_utils CGC_substrate_abund -bt HTOOS.depth.txt -i ../HTOOS.dbCAN -a TPM
    cd ..

We developed a set of Python scripts as ``dbcan_utils`` (included in the ``run_dbcan`` package) to take the raw read counts for all genes as input and output the normalized abundances (refer to Box 7) of CAZyme families, subfamilies, CGCs, and substrates (see Fig. 4). The parameter ``-a TPM`` can also be set to two other metrics: RPM, or RPKM61.

- **RPKM** is calculated as the number of mapped reads to a gene G divided by [(total number of mapped reads to all genes / 10^6) x (gene G length / 1000)].
- **RPM** is the number of mapped reads to a gene G divided by (total number of mapped reads to all genes / 10^6).
- **TPM** is calculated as [number of mapped reads to a gene G / (gene G length / 1000)] divided by the sum of [number of mapped reads to each gene / (the gene length / 1000)].


Box 7. Example output of dbcan_utils
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As an example, `HCZGU_abund  <https://bcb.unl.edu/dbCAN_tutorial/dataset5-Amelia2024/HCZGU_abund/>_` folder has 7 TSV files:

.. code-block:: shell

    -rw-rw-r--  1 jinfang jinfang 178K Jan  2 04:08 CGC_abund.out
    -rw-rw-r--  1 jinfang jinfang 3.3K Jan  2 04:08 CGC_substrate_majority_voting.out
    -rw-rw-r--  1 jinfang jinfang  12K Jan  2 04:08 CGC_substrate_PUL_homology.out
    -rw-rw-r--  1 jinfang jinfang 2.5K Jan  2 04:08 EC_abund.out
    -rw-rw-r--  1 jinfang jinfang 4.1K Jan  2 04:08 fam_abund.out
    -rw-rw-r--  1 jinfang jinfang  42K Jan  2 04:08 fam_substrate_abund.out
    -rw-rw-r--  1 jinfang jinfang  26K Jan  2 04:08 subfam_abund.out

Explanation of columns in these TSV files is as follows:

    - ``fam_abund.out``: CAZy family (from HMMER vs dbCAN HMMdb), sum of TPM, number of CAZymes in the family.
    - ``subfam_abund.out``: eCAMI subfamily (from HMMER vs dbCAN-sub HMMdb), sum of TPM, number of CAZymes in the subfamily.
    - ``EC_abund.out``: EC number (extracted from dbCAN-sub subfamily), sum of TPM, number of CAZymes with the EC.
    - ``fam_substrate_abund.out``: Substrate (from HMMER vs dbCAN-sub HMMdb), sum of TPM (all CAZymes in this substrate group), GeneID (all CAZyme IDs in this substrate group).
    - ``CGC_abund.out``: CGC_ID (e.g., k141_338400|CGC1), mean of TPM (all genes in the CGC), Seq_IDs (IDs of all genes in the CGC), TPM (of all genes in the CGC), Families (CAZyme family or other signature gene type of all genes in the CGC).
    - ``CGC_substrate_PUL_homology.out``: Substrate (from dbCAN-PUL blast search), sum of TPM, CGC_IDs (all CGCs predicted to have the substrate from dbCAN-PUL blast search), TPM (of CGCs in this substrate group).
    - ``CGC_substrate_majority_voting.out``: Substrate (from dbCAN-sub majority voting), sum of TPM, CGC_IDs (all CGCs predicted to have the substrate from dbCAN-sub majority voting), TPM (of CGCs in this substrate group).


Module 4: dbcan_plot for data visualization (Fig. 3) of abundances of CAZymes, CGCs, and substrates (TIMING variable)
`````````````````````````````````````````````````````````````````````````````````````````````````````````````````````

**CRITICAL STEP**

To visualize the CAZyme annotation result, we provide a set of Python scripts as dbcan_plot to make publication quality plots with the dbcan_utils results as the input. The dbcan_plot scripts are included in the run_dbcan package. Once the plots are made in PDF format, they can be transferred to users' Windows or Mac computers for visualization.

Five data folders will be needed as the input for ``dbcan_plot``:

1. two abundance folders ``HCZGU_abund`` and ``HTOOS_abund``,
2. two CAZyme annotation ``HCZGU.dbCAN`` and ``HTOOS.dbCAN``, and
3. the ``dbCAN-PUL folder`` (under the db folder, released from ``dbCAN-PUL.tar.gz``).

P13. Heatmap for CAZyme substrate abundance across samples (TIMING 1min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    dbcan_plot heatmap_plot --samples HCZGU,HTOOS -i HCZGU_abund/fam_substrate_abund.out,HTOOS_abund/fam_substrate_abund.out --show_abund --top 20

Here we plot the top 20 substrates in the two samples. The input files are the two CAZyme substrate abundance files calculated based on dbCAN-sub result. The default heatmap is ranked by substrate abundances. To rank the heatmap according to abundance profile using the function clustermap of seaborn package, users can invoke the ``--cluster_map`` parameter.

P14. Barplot for CAZyme family/subfamily/EC abundance across samples (Fig. S4C) (TIMING 1min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    dbcan_plot bar_plot --samples HCZGU,HTOOS --vertical_bar --top 20 -i HCZGU_abund/fam_abund.out,HTOOS_abund/fam_abund.out --pdf fam.pdf --db_dir db
    dbcan_plot bar_plot --samples HCZGU,HTOOS --vertical_bar --top 20 -i HCZGU_abund/subfam_abund.out,HTOOS_abund/subfam_abund.out --pdf subfam.pdf --db_dir db
    dbcan_plot bar_plot --samples HCZGU,HTOOS --vertical_bar --top 20 -i HCZGU_abund/EC_abund.out,HTOOS_abund/EC_abund.out --pdf ec.pdf --db_dir db


Users can choose to generate a barplot instead of heatmap using the ``bar_plot`` method.

P15. Synteny plot between a CGC and its best PUL hit with read mapping coverage to CGC (Fig. S4A) (TIMING 1min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    dbcan_plot CGC_synteny_coverage_plot -i HCZGU.dbCAN --readscount HCZGU_abund/HCZGU.cgc.depth.txt --cgcid '‘k141_1169008|CGC1'

The ``HCZGU.dbCAN`` folder contains the ``PUL_blast.out`` file. Using this file, the ``cgc_standard.out`` file,
and the best PUL's ``gff`` file in ``dbCAN-PUL.tar.gz``, the CGC_synteny_plot method will create the ``CGC-PUL synteny plot``.
The ``-cgcid`` parameter is required to specify which CGC to be plotted (``'k141_1169008|CGC1'`` in this example).
The ``HCZGU.cgc.depth.txt`` file is used to plot the read mapping coverage.

If users only want to plot the CGC structure:

.. code-block:: shell

    dbcan_plot CGC_plot -i HCZGU.dbCAN --cgcid 'k141_1169008|CGC1'

If users only want to plot the CGC structure plus the read mapping coverage:

.. code-block:: shell

    dbcan_plot CGC_coverage_plot -i HCZGU.dbCAN --cgcid 'k141_1169008|CGC1'  --readscount HCZGU_abund/HCZGU.cgc.depth.txt

If users only want to plot the synteny between the CGC and PUL:

.. code-block:: shell

    dbcan_plot CGC_synteny_plot -i HCZGU.dbCAN --cgcid 'k141_1169008|CGC1'


.. warning::

    The CGC IDs in different samples do not match each other. For example, specifying -i HCZGU.dbCAN is to plot the 'k141_1169008|CGC1' in the HCZGU sample.
    The 'k141_1169008|CGC1' in the HTOOS sample most likely does not exist, and even it does, the CGC has a different sequence even if the ID is the same.
