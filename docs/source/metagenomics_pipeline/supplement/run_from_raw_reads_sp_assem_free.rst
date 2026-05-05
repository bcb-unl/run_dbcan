Run from Raw Reads(Cater 2023): Supplementary Protocol for assembly-free
===============================================================================

Introduction
------------

Overview
````````

There are three sections in this document: 

Procedure of `assembly-free`


.. _cater_2023: https://www.sciencedirect.com/science/article/pii/S0092867423005974

Example : Carter 2023 Dataset ``Carter 2023``
----------------------------------------------

S1. Download Carter2023 raw reads (~10min)
``````````````````````````````````````````
To download the required raw reads, use the following wget commands:

.. code-block:: shell

    wget https://bcb.unl.edu/dbCAN_tutorial/dataset1-Carter2023/individual_assembly/Dry2014_1.fastq.gz
    wget https://bcb.unl.edu/dbCAN_tutorial/dataset1-Carter2023/individual_assembly/Dry2014_2.fastq.gz
    wget https://bcb.unl.edu/dbCAN_tutorial/dataset1-Carter2023/individual_assembly/Wet2014_1.fastq.gz
    wget https://bcb.unl.edu/dbCAN_tutorial/dataset1-Carter2023/individual_assembly/Wet2014_2.fastq.gz

These raw data were originally downloaded from
https://www.ncbi.nlm.nih.gov/sra/?term=ERR7745896
and https://www.ncbi.nlm.nih.gov/sra/?term=ERR7738162
and renamed to indicate their collected seasons.

S2. Install Anaconda (~3min)
````````````````````````````

Download and install the latest version of Anaconda for Linux from
https://www.anaconda.com/download#downloads. Once Anaconda is
successfully installed, proceed to create a dedicated conda environment
named ``CAZyme_annotation`` and activate it.
Subsequently, all the required tools can be seamlessly installed within
this environment (yml file is provided on github).

.. code-block:: shell

    conda env create -f CAZyme_annotation_main.yml
    conda activate CAZyme_annotation_main



S4. Configure databases required by run_dbcan (~2h)
```````````````````````````````````````````````````
To install the databases, execute the following commands:

.. include:: prepare_the_database.rst


Download database required by Kraken2 (very slow; can be skipped
if users do not intend to run Kraken2):

.. code-block:: shell

    kraken2-build --standard --db K2

or could download database directly from the https://benlangmead.github.io/aws-indexes/k2,
and extract the tarball to the db folder, here we use standard database as an example:

.. code-block:: shell

    wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20240605.tar.gz
    mkdir K2
    tar -xvf k2_standard_20240605.tar.gz -C K2


**CRITICAL STEP**

    The downloaded files must be all in the right location (the db folder).

    The CAZy.dmnd file is needed for DIAMOND search.

    The dbCAN.hmm and dbCAN-sub.hmm files are for HMMER search.

    The tcdb.dmnd, TF.HMM, and STP.hmm files are for CGC prediction.

    The PUL.faa file consists of protein sequences from experimentally
    validated PULs for BLAST search to predict substrates for CGCs.

    The dbCAN-PUL_12-12-2023.txt and dbCAN-PUL_12-12-2023.xlsx files contain
    PUL-substrate mapping curated from literature.

    Lastly, the
    fam-substrate-mapping.tsv file is the family-EC-substrate
    mapping table for the prediction of CAZyme substrates.

.. warning::

    Users should use a clean version of Anaconda. If the above steps failed, we suggest users reinstall their Anaconda.
    The Anaconda installation and configuration step may experience
    prolonged time while resolving environment dependencies.
    Users should be patient during this process. Alternatively,
    users may consider "mamba", another Python package manager
    that offers similar functionality to Anaconda. Information and
    access to mamba software can be found at
    https://github.com/mamba-org/mamba.

Procedure
--------------------------------------------

Module 1: Reads processing (Fig. 2) to obtain contigs
`````````````````````````````````````````````````````

P1. Contamination Check
^^^^^^^^^^^^^^^^^^^^^^^

Use `kraken2` to check for contaminated reads:

.. code-block:: shell

    kraken2 --threads 32 --quick --paired --db K2 --report Wet2014.kreport --output Wet2014.kraken.output Wet2014_1.fastq.gz Wet2014_2.fastq.gz
    kraken2 --threads 32 --quick --paired --db K2 --report Dry2014.kreport --output Dry2014.kraken.output Dry2014_1.fastq.gz Dry2014_2.fastq.gz

Kraken2 found very little contamination in the Carter2023 data. Consequently, there was no need for the contamination removal step.
For Kraken2 visualization, users can use pavian (https://github.com/fbreitwieser/pavian?tab=readme-ov-file).

If contamination is identified, users can align the reads to the reference genomes of potential contamination source organisms to remove
the aligned reads (Box 1). The most common source in human microbiome studies is from human hosts.


Box 1: Removing Contamination Reads from Humans
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Kraken2 will produce the following output files.

.. code-block:: shell

    -rw-rw-r-- 1 jinfang jinfang 2.0G Dec 12 10:24 Dry2014.kraken.output
    -rw-rw-r-- 1 jinfang jinfang 1.2M Dec 12 10:25 Dry2014.kreport
    -rw-rw-r-- 1 jinfang jinfang 5.1G Dec 12 09:47 Wet2014.kraken.output
    -rw-rw-r-- 1 jinfang jinfang 1.1M Dec 12 09:48 Wet2014.kreport


Suppose from these files, we have identified humans as the contamination source, we can use the following commands to remove the contamination reads by aligning reads to the human reference genome.

.. code-block:: shell

    wget https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
    bwa index -p hg38 Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
    bwa mem hg38 Wet2014_1.fastq.gz Wet2014_2.fastq.gz -t 32 -o Wet2014.hg38.sam
    bwa mem hg38 Dry2014_1.fastq.gz Dry2014_2.fastq.gz -t 32 -o Dry2014.hg38.sam
    samtools view -f 12 Wet2014.hg38.sam > Wet2014.hg38.unmap.bam
    samtools view -f 12 Dry2014.hg38.sam > Dry2014.hg38.unmap.bam
    samtools fastq -1 Wet2014_1.clean.fq.gz -2 Wet2014_2.clean.fq.gz Wet2014.hg38.unmap.bam
    samtools fastq -1 Dry2014_1.clean.fq.gz -2 Dry2014_2.clean.fq.gz Dry2014.hg38.unmap.bam


KrakenTools could also extract host reads quickly and easied which is recommended. We use tax 2759 (plant) as an example.
Please read KrakenTools README for more information (https://github.com/jenniferlu717/KrakenTools?tab=readme-ov-file).

.. code-block:: shell

    extract_kraken_reads.py \
    -k Dry2014.kraken.output \
    -s1 Dry2014_1.fastq.gz -s2 Dry2014_2.fastq.gz \
    --fastq-output --exclude \
    --taxid 9606 \
    -o Dry2014_1.clean.fq.gz -o2 Dry2014_2.clean.fq.gz




P2. Trim adapter and low-quality reads (TIMING ~20min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    trim_galore --paired Wet2014_1.fastq.gz Wet2014_2.fastq.gz --illumina -j 36
    trim_galore --paired Dry2014_1.fastq.gz Dry2014_2.fastq.gz --illumina -j 36

We specified `--illumina` to indicate that the reads were generated using the Illumina sequencing platform.
Nonetheless, trim_galore can automatically detect adapters, providing flexibility for users who may know the specific sequencing platform.
Details of trimming are available in the trimming report file (Box 2).

Box 2: Example output of `trim_galore`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    In addition to the trimmed read files, `Trim_galore`` also generates a trimming report file.
    The trimming report contains details on read trimming, such as the number of trimmed reads.

.. code-block:: shell

    -rw-rw-r-- 1 jinfang jinfang 4.2K Dec 13 01:48 Dry2014_1.fastq.gz_trimming_report.txt
    -rw-rw-r-- 1 jinfang jinfang 2.0G Dec 13 01:55 Dry2014_1_val_1.fq.gz
    -rw-rw-r-- 1 jinfang jinfang 4.4K Dec 13 01:55 Dry2014_2.fastq.gz_trimming_report.txt
    -rw-rw-r-- 1 jinfang jinfang 2.4G Dec 13 01:55 Dry2014_2_val_2.fq.gz
    -rw-rw-r-- 1 jinfang jinfang 4.4K Dec 13 01:30 Wet2014_1.fastq.gz_trimming_report.txt
    -rw-rw-r-- 1 jinfang jinfang 3.4G Dec 13 01:46 Wet2014_1_val_1.fq.gz
    -rw-rw-r-- 1 jinfang jinfang 4.6K Dec 13 01:46 Wet2014_2.fastq.gz_trimming_report.txt
    -rw-rw-r-- 1 jinfang jinfang 3.7G Dec 13 01:46 Wet2014_2_val_2.fq.gz

.. warning::

    During the trimming process, certain reads may be entirely removed due to low quality in its entirety.
    Using the ``--retain_unpaired`` parameter in ``trim_galore`` allows for the preservation of single-end reads.
    In this protocol, this option was not selected, so that both reads of a forward-revise pair were removed.


P3. Covert the fastq file to fasta (~9min) 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Need to install `seqtk`` first


.. code-block:: shell

    seqtk seq -a Dry2014_1_val_1.fq.gz > Dry2014_1.fna 

    seqtk seq -a Dry2014_2_val_2.fq.gz > Dry2014_2.fna 

    seqtk seq -a Wet2014_1_val_1.fq.gz > Wet2014_1.fna 

    seqtk seq -a Wet2014_2_val_2.fq.gz > Wet2014_2.fna 



Module 2. Read mapping (Fig. 2) to calculate abundance for CAZyme families, subfamilies, CGCs, and substrates
``````````````````````````````````````````````````````````````````````````````````````````````````````````````
P4. Map reads to proteins in CAZyDB by using DIAMOND blastx (~5h10min) 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    ln -s db/CAZy.dmnd . 

    diamond blastx --db CAZy --out Dry2014_1.blastx --threads 32 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp qlen slen --id 80 --query-cover 90 --query Dry2014_1.fna --max-target-seqs 1 --quiet 

    diamond blastx --db CAZy --out Dry2014_2.blastx --threads 32 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp qlen slen --id 80 --query-cover 90 --query Dry2014_2.fna --max-target-seqs 1 --quiet 

    diamond blastx --db CAZy --out Wet2014_1.blastx --threads 32 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp qlen slen --id 80 --query-cover 90 --query Wet2014_1.fna --max-target-seqs 1 --quiet 

    diamond blastx --db CAZy --out Wet2014_2.blastx --threads 32 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp qlen slen --id 80 --query-cover 90 --query Wet2014_2.fna --max-target-seqs 1 --quiet 

The parameters `--query-cover 90`` and `--id 80`` are required to align the reads to CAZy database. The thresholds of these parameters will significantly impact the results. The `--max-target-seqs 1`` parameter specify to only keep the best CAZy hit, which means each read will be mapped to only one CAZyme protein in the CAZyDB. 

P5. dbcan_asmfree to calculate the abundance of CAZyme families, subfamily, EC and substrate (~11min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We developed a Python script as dbcan_asmfree (included in the run_dbcan package) to take the DIAMOND blastx results from pair-end reads and the raw reads as inputs. With the diamond_fam_abund method, it will output the CAZyme family abundance, with the diamond_subfam_abund method, it will output the eCAMI subfamily abundance, with the diamond_EC_abund method, it will output the EC abundance and with the diamond_substrate_abund method, it will output the substrate abundance. 


.. code-block:: shell

    dbcan_asmfree diamond_fam_abund -paf1 Dry2014_1.blastx -paf2 Dry2014_2.blastx --raw_reads Dry2014_1_val_1.fq.gz  -n FPKM -o Dry2014_fam_abund 

    dbcan_asmfree diamond_fam_abund -paf1 Wet2014_1.blastx -paf2 Wet2014_2.blastx --raw_reads Wet2014_1_val_1.fq.gz -n FPKM -o Wet2014_fam_abund 

    

    dbcan_asmfree diamond_subfam_abund -paf1 Dry2014_1.blastx -paf2 Dry2014_2.blastx --raw_reads Dry2014_1_val_1.fq.gz -o Dry2014_subfam_abund -n FPKM 

    dbcan_asmfree diamond_subfam_abund -paf1 Wet2014_1.blastx -paf2 Wet2014_2.blastx --raw_reads Wet2014_1_val_1.fq.gz -o Wet2014_subfam_abund -n FPKM 

    

    dbcan_asmfree diamond_EC_abund -i Dry2014_subfam_abund -o Dry2014_EC_abund 

    dbcan_asmfree diamond_EC_abund -i Wet2014_subfam_abund -o Wet2014_EC_abund 

    

    dbcan_asmfree diamond_substrate_abund -i Dry2014_subfam_abund -o Dry2014_substrate_abund 

    dbcan_asmfree diamond_substrate_abund -i Wet2014_subfam_abund -o Wet2014_substrate_abund 

The substrate prediction in the assemble-free route is based on the fam-substrate-mapping-08012023.tsv mapping table. Specifically, from the blastx result files, reads are first grouped according to their best CAZyDB proteins. Then, reads are further grouped according to the CAZyme families and eCAMI subfamilies where the best CAZyDB proteins belong to. Lastly, the fam-substrate-mapping-08012023.tsv mapping table is used to map CAZyme families and ECs (from eCAMI subfamilies) to substrates. 
Therefore, the substrate abundance is calculated in a similar way as implemented in dbCAN-sub. 



Module 4: dbcan_plot for data visualization (Fig. 2) of abundances of CAZymes, CGCs, and substrates (TIMING variable)
`````````````````````````````````````````````````````````````````````````````````````````````````````````````````````



P6. Barplot for CAZyme family, subfamily, and EC abundance across samples (TIMING 1min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    dbcan_plot bar_plot --samples Wet2014,Dry2014 --vertical_bar --top 20 -i Wet2014_fam_abund,Dry2014_fam_abund --pdf fam.pdfDry2014.blastx.CAZy.FPKM.tsv,Wet2014.blastx.CAZy.FPKM.tsv 

    dbcan_plot bar_plot --samples Wet2014,Dry2014 --vertical_bar --top 20 -i Wet2014_subfam_abund,Dry2014_subfam_abund --pdf subfam.pdf 

    dbcan_plot bar_plot --samples Wet2014,Dry2014 --vertical_bar --top 20 -i Wet2014_fam_abund,Dry2014_fam_abund --pdf ec.pdf 

Here we plot the top 20 substrates in the two samples. The input files are the two CAZyme substrate abundance files calculated based on
dbCAN-sub result. The default heatmap is ranked by substrate abundances. To rank the heatmap according to abundance profile using
the clustermap function of the seaborn package (https://github.com/mwaskom/seaborn), users can invoke the ``--cluster_map`` parameter.

P7. Heatmap for CAZyme substrate abundance across samples (Fig. S2D) (TIMING 1min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    dbcan_plot heatmap_plot --samples Wet2014,Dry2014 --show_abund --top 20 -i  Wet2014_substrate_abund,Dry2014_substrate_abund 
