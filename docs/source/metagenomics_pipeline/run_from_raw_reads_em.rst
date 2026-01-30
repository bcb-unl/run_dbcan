Run from Raw Reads(Emilson 2024): Automated CAZyme and Glycan Substrate Annotation in Microbiomes: A Step-by-Step Protocol
===========================================================================================================================



.. _emilson_2024: https://www.nature.com/articles/s41467-023-44431-4

Example 2: Emilson 2024 Dataset
-------------------------------

The Emilson 2024 dataset `emilson_2024`_

This paper explores the distribution of carbohydrate-related genes by exploring soil microbiomes at different depths and hillsides.

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

    conda env create -f CAZyme_annotation_em.yml
    conda activate CAZyme_annotation_em


Module 1: Reads processing to obtain contigs
`````````````````````````````````````````````````````

P1. Contamination Check (TIMING ~10min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Download the Emilson 2024 dataset:

.. code-block:: shell

    wget https://bcb.unl.edu/dbCAN_tutorial/dataset6-Emilson2024/MAG_HS-05_1.fastq.gz 
    wget https://bcb.unl.edu/dbCAN_tutorial/dataset6-Emilson2024/MAG_HS-05_2.fastq.gz  
    wget https://bcb.unl.edu/dbCAN_tutorial/dataset6-Emilson2024/MAG_HS-60_1.fastq.gz  
    wget https://bcb.unl.edu/dbCAN_tutorial/dataset6-Emilson2024/MAG_HS-60_2.fastq.gz  
    wget https://bcb.unl.edu/dbCAN_tutorial/dataset6-Emilson2024/RNA_HS-05_1.fastq.gz 
    wget https://bcb.unl.edu/dbCAN_tutorial/dataset6-Emilson2024/RNA_HS-05_2.fastq.gz  
    wget https://bcb.unl.edu/dbCAN_tutorial/dataset6-Emilson2024/RNA_HS-60_1.fastq.gz  
    wget https://bcb.unl.edu/dbCAN_tutorial/dataset6-Emilson2024/RNA_HS-60_2.fastq.gz  

These MAG DNA-seq raw data were originally downloaded from https://www.ncbi.nlm.nih.gov/sra/?term=SRR20140034 and https://www.ncbi.nlm.nih.gov/sra/?term=SRR20140038. The RNA-seq raw data were originally downloaded from https://www.ncbi.nlm.nih.gov/sra/?term=SRR20140125 and https://www.ncbi.nlm.nih.gov/sra/?term=SRR20140128. They are renamed to indicate their collected samples

Use `kraken2` to check for contaminated reads:

.. code-block:: shell

    kraken2 --threads 32 --quick --paired --db K2 --report MAG_HS-05.kreport --output MAG_HS-05.kraken.output MAG_HS-05_1.fastq.gz MAG_HS-05_2.fastq.gz 

    kraken2 --threads 32 --quick --paired --db K2 --report MAG_HS-60.kreport --output MAG_HS-60.kraken.output MAG_HS-60_1.fastq.gz MAG_HS-60_2.fastq.gz 

    kraken2 --threads 32 --quick --paired --db K2 --report RNA_HS-05.kreport --output RNA_HS-05.kraken.output RNA_HS-05_1.fastq.gz RNA_HS-05_2.fastq.gz 

    kraken2 --threads 32 --quick --paired --db K2 --report RNA_HS-60.kreport --output RNA_HS-60.kraken.output RNA_HS-60_1.fastq.gz RNA_HS-60_2.fastq.gz 



Kraken2 found very little contamination in the Emilson2024 data. Consequently, there was no need for the contamination removal step.  

If contamination is identified, users can align the reads to the reference genomes of potential contamination source organisms to remove the aligned reads (Box 1). The most common source in human microbiome studies is from human hosts. 

Box 1: Example to remove contamination reads from human
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Kraken2 will produce the following output files:

.. code-block:: shell

    -rw-r--r-- 1 jinfang yinlab 702M May 30 00:23 MAG_HS-05.kraken.output 
    -rw-r--r-- 1 jinfang yinlab 819K May 30 00:23 MAG_HS-05.kreport 
    -rw-r--r-- 1 jinfang yinlab 1.5G May 30 00:22 MAG_HS-60.kraken.output 
    -rw-r--r-- 1 jinfang yinlab 903K May 30 00:22 MAG_HS-60.kreport
    -rw-r--r-- 1 jinfang yinlab 702M May 30 00:23 RNA_HS-05.kraken.output 
    -rw-r--r-- 1 jinfang yinlab 819K May 30 00:23 RNA_HS-05.kreport 
    -rw-r--r-- 1 jinfang yinlab 1.5G May 30 00:22 RNA_HS-60.kraken.output 
    -rw-r--r-- 1 jinfang yinlab 903K May 30 00:22 RNA_HS-60.kreport 

Suppose from these files, we have identified humans as the contamination source, we can use the following commands to remove the contamination reads by aligning reads to the human reference genome.

.. code-block:: shell

    wget https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz 

    bwa index -p hg38 Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz 

    bwa mem hg38 MAG_HS-05_1.fastq.gz MAG_HS-05_2.fastq.gz -t 32 -o MAG_HS-05.hg38.sam 

    bwa mem hg38 MAG_HS-60_1.fastq.gz MAG_HS-60_2.fastq.gz -t 32 -o MAG_HS-60.hg38.sam 

    samtools view -f 12 MAG_HS-05.hg38.sam > MAG_HS-05.hg38.unmap.bam 

    samtools view -f 12 MAG_HS-60.hg38.sam > MAG_HS-60.hg38.unmap.bam 

    samtools fastq -1 MAG_HS-05_1.clean.fq.gz -2 MAG_HS-05_2.clean.fq.gz MAG_HS-05.hg38.unmap.bam 

    samtools fastq -1 MAG_HS-60_1.clean.fq.gz -2 MAG_HS-60_2.clean.fq.gz MAG_HS-60.hg38.unmap.bam 


Or by KrakenTool

.. code-block:: shell

    extract_kraken_reads.py -k MAG_HS-05.kraken.output -taxid 9606 -exclude -s1 MAG_HS-05_1.fastq.gz -s2 MAG_HS-05_2.fastq.gz -o MAG_HS-05_1.clean.fastq -o2 MAG_HS-05_2.clean.fastq 

    gzip MAG_HS-05_1.clean.fastq 

    gzip MAG_HS-05_2.clean.fastq 

    extract_kraken_reads.py -k MAG_HS-60.kraken.output -taxid 9606 -exclude -s1 MAG_HS-60_1.fastq.gz -s2 MAG_HS-60_2.fastq.gz -o MAG_HS-60_1.clean.fastq -o2 MAG_HS-60_2.clean.fastq 
    
    gzip MAG_HS-60_1.clean.fastq 

    gzip MAG_HS-60_2.clean.fastq 


P2. Trim adapter and low-quality reads (TIMING ~20min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    trim_galore --paired MAG_HS-60_1.fastq.gz MAG_HS-60_2.fastq.gz --illumina -j 36 

    trim_galore --paired MAG_HS-05_1.fastq.gz MAG_HS-05_2.fastq.gz --illumina -j 36 

    trim_galore --paired RNA_HS-60_1.fastq.gz RNA_HS-60_2.fastq.gz --illumina -j 36

    trim_galore --paired RNA_HS-05_1.fastq.gz RNA_HS-05_2.fastq.gz --illumina -j 36

We specified --illumina to indicate that the reads were generated using the Illumina sequencing platform.
Nonetheless, trim_galore can automatically detect adapters, providing flexibility for users who may know the specific sequencing platform.
Details of trimming are available in the trimming report file (Box 2).

Box 2: Example output of `trim_galore`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    In addition to the trimmed read files, `Trim_galore`` also generates a trimming report file.
    The trimming report contains details on read trimming, such as the number of trimmed reads.

.. code-block:: shell

    -rw-r--r-- 1 jinfang yinlab 4.3K May 30 00:25 MAG_HS-05_1.fastq.gz_trimming_report.txt 

    -rw-r--r-- 1 jinfang yinlab 1.4G May 30 00:28 MAG_HS-05_1_val_1.fq.gz 

    -rw-r--r-- 1 jinfang yinlab 4.5K May 30 00:28 MAG_HS-05_2.fastq.gz_trimming_report.txt 

    -rw-r--r-- 1 jinfang yinlab 1.4G May 30 00:28 MAG_HS-05_2_val_2.fq.gz 

    -rw-r--r-- 1 jinfang yinlab 4.3K May 30 00:24 MAG_HS-60_1.fastq.gz_trimming_report.txt 

    -rw-r--r-- 1 jinfang yinlab 2.7G May 30 00:30 MAG_HS-60_1_val_1.fq.gz 

    -rw-r--r-- 1 jinfang yinlab 4.5K May 30 00:30 MAG_HS-60_2.fastq.gz_trimming_report.txt 

    -rw-r--r-- 1 jinfang yinlab 2.7G May 30 00:30 MAG_HS-60_2_val_2.fq.gz 

    -rw-r--r-- 1 jinfang yinlab 4.3K May 30 00:25 RNA_HS-05_1.fastq.gz_trimming_report.txt 

    -rw-r--r-- 1 jinfang yinlab 1.4G May 30 00:28 RNA_HS-05_1_val_1.fq.gz 

    -rw-r--r-- 1 jinfang yinlab 4.5K May 30 00:28 RNA_HS-05_2.fastq.gz_trimming_report.txt 

    -rw-r--r-- 1 jinfang yinlab 1.4G May 30 00:28 RNA_HS-05_2_val_2.fq.gz 

    -rw-r--r-- 1 jinfang yinlab 4.3K May 30 00:24 RNA_HS-60_1.fastq.gz_trimming_report.txt 

    -rw-r--r-- 1 jinfang yinlab 2.7G May 30 00:30 RNA_HS-60_1_val_1.fq.gz 

    -rw-r--r-- 1 jinfang yinlab 4.5K May 30 00:30 RNA_HS-60_2.fastq.gz_trimming_report.txt 

    -rw-r--r-- 1 jinfang yinlab 2.7G May 30 00:30 RNA_HS-60_2_val_2.fq.gz 

.. warning::

    During the trimming process, certain reads may be entirely removed due to low quality in its entirety.
    Using the ``--retain_unpaired`` parameter in ``trim_galore`` allows for the preservation of single-end reads.
    In this protocol, this option was not selected, so that both reads of a forward-revise pair were removed.

P3. Assemble reads into contigs (TIMING ~84min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Use Megahit for assembling reads into contigs:

.. code-block:: shell

    megahit -m 0.5 -t 32 -o megahit_MAG_HS-60 -1 MAG_HS-60_1_val_1.fq.gz -2 MAG_HS-60_2_val_2.fq.gz --out-prefix MAG_HS-60 --min-contig-len 1000 

    megahit -m 0.5 -t 32 -o megahit_MAG_HS-05 -1 MAG_HS-05_1_val_1.fq.gz -2 MAG_HS-05_2_val_2.fq.gz --out-prefix MAG_HS-05 --min-contig-len 1000 


MEGAHIT generates two output folders: `megahit_MAG_HS-05`` and `megahit_MAG_HS-60`. Each contains five files and one sub-folder (Box 3). MAG_HS-60.contigs.fa is the final contig sequence file. We set --min-contig-len 1000, a common practice to retain all contigs longer than 1,000 base pairs.


Box 3: Example output of `MEGAHIT`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: shell

    -rw-r--r-- 1 jinfang yinlab  262 May 30 03:09 checkpoints.txt 

    -rw-r--r-- 1 jinfang yinlab    0 May 30 03:09 done 

    drwxr-xr-x 2 jinfang yinlab  41K May 30 03:09 intermediate_contigs 

    -rw-r--r-- 1 jinfang yinlab 135M May 30 03:09 MAG_HS-60.contigs.fa 

    -rw-r--r-- 1 jinfang yinlab 205K May 30 03:09 MAG_HS-60.log 

    -rw-r--r-- 1 jinfang yinlab 1.1K May 30 00:30 options.json 


P4. Predict genes by `pyrodigal` (TIMING ~1h)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    pyrodigal -p meta -i megahit_MAG_HS-60/MAG_HS-60.contigs.fa -d MAG_HS-60.cds -a MAG_HS-60.faa -f gff -o MAG_HS-60.gff -j 36

    pyrodigal -p meta -i megahit_MAG_HS-05/MAG_HS-05.contigs.fa -d MAG_HS-05.cds -a MAG_HS-05.faa -f gff -o MAG_HS-05.gff -j 36

The `pyrodigal` tool predicts genes from contigs. The `-p meta` parameter indicates that the input is metagenomic contigs. The output files include the coding sequence (CDS) file, protein sequence file, and gene location file (GFF) (Box 4).

Box 4: Example output of `pyrodigal`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: shell

    -rw-r--r-- 1 jinfang yinlab 138M May 30 04:47 MAG_HS-60.cds 

    -rw-r--r-- 1 jinfang yinlab  56M May 30 04:47 MAG_HS-60.gff 

    -rw-r--r-- 1 jinfang yinlab  61M May 30 04:47 MAG_HS-60.faa

if applied prodigal, please use the following command to fix id mapping:

.. code-block:: shell

    dbcan_utils gff_fix -i MAG_HS-60.faa -g MAG_HS-60.gff 

    dbcan_utils gff_fix -i MAG_HS-05.faa -g MAG_HS-05.gff 


Module 2. run_dbcan annotation to obtain CAZymes, CGCs, and substrates
```````````````````````````````````````````````````````````````````````````````

**CRITICAL STEP**

Users can skip P5 and P6, and directly run P7 (much slower though), if they want to predict not only CAZymes and CGCs, but also substrates.

P5. CAZyme annotation at the CAZyme family level (TIMING ~10min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    run_dbcan easy_CAZyme --input_raw_data MAG_HS-60.faa --mode protein --output_dir MAG_HS-60.CAZyme --db_dir db --input_format NCBI
    run_dbcan easy_CAZyme --input_raw_data MAG_HS-05.faa --mode protein --output_dir MAG_HS-05.CAZyme --db_dir db --input_format NCBI

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

    run_dbcan easy_CGC --input_raw_data MAG_HS-60.faa --mode protein --output_dir MAG_HS-60.PUL --input_format NCBI --gff_type prodigal --input_gff MAG_HS-60.gff --db_dir db
    run_dbcan easy_CGC --input_raw_data MAG_HS-05.faa --mode protein --output_dir MAG_HS-05.PUL --input_format NCBI --gff_type prodigal --input_gff MAG_HS-05.gff --db_dir db

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

    run_dbcan easy_substrate --input_raw_data MAG_HS-60.faa --mode protein --output_dir MAG_HS-60.PUL --input_format NCBI --gff_type prodigal --input_gff MAG_HS-60.gff --db_dir db
    run_dbcan easy_substrate --input_raw_data MAG_HS-05.faa --mode protein --output_dir MAG_HS-05.PUL --input_format NCBI --gff_type prodigal --input_gff MAG_HS-05.gff --db_dir db


Box 6. Example output folder content of run_dbcan substrate prediction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the `MAG_HS-05.dbCAN_ <https://bcb.unl.edu/dbCAN_tutorial/dataset5-Emilson2024/MAG_HS-05.dbCAN/>_` directory, a total of 17 files and 1 folder are generated:

    .. code-block:: shell

        -rw-rw-r--  1 jinfang jinfang  39M Nov  1 22:18 PUL_blast.out
        -rw-rw-r--  1 jinfang jinfang 3.1M Nov  1 22:15 CGC.faa
        -rw-rw-r--  1 jinfang jinfang 6.9M Nov  1 22:15 cgc.gff
        -rw-rw-r--  1 jinfang jinfang 702K Nov  1 22:15 cgc.out
        -rw-rw-r--  1 jinfang jinfang 321K Nov  1 22:15 cgc_standard.out
        -rw-rw-r--  1 jinfang jinfang 1.5M Nov  1 22:15 cgc_standard.out.json
        -rw-rw-r--  1 jinfang jinfang 556K Nov  1 22:14 dbcan-sub.hmm.out
        -rw-rw-r--  1 jinfang jinfang 345K Nov  1 22:14 diamond.out
        -rw-rw-r--  1 jinfang jinfang 455K Nov  1 22:14 dtemp.out
        -rw-rw-r--  1 jinfang jinfang 298K Nov  1 22:14 hmmer.out
        -rw-rw-r--  1 jinfang jinfang 270K Nov  1 22:15 overview.txt
        -rw-rw-r--  1 jinfang jinfang 1.1M Nov  1 22:15 stp.out
        -rw-rw-r--  1 jinfang jinfang  54K Nov  1 22:18 substrate.out
        drwxrwxr-x  2 jinfang jinfang  32K Nov  2 09:48 synteny.pdf
        -rw-rw-r--  1 jinfang jinfang 288K Nov  1 22:14 tf-1.out
        -rw-rw-r--  1 jinfang jinfang 237K Nov  1 22:14 tf-2.out
        -rw-rw-r--  1 jinfang jinfang 804K Nov  1 22:15 tp.out
        -rw-rw-r--  1 jinfang jinfang  31M Nov  1 21:07 uniInput


Module 3. Read mapping (Fig. 3) to calculate abundance for CAZyme families, subfamilies, CGCs, and substrates
``````````````````````````````````````````````````````````````````````````````````````````````````````````````

P8. Read mapping to all contigs of each sample (TIMING ~10 min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    mkdir samfiles 

    bwa index MAG_HS-05/MAG_HS-05.contigs.fa 

    bwa index megahit_MAG_HS-60/MAG_HS-60.contigs.fa 

    bwa mem -t 32 -o samfiles/MAG_HS-05.sam megahit_MAG_HS-05/MAG_HS-05.contigs.fa MAG_HS-05_1_val_1.fq.gz MAG_HS-05_2_val_2.fq.gz 

    bwa mem -t 32 -o samfiles/MAG_HS-60.sam megahit_MAG_HS-60/MAG_HS-60.contigs.fa MAG_HS-60_1_val_1.fq.gz MAG_HS-60_2_val_2.fq.gz 

    bwa mem -t 32 -o samfiles/RNA_HS-05.sam megahit_MAG_HS-05/MAG_HS-05.contigs.fa RNA_HS-05_1_val_1.fq.gz RNA_HS-05_2_val_2.fq.gz 

    bwa mem -t 32 -o samfiles/RNA_HS-60.sam megahit_MAG_HS-60/MAG_HS-60.contigs.fa RNA_HS-60_1_val_1.fq.gz RNA_HS-60_2_val_2.fq.gz 


P9. Read mapping to all contigs of each sample (TIMING ~10min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    cd samfiles 

    samtools sort -@ 32 -o MAG_HS-05.bam MAG_HS-05.sam 

    samtools sort -@ 32 -o MAG_HS-60.bam MAG_HS-60.sam 

    samtools sort -@ 32 -o RNA_HS-05.bam RNA_HS-05.sam 

    samtools sort -@ 32 -o RNA_HS-60.bam RNA_HS-60.sam 

    samtools index MAG_HS-05.bam 

    samtools index MAG_HS-60.bam 

    samtools index RNA_HS-05.bam 

    samtools index RNA_HS-60.bam 

    rm -rf *sam 

    cd .. 

P10. Read count calculation for all proteins of each sample using dbcan_utils (TIMING ~2min) 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell


    mkdir MAG_HS-05_abund && cd MAG_HS-05_abund 

    dbcan_utils cal_coverage -g ../MAG_HS-05.fix.gff -i ../samfiles/MAG_HS-05.bam -o MAG_HS-05.depth.txt -t 6 --overlap_base_ratio 0.2 --mapping_quality 30 ----identity 0.98  

    cd .. && mkdir MAG_HS-60_abund && cd MAG_HS-60_abund 

    dbcan_utils cal_coverage -g ../MAG_HS-60.fix.gff -i ../samfiles/MAG_HS-60.bam -o MAG_HS-60.depth.txt -t 6 --overlap_base_ratio 0.2 --mapping_quality 30 --–-identity 0.98  

    cd .. && mkdir RNA_HS-05_abund && cd RNA_HS-05_abund 

    dbcan_utils cal_coverage -g ../MAG_HS-05.fix.gff -i ../samfiles/RNA_HS-05.bam -o RNA_HS-05.depth.txt -t 6 --overlap_base_ratio 0.2 --mapping_quality 30 --–-identity 0.98 

    cd .. && mkdir RNA_HS-60_abund && cd MAG_HS-60_abund 

    dbcan_utils cal_coverage -g ../MAG_HS-60.fix.gff -i ../samfiles/RNA_HS-60.bam -o RNA_HS-60.depth.txt -t 6 --overlap_base_ratio 0.2 --mapping_quality 30 -–--identity 0.98 

    cd .. 

**Note**

How to filter the mapped reads are the key step of the abundance calculation. Several parameters can be set to filter the reads. The command dbcan_utils provides --overlap_base_ratio, --mapping_quality and --identity to filter the mapped reads which will may produce more reliable abundance. The default values for these three parameters are 0.2, 30 and 0.98, respectively. Read counts are saved in depth.txt files of each sample.  


P11. Read count calculation for a given region of contigs using Samtools (TIMING ~2min) 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    cd MAG_HS-05_abund 

    samtools depth -r k141_282901:2644-8574 ../samfiles/MAG_HS-05.bam > MAG_HS-05.cgc.depth.txt 


The parameter `-r k141_282901:2644-8574` specifies a region in a contig. For any CGC, its positional range can be found in the file cgc_standard.out produced by run_dbcan. The cgc.depth.txt files contain the raw read counts for the specified region. 


.. warning::

    The contig IDs are automatically generated by MEGAHIT. There is a small chance that a same contig ID appears in both samples. However, the two contigs in the two samples do not match each other even the ID is the same. For example, the contig ID k141_282901 is most likely only found in the MAG_HS-05 sample. Even if there is a k141_282901 in MAG_HS-60, the actual contigs in two samples are different.  

P12. dbcan_utils to calculate the abundance of CAZyme families, subfamilies, CGCs, and substrates (TIMING ~1min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    dbcan_utils fam_abund -bt MAG_HS-05.depth.txt -i ../MAG_HS-05.dbCAN -a TPM 

    dbcan_utils fam_substrate_abund -bt MAG_HS-05.depth.txt -i ../MAG_HS-05.dbCAN -a TPM 

    dbcan_utils CGC_abund -bt MAG_HS-05.depth.txt -i ../MAG_HS-05.dbCAN -a TPM 

    dbcan_utils CGC_substrate_abund -bt MAG_HS-05.depth.txt -i ../MAG_HS-05.dbCAN -a TPM 


    cd .. && cd RNA_HS-60_abund 

    dbcan_utils fam_abund -bt RNA_HS-05.depth.txt -i ../MAG_HS-05.dbCAN -a TPM 

    dbcan_utils fam_substrate_abund -bt RNA_HS-05.depth.txt -i ../MAG_HS-05.dbCAN -a TPM 

    dbcan_utils CGC_abund -bt RNA_HS-05.depth.txt -i ../MAG_HS-05.dbCAN -a TPM 

    dbcan_utils CGC_substrate_abund -bt RNA_HS-05.depth.txt -i ../MAG_HS-05.dbCAN -a TPM 


    cd .. && cd MAG_HS-60_abund 

    dbcan_utils fam_abund -bt MAG_HS-60.depth.txt -i ../MAG_HS-60.dbCAN -a TPM 

    dbcan_utils fam_substrate_abund -bt MAG_HS-60.depth.txt -i ../MAG_HS-60.dbCAN -a TPM 

    dbcan_utils CGC_abund -bt MAG_HS-60.depth.txt -i ../MAG_HS-60.dbCAN -a TPM 

    dbcan_utils CGC_substrate_abund -bt MAG_HS-60.depth.txt -i ../MAG_HS-60.dbCAN -a TPM 

    cd .. && cd RNA_HS-60_abund 

    dbcan_utils fam_abund -bt RNA_HS-60.depth.txt -i ../MAG_HS-60.dbCAN -a TPM 

    dbcan_utils fam_substrate_abund -bt RNA_HS-60.depth.txt -i ../MAG_HS-60.dbCAN -a TPM 

    dbcan_utils CGC_abund -bt RNA_HS-60.depth.txt -i ../MAG_HS-60.dbCAN -a TPM 

    dbcan_utils CGC_substrate_abund -bt RNA_HS-60.depth.txt -i ../MAG_HS-60.dbCAN -a TPM 


We developed a set of Python scripts as ``dbcan_utils`` (included in the ``run_dbcan`` package) to take the raw read counts for all genes as input and output the normalized abundances (refer to Box 7) of CAZyme families, subfamilies, CGCs, and substrates (see Fig. 4). The parameter ``-a TPM`` can also be set to two other metrics: RPM, or RPKM61.

- **RPKM** is calculated as the number of mapped reads to a gene G divided by [(total number of mapped reads to all genes / 10^6) x (gene G length / 1000)].
- **RPM** is the number of mapped reads to a gene G divided by (total number of mapped reads to all genes / 10^6).
- **TPM** is calculated as [number of mapped reads to a gene G / (gene G length / 1000)] divided by the sum of [number of mapped reads to each gene / (the gene length / 1000)].


Box 7. Example output of dbcan_utils
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As an example, `RNA_HS-05_abund <https://bcb.unl.edu/dbCAN_tutorial/dataset6-Emilson2024/RNA_HS-05_abund/>_` folder has 7 TSV files:

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


P13. Heatmap for CAZyme substrate abundance across samples (Fig. S4B) (TIMING 1min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    dbcan_plot heatmap_plot --samples MAG_HS-05,MAG_HS-60,RNA_HS-05,RNA_HS-60 -i MAG_HS-05_abund/fam_substrate_abund.out,MAG_HS-60_abund/fam_substrate_abund.out,RNA_HS-05_abund/fam_substrate_abund.out,RNA_HS-60_abund/fam_substrate_abund.out --show_abund --top 20 

Here we plot the top 20 substrates in the two samples. The input files are the two CAZyme substrate abundance files calculated based on dbCAN-sub result. The default heatmap is ranked by substrate abundances. To rank the heatmap according to abundance profile using the function clustermap of seaborn package, users can invoke the ``--cluster_map`` parameter.

P14. Barplot for CAZyme family/subfamily/EC abundance across samples (TIMING 1min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    dbcan_plot bar_plot --samples MAG_HS-05,MAG_HS-60,RNA_HS-05,RNA_HS-60 --vertical_bar --top 20 -i MAG_HS-05_abund/fam_abund.out,MAG_HS-60_abund/fam_abund.out,RNA_HS-05_abund/fam_abund.out,RNA_HS-60_abund/fam_abund.out --pdf fam.pdf 

    dbcan_plot bar_plot --samples MAG_HS-05,MAG_HS-60,RNA_HS-05,RNA_HS-60 --vertical_bar --top 20 -i MAG_HS-05_abund/subfam_abund.out,MAG_HS-60_abund/subfam_abund.out,RNA_HS-05_abund/subfam_abund.out,RNA_HS-60_abund/subfam_abund.out --pdf subfam.pdf 

    dbcan_plot bar_plot --samples MAG_HS-05,MAG_HS-60,RNA_HS-05,RNA_HS-60 --vertical_bar --top 20 -i MAG_HS-05_abund/EC_abund.out,MAG_HS-60_abund/EC_abund.out,RNA_HS-05_abund/EC_abund.out,RNA_HS-60_abund/EC_abund.out --pdf ec.pdf 


Users can choose to generate a barplot instead of heatmap using the ``bar_plot`` method.

P15. Synteny plot between a CGC and its best PUL hit with read mapping coverage to CGC (TIMING 1min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    dbcan_plot CGC_synteny_coverage_plot -i MAG_HS-05.dbCAN --readscount MAG_HS-05_abund/MAG_HS-05.cgc.depth.txt --cgcid 'k141_282901|CGC1'  

The `MAG_HS-05.dbCAN`` folder contains the `PUL.out`` file. Using this file, the cgc_standard.out file, and the best PUL's gff file in `dbCAN-PUL.tar.gz`, the CGC_synteny_plot method will create the CGC-PUL synteny plot. The --cgcid parameter is required to specify which CGC to plot (k141_282901|CGC1 in this example). The MAG_HS-05.cgc.depth.txt file is used to plot the read mapping coverage. 


If users only want to plot the CGC structure:

.. code-block:: shell

    dbcan_plot CGC_plot -i MAG_HS-05.dbCAN --cgcid 'k141_282901|CGC1' 

If users only want to plot the CGC structure plus the read mapping coverage:

.. code-block:: shell

    dbcan_plot CGC_coverage_plot -i MAG_HS-05.dbCAN --cgcid 'k141_282901|CGC1'  --readscount MAG_HS-05_abund/MAG_HS-05.cgc.depth.txt

If users only want to plot the synteny between the CGC and PUL:

.. code-block:: shell

    dbcan_plot CGC_synteny_plot -i MAG_HS-05.dbCAN --cgcid 'k141_282901|CGC1'  


.. warning::

    The CGC IDs in different samples do not match each other. For example, specifying ``-i MAG_HS-05.dbCAN`` is to plot
    the ``'k141_282901|CGC1'`` in the fefifo_8022_1 sample. The ``'k141_282901|CGC1'`` in the fefifo_8022_7 sample most likely does not exist,
    and even it does, the CGC has a different sequence even if the ID is the same.
