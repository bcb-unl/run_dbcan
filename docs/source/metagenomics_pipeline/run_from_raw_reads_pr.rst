Run from Raw Reads(Priest 2023): Automated CAZyme and Glycan Substrate Annotation in Microbiomes: A Step-by-Step Protocol
=========================================================================================================================


.. _priest_2023: https://www.nature.com/articles/s43705-023-00324-7

Example 3: Priest2023 Dataset `priest_2023`_
-------------------------------------------------

The Priest2023 dataset `priest_2023`_ was published in 2023 from a study of glycan utilization in the particulate and dissolved organic matter pools (POM and DOM) in the ocean waters of the Arctic during late summer1. This study applied PacBio HiFi long read sequencing for metagenomes (MGs) of 8 samples (4 sites of two depth) and Illumina short read sequencing for metatranscriptomes (MTs) of 4 samples (4 sites). It investigated CAZyme abundance and transcription profiles as well as inferred glycan degradation activities by comparing MGs and MTs from samples of different sites and depth. The sites are categorized into above-slope, above-shelf and open-ocean groups, and the water depth includes surface water (SRF) and bottom of the surface mixed layer (BML). CAZyme abundance and expression showed different profiles for the degradation of complex algal glycans such as sulfated fucan, laminarin and xylan in different site groups and water depths.

For this protocol, we will select two MG samples (HiFi DNA long reads, FRAM_WSC20_S25_SRF_MG and FRAM_WSC20_S25_BML_MG) and two MT samples (FRAM_WSC20_S25_SRF_MT and FRAM_WSC20_S25_BML_MT, paired-end 2x151bp reads) for CAZyme annotation. We will assemble HiFi long MG reads into contigs, which will be used to predict CAZymes and glycan substrates. We will map the MT reads to the contigs to calculate expression levels.
Procedure
---------


Software and versions
`````````````````````

- **Anaconda** (`Anaconda <https://www.anaconda.com>`_, version 23.7.3)
- **Flye** (`Flye <https://github.com/mikolmogorov/Flye>`_, version 2.9.5)
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

    conda env create -f CAZyme_annotation_pr.yml
    conda activate CAZyme_annotation_pr

To install the databases, execute the following commands:

.. include:: prepare_the_database.rst





Module 1: Reads processing (Fig. 3) to obtain contigs
`````````````````````````````````````````````````````

P1. Contamination Check (TIMING ~10min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Download the Priest2023 dataset:


.. code-block:: shell

    wget https://bcb.unl.edu/dbCAN_tutorial/dataset3-Priest2023/BML_MG.fastq.gz
    wget https://bcb.unl.edu/dbCAN_tutorial/dataset3-Priest2023/SRF_MG.fastq.gz
    wget https://bcb.unl.edu/dbCAN_tutorial/dataset3-Priest2023/BML_MT_1.fastq.gz
    wget https://bcb.unl.edu/dbCAN_tutorial/dataset3-Priest2023/BML_MT_2.fastq.gz
    wget https://bcb.unl.edu/dbCAN_tutorial/dataset3-Priest2023/SRF_MT_1.fastq.gz
    wget https://bcb.unl.edu/dbCAN_tutorial/dataset3-Priest2023/SRF_MT_2.fastq.gz

Use `kraken2` to check for contaminated reads:


.. code-block:: shell

    kraken2 --threads 32 --quick --paired --db K2 --report SRF_MT.kreport --output SRF_MT.kraken.output SRF_MT_1.fastq.gz SRF_MT_2.fastq.gz
    kraken2 --threads 32 --quick --paired --db K2 --report BML_MT.kreport --output BML_MT.kraken.output BML_MT_1.fastq.gz BML_MT_2.fastq.gz

    kraken2 --threads 32 --quick --paired --db K2 --report SRF_MG.kreport --output SRF_MG.kraken.output SRF_MG.fastq.gz
    kraken2 --threads 32 --quick --paired --db K2 --report BML_MG.kreport --output BML_MG.kraken.output BML_MG.fastq.gz

Kraken2 found much contamination (Box 1) from human in the Priest2023 data. Consequently, human reads need to be removed before assembly.

Reads can be aligned to the reference genomes of potential contamination source organisms to remove the aligned reads. The most common source in microbiome studies is from human.

Box 1: Example of Kraken2 output files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    The `kreport` files can be examined to identify potential contamination source organisms.

    .. code-block:: shell

        -rw-rw-r-- 1 jinfang jinfang   54M Dec 27 07:42 BML_MG.kraken.output
        -rw-rw-r-- 1 jinfang jinfang  1.2M Dec 27 07:42 BML_MG.kreport
        -rw-rw-r-- 1 jinfang jinfang  3.4G Dec 27 08:01 BML_MT.kraken.output
        -rw-rw-r-- 1 jinfang jinfang 1023K Dec 27 08:02 BML_MT.kreport
        -rw-rw-r-- 1 jinfang jinfang   61M Dec 27 07:39 SRF_MG.kraken.output
        -rw-rw-r-- 1 jinfang jinfang  1.2M Dec 27 07:39 SRF_MG.kreport
        -rw-rw-r-- 1 jinfang jinfang  2.6G Dec 27 07:50 SRF_MT.kraken.output
        -rw-rw-r-- 1 jinfang jinfang  1.1M Dec 27 07:51 SRF_MT.kreport


P2. Remove contamination reads from human (TIMING ~40min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

From the Kraken2 output files, we identified humans as the contamination source, we can use the following commands to remove the contamination reads by aligning reads to the human reference genome.

.. code-block:: shell

    mkdir hg38 && cd hg38 && wget https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
    cd .. && mkdir contamination && cd contamination
    minimap2 -a -x map-hifi -MD -t 32 -o SRF_MG.hg38.sam ../hg38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz ../SRF_MG.fastq.gz
    minimap2 -a -x map-hifi -MD -t 32 -o BML_MG.hg38.sam ../hg38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz ../BML_MG.fastq.gz
    samtools fastq -f 4 -@ 32 -0 ../SRF_MG.clean.fq.gz SRF_MG.hg38.sam
    samtools fastq -f 4 -@ 32 -0 ../BML_MG.clean.fq.gz BML_MG.hg38.sam
    bwa mem ../hg38/hg38 ../SRF_MT_1.fastq.gz ../SRF_MT_2.fastq.gz -t 32 -o SRF_MT.hg38.sam
    bwa mem ../hg38/hg38 ../BML_MT_1.fastq.gz ../BML_MT_2.fastq.gz -t 32 -o BML_MT.hg38.sam
    samtools fastq -f 12 -@ 32 -1 ../SRF_MT_1.clean.fq.gz -2 ../SRF_MT_2.clean.fq.gz SRF_MT.hg38.sam
    samtools fastq -f 12 -@ 32 -1 ../BML_MT_1.clean.fq.gz -2 ../BML_MT_2.clean.fq.gz BML_MT.hg38.sam
    cd ..

Or we could use KrakenTools to extract reads.

.. code-block:: shell

    extract_kraken_reads.py -k SRF_MT.kraken.output -taxid 9606 -exclude -s1 SRF_MT_1.fastq.gz -s2 SRF_MT_2.fastq.gz -o SRF_MT_1.clean.fastq -o2 SRF_MT_2.clean.fastq
    gzip SRF_MT_1.clean.fastq
    gzip SRF_MT_2.clean.fastq

    extract_kraken_reads.py -k BML_MT.kraken.output -taxid 9606 -exclude -s1 BML_MT_1.fastq.gz -s2 BML_MT_2.fastq.gz -o BML_MT_1.clean.fastq -o2 BML_MT_2.clean.fastq
    gzip BML_MT_1.clean.fastq
    gzip BML_MT_2.clean.fastq

    extract_kraken_reads.py -k BML_MG.kraken.output -taxid 9606 -exclude -s BML_MG_.fastq.gz -o BML_MG_.clean.fastq
    gzip BML_MG.clean.fastq

    extract_kraken_reads.py -k SRF_MG.kraken.output -taxid 9606 -exclude -s SRF_MG.fastq.gz -o SRF_MG.clean.fastq
    gzip SRF_MG.clean.fastq


P3| Trim adapter and low-quality reads (TIMING ~20min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


.. code-block:: shell

    trim_galore --illumina -j 8 --paired BML_MT_1.clean.fastq.gz BML_MT_2.clean.fastq.gz
    trim_galore --illumina -j 8 --paired SRF_MT_1.clean.fastq.gz SRF_MT_2.clean.fastq.gz

The HiFi long reads do not need to be trimmed. Hence, this step only applies to MT illumina short read data. We specified --illumina to indicate that the reads were generated using the Illumina sequencing platform. Nonetheless, trim_galore possesses the ability to automatically detect the adapter, providing flexibility in adapter handling for users who may know the specific sequencing platform. Details of trimming are available in the trimming report file (Box 2).

Box 2: Example output of trim_galore
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    In addition to the trimmed read files, Trim_galore also generates a trimming report file. The trimming report contains details on read trimming, such as the number of trimmed reads.

    .. code-block:: shell

        -rw-rw-r-- 1 jinfang jinfang 4.2K Dec 28 21:56 BML_MT_1.clean.fq.gz_trimming_report.txt
        -rw-rw-r-- 1 jinfang jinfang 2.3G Dec 28 22:05 BML_MT_1.clean_val_1.fq.gz
        -rw-rw-r-- 1 jinfang jinfang 4.7K Dec 28 22:05 BML_MT_2.clean.fq.gz_trimming_report.txt
        -rw-rw-r-- 1 jinfang jinfang 3.0G Dec 28 22:05 BML_MT_2.clean_val_2.fq.gz
        -rw-rw-r-- 1 jinfang jinfang 4.9K Dec 28 10:07 SRF_MT_1.clean.fq.gz_trimming_report.txt
        -rw-rw-r-- 1 jinfang jinfang 2.7G Dec 28 10:19 SRF_MT_1.clean_val_1.fq.gz
        -rw-rw-r-- 1 jinfang jinfang 5.1K Dec 28 10:19 SRF_MT_2.clean.fq.gz_trimming_report.txt
        -rw-rw-r-- 1 jinfang jinfang 3.3G Dec 28 10:19 SRF_MT_2.clean_val_2.fq.gz

.. warning::

    During the trimming process, certain reads may be entirely removed due to low quality in its entirety. Using the `--retain_unpaired` parameter in trim_galore allows for the preservation of single-end reads. In this protocol, this option was not selected, so that both reads of a forward-revise pair were removed.



P4. Assemble HiFi reads into metagenome (TIMING ~4h20min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Flye was used to assemble the HiFi long reads into contigs.

.. code-block:: shell

    flye --threads 32 --meta --pacbio-hifi BML_MG.clean.fq.gz --hifi-error 0.01 --keep-haplotypes --out-dir flye_BML_MG
    flye --threads 32 --meta --pacbio-hifi SRF_MG.clean.fq.gz --hifi-error 0.01 --keep-haplotypes --out-dir flye_SRF_MG

Flye generates two folders `flye_BML_MG` and `flye_SRF_MG`. Each folder
contains 6 files and 5 sub-folders (Box 3), among them `assembly.fasta` is the final contig sequence file.
We set `--hifi-error` 0.01, a generally accepted error rate of HiFi sequencing.
Parameter `--meta` is set to assemble reads into metagenomes.

Box 3: Example output of Flye
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    .. code-block:: shell

        drwxrwxr-x  2 jinfang jinfang 4.0K Dec 27 20:15 00-assembly
        drwxrwxr-x  2 jinfang jinfang 4.0K Dec 27 20:43 10-consensus
        drwxrwxr-x  2 jinfang jinfang 4.0K Dec 27 21:14 20-repeat
        drwxrwxr-x  2 jinfang jinfang 4.0K Dec 27 21:16 30-contigger
        drwxrwxr-x  2 jinfang jinfang 4.0K Dec 27 22:06 40-polishing
        -rw-rw-r--  1 jinfang jinfang 314M Dec 27 22:06 assembly.fasta
        -rw-rw-r--  1 jinfang jinfang 311M Dec 27 22:06 assembly_graph.gfa
        -rw-rw-r--  1 jinfang jinfang 6.6M Dec 27 22:06 assembly_graph.gv
        -rw-rw-r--  1 jinfang jinfang 867K Dec 27 22:06 assembly_info.txt
        -rw-rw-r--  1 jinfang jinfang  61M Dec 27 22:06 flye.log
        -rw-rw-r--  1 jinfang jinfang   92 Dec 27 22:06 params.json


P5. Predict genes by Prodigal (~14min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    pyrodigal -p meta -i flye_BML_MG/assembly.fasta -d BML_MG.cds -a BML_MG.faa -f gff -o BML_MG.gff -j 36
    pyrodigal -p meta -i flye_SRF_MG/assembly.fasta -d SRF_MG.cds -a SRF_MG.faa -f gff -o SRF_MG.gff -j 36


Box 3: Example output of Prokka
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: shell

    -rw-rw-r--  1 jinfang jinfang 105M Dec 27 23:26 BML_MG.faa
    -rw-rw-r--  1 jinfang jinfang 314M Dec 27 22:06 BML_MG.cds
    -rw-rw-r--  1 jinfang jinfang 467M Dec 27 23:26 BML_MG.gff


Module 1: run_dbcan annotation to obtain CAZymes, CGCs, and substrates
```````````````````````````````````````````````````````````````````````````````

Hint: If you apply ``prodigal`` not ``pyrodigal``, use this to fix gene ID in gff files by dbcan_utils (TIMING ~1min)

.. code-block:: shell

    dbcan_utils gff_fix -i BML_MG.faa -g BML_MG.gff
    dbcan_utils gff_fix -i SRF_MG.faa -g SRF_MG.gff


Users can skip P6 and P7, and directly run P8 (much slower though), if they want to predict not only CAZymes and CGCs, but also substrates.

P6. CAZyme annotation at family level (TIMING ~10min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    run_dbcan easy_CAZyme --input_raw_data BML_MG.faa --mode protein --output_dir BML_MG.CAZyme --db_dir db --input_format NCBI
    run_dbcan easy_CAZyme --input_raw_data SRF_MG.faa --mode protein --output_dir SRF_MG.CAZyme --db_dir db --input_format NCBI


Two arguments are required for ``run_dbcan``: the input sequence file (faa files) and the sequence type (protein).
By default, ``run_dbcan`` will use three methods (``pyHMMER`` vs ``dbCAN HMMdb``, ``DIAMOND`` vs ``CAZy``, ``pyHMMER`` vs ``dbCAN-sub HMMdb``) for
CAZyme annotation.


The sequence type can be `protein`, `prok`, `meta`. If the input sequence file contains metagenomic contig sequences (`fna` file),
the sequence type has to be `meta`, and `prodigal` will be called to predict genes.


P7. CGC prediction (TIMING ~15 min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following commands will re-run run_dbcan to not only predict CAZymes but also CGCs with protein `faa` and gene location `gff` files.

.. code-block:: shell

    run_dbcan easy_CGC --input_raw_data BML_MG.faa --mode protein  --output_dir BML_MG.PUL --input_format NCBI --gff_type prodigal --input_gff BML_MG.gff --db_dir db
    run_dbcan easy_CGC --input_raw_data SRF_MG.faa --mode protein  --output_dir SRF_MG.PUL --input_format NCBI --gff_type prodigal --input_gff SRF_MG.gff --db_dir db

.. warning::

    **Creating own gff file**
    If the users would like to create their own ``gff`` file (instead of using Prokka or Prodigal),
    it is important to make sure the value of ID attribute in the ``gff`` file matches the protein ID in the protein ``faa`` file.

    **[Troubleshooting]CGC not found**
    If no result is found in CGC output file, it is most likely because the sequence IDs in ``gff`` file and ``faa`` file do not match.
    Another less likely reason is that the contigs are too short and fragmented and not suitable for CGC prediction.


P8. Substrate prediction for CAZymes and CGCs (TIMING ~5h)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following commands will re-run run_dbcan to predict CAZymes, CGCs.

.. code-block:: shell

    run_dbcan easy_substrate --input_raw_data BML_MG.faa --mode protein  --output_dir BML_MG.PUL --input_format NCBI --gff_type prodigal --input_gff BML_MG.gff --db_dir db
    run_dbcan easy_substrate --input_raw_data SRF_MG.faa --mode protein  --output_dir SRF_MG.PUL --input_format NCBI --gff_type prodigal --input_gff SRF_MG.gff --db_dir db


Box 6: Example output folder content of run_dbcan substrate prediction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    In the output directory (https://bcb.unl.edu/dbCAN_tutorial/dataset3-Priest2023/BML_MG.dbCAN/), a total of 17 files and 1 folder are generated:

    .. code-block:: shell

        -rw-rw-r--  1 jinfang jinfang  9.6M Dec 28 10:18 PUL_blast.out
        -rw-rw-r--  1 jinfang jinfang  1.8M Dec 28 10:18 CGC.faa
        -rw-rw-r--  1 jinfang jinfang   26M Dec 28 10:18 cgc.gff
        -rw-rw-r--  1 jinfang jinfang  450K Dec 28 10:18 cgc.out
        -rw-rw-r--  1 jinfang jinfang  212K Dec 28 10:18 cgc_standard.out
        -rw-rw-r--  1 jinfang jinfang 1005K Dec 28 10:18 cgc_standard.out.json
        -rw-rw-r--  1 jinfang jinfang  406K Dec 28 10:11 dbcan-sub.hmm.out
        -rw-rw-r--  1 jinfang jinfang  325K Dec 28 10:11 diamond.out
        -rw-rw-r--  1 jinfang jinfang  332K Dec 28 10:11 dtemp.out
        -rw-rw-r--  1 jinfang jinfang  220K Dec 28 10:11 hmmer.out
        -rw-rw-r--  1 jinfang jinfang  240K Dec 28 10:18 overview.txt
        -rw-rw-r--  1 jinfang jinfang  1.7M Dec 28 10:17 stp.out
        -rw-rw-r--  1 jinfang jinfang   17K Dec 28 10:18 substrate.out
        drwxrwxr-x  2 jinfang jinfang   12K Dec 28 10:19 synteny.pdf
        -rw-rw-r--  1 jinfang jinfang  293K Dec 28 10:13 tf-1.out
        -rw-rw-r--  1 jinfang jinfang  222K Dec 28 10:15 tf-2.out
        -rw-rw-r--  1 jinfang jinfang  1.7M Dec 28 10:17 tp.out
        -rw-rw-r--  1 jinfang jinfang  105M Dec 28 05:57 uniInput


Descriptions of Output Files:

    - ``PUL_blast.out``: BLAST results between CGCs and PULs.
    - ``CGC.faa``: CGC Fasta sequences.
    - ``cgc.gff``: reformatted from the user input gff file by marking CAZymes, TFs, TCs, and STPs.
    - ``cgc.out``: raw output of CGC predictions.

Each entry in cgc.out includes:

1.	CGC_id: CGC1
2.	type: CAZyme
3.	contig_id: contig_10157
4.	gene_id: BML_MG_01992
5.	start: 33003
6.	end: 36077
7.	strand: +
8.	annotation: GH2

Explanation: the gene BML_MG_01992 encodes a GH2 CAZyme in the CGC1 of the contig contig_10157. CGC1 also has other genes, which are provided in other rows. BML_MG_01992 is on the positive strand of contig_10157 from 33003 to 36077. The type can be one of the four signature gene types (CAZymes, TCs, TFs, STPs) or the null type (not annotated as one of the four signature genes).

`cgc_standard.out.json`: JSON format of cgc_standard.out.
`dbcan-sub.hmm.out`: HMMER search result against dbCAN-sub HMMdb, including a column with CAZyme substrates extracted from fam-substrate-mapping-08012023.tsv.
`diamond.out`: DIAMOND search result against the CAZy annotated protein sequences (CAZyDB.07262023.fa).
`dtemp.out`: temporary file.
`hmmer.out`: HMMER search result against dbCAN HMMdb.
`overview.txt`: summary of CAZyme annotation from three methods in TSV format. An example row has the following columns:
1.	Gene_ID: BML_MG_01761
2.	EC#: 2.4.99.-:5
3.	dbCAN: GT112(19-370)
4.	dbCAN_sub: GT112_e0
5.	DIAMOND: GT112
6.	#ofTools: 3
Explanation: the protein BML_MG_01761 is annotated by 3 tools to be a CAZyme: (1) GT112 (CAZy defined family GT112) by HMMER vs dbCAN HMMdb with a domain range from aa position 19 to 370, (2) GT112_e0 (eCAMI defined subfamily e0; e indicates it is from eCAMI not CAZy) by HMMER vs dbCAN-sub HMMdb (derived from eCAMI subfamilies), and (3) GT112 by DIAMOND vs CAZy annotated protein sequences. The second column 2.4.99.-:5 is extracted from eCAMI, meaning that the eCAMI subfamily GT112_e0 contains 5 member proteins which have an EC 2.4.99.- according to CAZy. In most cases, the 3 tools will have the same CAZyme family assignment. When they give different assignment. We recommend a preference order: dbCAN > eCAMI/dbCAN-sub > DIAMOND. See our dbCAN2 paper2, dbCAN3 paper3, and eCAMI4 for more details.
Note: If users invoked the --use_signalP parameter when running run_dbcan, there will be an additional column called signal in the overview.txt.
stp.out: HMMER search result against the MiST5 compiled signal transduction protein HMMs from Pfam.
tf-1.out: HMMER search result against the DBD6 compiled transcription factor HMMs from Pfam 7.
tf-2.out: HMMER search result against the DBD compiled transcription factor HMMs from Superfamily 8.
tp.out: DIAMOND search result against the TCDB 9 annotated protein sequences.
substrate.out: summary of substrate prediction results for CGCs in TSV format from two approaches3 (dbCAN-PUL blast search and dbCAN-sub majority voting). An example row has the following columns:
1.	CGC_ID: contig_10778|CGC2
2.	Best hit PUL_ID in dbCAN-PUL: PUL0400
3.	Substrate of the hit PUL: alginate
4.	Sum of bitscores for homologous gene pairs between CGC and PUL: 851.0
5.	Types of homologous gene pairs: CAZyme-CAZyme;CAZyme-CAZyme;CAZyme-CAZyme;CAZyme-CAZyme
6.	Substrate predicted by majority voting of CAZymes in CGC: alginate
7.	Voting score: 2.0
Explanation: The CGC2 of contig_10778 has its best hit PUL0400 (from PUL_blast.out) with alginate as substrate (from dbCAN-PUL_12-12-2023.xlsx). Four signature genes are matched between contig_10778|CGC2 and PUL0400 (from PUL_blast.out): all the four are CAZymes. The sum of blast bitscores of the 4 homologous pairs (CAZyme-CAZyme;CAZyme-CAZyme;CAZyme-CAZyme;CAZyme-CAZyme) is 851.0. Hence, the substrate of contig_10778|CGC2 is predicted to be alginate according to dbCAN-PUL blast search. The last two columns are based on the dbCAN-sub result (dbcan-sub.hmm.out), according to which two CAZymes in contig_10778|CGC2 are predicted to have alginate substrate. The voting score is thus 2.0, so that according to the majority voting rule, contig_10778|CGC2 is predicted to have an alginate substrate.
Note: for many CGCs, only one of the two approaches produces substrate prediction. In some cases, the two approaches produce different substrate assignments. We recommend a preference order: dbCAN-PUL blast search > dbCAN-sub majority voting. See our dbCAN3 paper3 for more details.
synteny.pdf: a folder with syntenic block alignment plots between all CGCs and PULs.
uniInput: renamed Fasta file from input protein sequence file.





Module 3. Read mapping (Fig. 3) to calculate abundance for CAZyme families, subfamilies, CGCs, and substrates
``````````````````````````````````````````````````````````````````````````````````````````````````````````````

P9. Read mapping to all contigs of each sample (TIMING ~10 min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    mkdir samfiles
    bwa index flye_BML_MG/assembly.fasta 
    bwa index flye_SRF_MG/assembly.fasta
    
    bwa mem -t 32 -o samfiles/HCZGU.sam flye_BML_MG/assembly.fasta HCZGU_1_val_1.fq.gz HCZGU_2_val_2.fq.gz 
    bwa mem -t 32 -o samfiles/HTOOS.sam flye_SRF_MG/assembly.fasta HTOOS_1_val_1.fq.gz HTOOS_2_val_2.fq.gz 

    bwa mem -t 32 -o samfiles/BML_MT.sam flye_BML_MG/assembly.fasta BML_MT_1.clean_val_1.fq.gz BML_MT_2.clean_val_2.fq.gz 
    bwa mem -t 32 -o samfiles/SRF_MT.sam flye_SRF_MG/assembly.fasta SRF_MT_1.clean_val_1.fq.gz SRF_MT_2.clean_val_2.fq.gz 


P10. HiFi MG long read mapping to all contigs of each sample (TIMING ~20min) 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. code-block:: shell

    minimap2 -a -x map-hifi --MD -t 32 -o samfiles/BML_MG.sam flye_BML_MG/assembly.fasta BML_MG.clean.fq.gz 

    minimap2 -a -x map-hifi --MD -t 32 -o samfiles/SRF_MG.sam flye_SRF_MG/assembly.fasta SRF_MG.clean.fq.gz 

The --MD parameter is required in this command, because step P12 will access the MD tag in SAM files generated in this step to calculate the coverage.  

P11. Convert SAM to BAM and sort (TIMING ~5min)

.. code-block:: shell

    cd samfiles 

    samtools sort -@ 32 -o SRF_MG.bam SRF_MG.sam 

    samtools sort -@ 32 -o BML_MG.bam BML_MG.sam 

    samtools sort -@ 32 -o BML_MT.bam BML_MT.sam 

    samtools sort -@ 32 -o SRF_MT.bam SRF_MT.sam 

    samtools index BML_MT.bam 

    samtools index SRF_MT.bam 

    samtools index BML_MG.bam 

    samtools index SRF_MG.bam 

    rm -rf *sam 

    cd .. 


P12. Read count calculation for all proteins of each sample using dbcan_utilsBedtools (TIMING ~6min12min) 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Calculate HiFi MG read abundance for proteins in the two samples: 

.. code-block:: shell


    mkdir SRF_MG_abund && cd SRF_MG_abund 

    dbcan_utils cal_coverage -g ../SRF_MG.fix.gff -i ../samfiles/SRF_MG.bam -o SRF_MG.depth.txt -t 6 --overlap_base_ratio 0.2 --mapping_quality 30 –-identity 0.98 --hifi 

    cd .. 



    mkdir BML_MG_abund && cd BML_MG_abund 

    dbcan_utils cal_coverage -g ../BML_MG.fix.gff -i ../samfiles/BML_MG.bam -o BML_MG.depth.txt -t 6 --overlap_base_ratio 0.2 --mapping_quality 30 –-identity 0.98 --hifi 

    cd .. 



Calculate Illumina MT read abundance for proteins in the two samples: 

.. code-block:: shell


    mkdir BML_MT_abund && cd BML_MT_abund 

    dbcan_utils cal_coverage -g ../BML_MG.fix.gff -i ../samfiles/BML_MT.bam -o BML_MT.depth.txt -t 6 --overlap_base_ratio 0.2 --mapping_quality 30 –-identity 0.98 

    cd .. 



    mkdir SRF_MT_abund && cd SRF_MT_abund 

    dbcan_utils cal_coverage -g ../SRF_MG.fix.gff -i ../samfiles/SRF_MT.bam -o SRF_MT.depth.txt -t 6 --overlap_base_ratio 0.2 --mapping_quality 30 –-identity 0.98 

    cd .. 



Read counts are saved in depth.txt files of each sample.  

P13. Read count calculation for a given region of contigs using Samtools (TIMING ~1min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    cd SRF_MG_abund 

    samtools depth -r contig_12725:10309-25115 ../samfiles/SRF_MT.bam > SRF_MT.cgc.depth.txt 

    cd .. && cd SRF_MT_abund 

    samtools depth -r contig_12725:10309-25115 ../samfiles/SRF_MG.bam > SRF_MG.cgc.depth.txt 

    cd .. 

The parameter `-r contig_12725:10309-25115contig_6487:2406-6727`` specifies a region in a contig. For any CGC, its positional range can be found in the file cgc_standard output produced by run_dbcan. The depth.txt files contain the raw read counts for the specified region. 



P14. dbcan_utils to calculate the abundance of CAZyme families, subfamilies, CGCs, and substrates (TIMING ~1min) 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
As the Priest2023 dataset has both MG and MT data, for each sample we will have two abundance folders, e.g., BML_MG_abund, BML_MT_abund for BML sample, and SRF_MG_abund, SRF_MT_abund for SRF sample. 

.. code-block:: shell

    cd BML_MT_abund 

    dbcan_utils fam_abund -bt BML_MT.depth.txt -i ../BML_MG.dbCAN -a TPM 

    dbcan_utils fam_substrate_abund -bt BML_MT.depth.txt -i ../BML_MG.dbCAN -a TPM 

    dbcan_utils CGC_abund -bt BML_MT.depth.txt -i ../BML_MG.dbCAN -a TPM 

    dbcan_utils CGC_substrate_abund -bt BML_MT.depth.txt -i ../BML_MG.dbCAN -a TPM 


    cd .. && cd SRF_MT_abund 

    dbcan_utils fam_abund -bt SRF_MT.depth.txt -i ../SRF_MG.dbCAN -a TPM 

    dbcan_utils fam_substrate_abund -bt SRF_MT.depth.txt -i ../SRF_MG.dbCAN -a TPM 

    dbcan_utils CGC_abund -bt SRF_MT.depth.txt -i ../SRF_MG.dbCAN -a TPM 

    dbcan_utils CGC_substrate_abund -bt SRF_MT.depth.txt -i ../SRF_MG.dbCAN -a TPM 


    cd .. && cd BML_MG_abund 

    dbcan_utils fam_abund -bt BML_MG.depth.txt -i ../BML_MG.dbCAN -a TPM 

    dbcan_utils fam_substrate_abund -bt BML_MG.depth.txt -i ../BML_MG.dbCAN -a TPM 

    dbcan_utils CGC_abund -bt BML_MG.depth.txt -i ../BML_MG.dbCAN -a TPM 

    dbcan_utils CGC_substrate_abund -bt BML_MG.depth.txt -i ../BML_MG.dbCAN -a TPM 


    cd .. && cd SRF_MG_abund 

    dbcan_utils fam_abund -bt SRF_MG.depth.txt -i ../SRF_MG.dbCAN -a TPM 

    dbcan_utils fam_substrate_abund -bt SRF_MG.depth.txt -i ../SRF_MG.dbCAN -a TPM 

    dbcan_utils CGC_abund -bt SRF_MG.depth.txt -i ../SRF_MG.dbCAN -a TPM 

    dbcan_utils CGC_substrate_abund -bt SRF_MG.depth.txt -i ../SRF_MG.dbCAN -a TPM 

    cd .. 

We developed a set of Python scripts as ``dbcan_utils`` (included in the ``run_dbcan`` package) to take the raw read counts for all genes as input and output the normalized abundances (refer to Box 7) of CAZyme families, subfamilies, CGCs, and substrates (see Fig. 4). The parameter ``-a TPM`` can also be set to two other metrics: RPM, or RPKM61.

- **RPKM** is calculated as the number of mapped reads to a gene G divided by [(total number of mapped reads to all genes / 10^6) x (gene G length / 1000)].
- **RPM** is the number of mapped reads to a gene G divided by (total number of mapped reads to all genes / 10^6).
- **TPM** is calculated as [number of mapped reads to a gene G / (gene G length / 1000)] divided by the sum of [number of mapped reads to each gene / (the gene length / 1000)].

Box 7. Example output of dbcan_utils
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As an example, `SRF_MG_abund  <https://bcb.unl.edu/dbCAN_tutorial/dataset3-Priest2023/SRF_MG_abund/>_` folder has 7 TSV files:

.. code-block:: shell
    
    -rw-rw-r--  1 jinfang jinfang 122K Jun 19 07:30 CGC_abund.out 

    -rw-rw-r--  1 jinfang jinfang  835 Jun 19 07:30 CGC_substrate_majority_voting.out 

    -rw-rw-r--  1 jinfang jinfang 4.5K Jun 19 07:30 CGC_substrate_PUL_homology.out 

    -rw-rw-r--  1 jinfang jinfang 2.6K Jun 19 07:30 EC_abund.out 

    -rw-rw-r--  1 jinfang jinfang 3.3K Jun 19 07:30 fam_abund.out 

    -rw-rw-r--  1 jinfang jinfang  18K Jun 19 07:30 fam_substrate_abund.out 

    -rw-rw-r--  1 jinfang jinfang 6.9M Jun 19 04:07 SRF_MG.depth.txt 

    -rw-rw-r--  1 jinfang jinfang  19K Jun 19 07:30 subfam_abund.out 


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

P15. Heatmap for CAZyme substrate abundance across samples (TIMING 1min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    dbcan_plot heatmap_plot --show_abund --top 20 --samples BML_MG,BML_MT,SRF_MG,SRF_MT -i BML_MG_abund/fam_substrate_abund.out,BML_MT_abund/fam_substrate_abund.out,SRF_MG_abund/fam_substrate_abund.out,SRF_MT_abund/fam_substrate_abund.out 





P16. Barplot for CAZyme family/subfamily/EC abundance across samples (Fig. S5A-C) (TIMING 1min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell


    dbcan_plot bar_plot --samples BML_MG,BML_MT,SRF_MG,SRF_MT --vertical_bar --top 20 -i BML_MG_abund/fam_abund.out,BML_MT_abund/fam_abund.out,SRF_MG_abund/fam_abund.out,SRF_MT_abund/fam_abund.out --pdf fam.pdf 

    dbcan_plot bar_plot --samples BML_MG,BML_MT,SRF_MG,SRF_MT --vertical_bar --top 20 -i BML_MG_abund/subfam_abund.out,BML_MT_abund/subfam_abund.out,SRF_MG_abund/subfam_abund.out,SRF_MT_abund/subfam_abund.out --pdf subfam.pdf 

    dbcan_plot bar_plot --samples BML_MG,BML_MT,SRF_MG,SRF_MT --vertical_bar --top 20 -i BML_MG_abund/EC_abund.out,BML_MT_abund/EC_abund.out,SRF_MG_abund/EC_abund.out,SRF_MT_abund/EC_abund.out --pdf ec.pdf 


P17.  Synteny plot between a CGC and its best PUL hit with read mapping coverage to CGC(use --db_dir if needed) (Fig. S5E) (TIMING 1min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    dbcan_plot CGC_synteny_coverage_plot -i SRF_MG.dbCAN --cgcid 'contig_12725|CGC1' --readscount SRF_MG_abund/SRF_MG.cgc.depth.txt 

    dbcan_plot CGC_synteny_coverage_plot -i SRF_MG.dbCAN --cgcid 'contig_12725|CGC1' --readscount SRF_MT_abund/SRF_MT.cgc.depth.txt 


If users only want to plot the CGC structure:

.. code-block:: shell

    dbcan_plot CGC_plot -i SRF_MG.dbCAN --cgcid 'contig_12725|CGC1' 

If users only want to plot the CGC structure plus the read mapping coverage:

.. code-block:: shell

    dbcan_plot CGC_coverage_plot -i SRF_MG.dbCAN --cgcid 'contig_12725|CGC1' --readscount SRF_MG_abund/SRF_MG.cgc.depth.txt 

If users only want to plot the synteny between the CGC and PUL:

.. code-block:: shell

    dbcan_plot CGC_synteny_plot -i SRF_MG.dbCAN --cgcid 'contig_12725|CGC1' 




