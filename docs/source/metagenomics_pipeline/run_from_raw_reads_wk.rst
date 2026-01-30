Run from Raw Reads(Wastyk 2021): Automated CAZyme and Glycan Substrate Annotation in Microbiomes: A Step-by-Step Protocol
===========================================================================================================================

.. _wastyk_2021: https://www.cell.com/cell/fulltext/S0092-8674(21)00754-6

Example 2: Wastyk2021 Dataset `wastyk_2021`_
-------------------------------------------------

The Wastyk2021 dataset `wastyk_2021`_ was published in 2021 from a human dietary intervention study. In the published paper, researchers studied how high-fermented and high-fiber diets influence the human microbiome metabolism and modulate the human immune status. Among various data analyses conducted in the paper1, CAZymes were mined from shotgun metagenomic reads of 18 healthy human participants, and each participant had four time points of stool samples for metagenome sequencing. CAZyme abundance profiles were compared before and after the high-fiber intervention (baseline vs high-fiber). One of the main findings from their CAZyme analysis was that high-fiber consumption increased the CAZyme abundance. For this protocol, we will select two samples (paired-end 2x146bp reads) of two time points (day 2 before high-fiber diet as baseline, and 10 weeks after high-fiber diet as intervention) from one participant. The protocol is for the individual sample route.

The raw read data, intermediate data from each analysis step, and final result data and visualization files are organized in nested folders available on our website https://bcb.unl.edu/dbCAN_tutorial/dataset2-Wastyk2021/) and https://dbcan.readthedocs.io.

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

    conda env create -f CAZyme_annotation_wk.yml
    conda activate CAZyme_annotation_wk

To install the databases, execute the following commands:

.. include:: prepare_the_database.rst



Module 1: Reads processing to obtain contigs
`````````````````````````````````````````````````````

P1. Contamination Check (TIMING ~10min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Download the Wastyk2021 dataset:

.. code-block:: shell

    wget https://bcb.unl.edu/dbCAN_tutorial/dataset2-Wastyk2021/fefifo_8022_1__shotgun_1.fastq.gz
    wget https://bcb.unl.edu/dbCAN_tutorial/dataset2-Wastyk2021/fefifo_8022_1__shotgun_2.fastq.gz
    wget https://bcb.unl.edu/dbCAN_tutorial/dataset2-Wastyk2021/fefifo_8022_7__shotgun_1.fastq.gz
    wget https://bcb.unl.edu/dbCAN_tutorial/dataset2-Wastyk2021/fefifo_8022_7__shotgun_2.fastq.gz

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


Use `kraken2` to check for contaminated reads:

.. code-block:: shell

    kraken2 --threads 32 --quick --paired --db K2 --report fefifo_8022_1.kreport --output fefifo_8022_1.kraken.output fefifo_8022_1__shotgun_1.fastq.gz fefifo_8022_1__shotgun_2.fastq.gz
    kraken2 --threads 32 --quick --paired --db K2 --report fefifo_8022_7.kreport --output fefifo_8022_7.kraken.output fefifo_8022_7__shotgun_1.fastq.gz fefifo_8022_7__shotgun_2.fastq.gz

Kraken2 found very little contamination in the data. Consequently, there was no need for the contamination removal step.

If contamination is identified, users can align the reads to the reference genomes of potential contamination source organisms to remove
the aligned reads (Box 1). The most common source in human microbiome studies is from human hosts.

Box 1: Example to remove contamination reads from human
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Kraken2 will produce the following output files:

    .. code-block:: shell

        -rw-rw-r-- 1 jinfang jinfang 1.1G Sep 21 23:17 fefifo_8022_1.kraken.output
        -rw-rw-r-- 1 jinfang jinfang 991K Sep 21 23:19 fefifo_8022_1.kreport
        -rw-rw-r-- 1 jinfang jinfang 574M Sep 21 23:21 fefifo_8022_7.kraken.output
        -rw-rw-r-- 1 jinfang jinfang 949K Sep 21 23:22 fefifo_8022_7.kreport


    Suppose from these files, we have identified humans as the contamination source, we can use the following commands to remove the contamination reads by aligning reads to the human reference genome.

    .. code-block:: shell

        wget https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
        bwa index -p hg38 Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
        bwa mem hg38 fefifo_8022_1__shotgun_1.fastq.gz fefifo_8022_1__shotgun_2.fastq.gz -t 32 -o fefifo_8022_1.hg38.sam
        bwa mem hg38 fefifo_8022_7__shotgun_1.fastq.gz fefifo_8022_7__shotgun_2.fastq.gz -t 32 -o fefifo_8022_7.hg38.sam
        samtools view -f 12 fefifo_8022_1.hg38.sam > fefifo_8022_1.hg38.unmap.bam
        samtools view -f 12 fefifo_8022_7.hg38.sam > fefifo_8022_7.hg38.unmap.bam
        samtools fastq -1 fefifo_8022_1_1.clean.fq.gz -2 fefifo_8022_1_2.clean.fq.gz fefifo_8022_1.hg38.unmap.bam
        samtools fastq -1 fefifo_8022_7_1.clean.fq.gz -2 fefifo_8022_7_2.clean.fq.gz fefifo_8022_7.hg38.unmap.bam


KrakenTools could also extract host reads quickly and easied which is recommended. We use tax 2759 (plant) as an example.
Please read KrakenTools README for more information (https://github.com/jenniferlu717/KrakenTools?tab=readme-ov-file).

.. code-block:: shell

    extract_kraken_reads.py \
    -k fefifo_8022_1.kraken.output \
    -s1 fefifo_8022_1__shotgun_1.fastq.gz -s2 fefifo_8022_1__shotgun_2.fastq.gz \
    --fastq-output --exclude \
    --taxid 9606 \
    -o fefifo_8022_1_1.clean.fq.gz -o2 fefifo_8022_1_2.clean.fq.gz

    extract_kraken_reads.py \
    -k fefifo_8022_7.kraken.output \
    -s1 fefifo_8022_7__shotgun_1.fastq.gz -s2 fefifo_8022_7__shotgun_2.fastq.gz \
    --fastq-output --exclude \
    --taxid 9606 \
    -o fefifo_8022_7_1.clean.fq.gz -o2 fefifo_8022_7_2.clean.fq.gz


P2. Trim adapter and low-quality reads (TIMING ~20min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    trim_galore --paired fefifo_8022_1__shotgun_1.fastq.gz fefifo_8022_1__shotgun_2.fastq.gz --illumina -j 36
    trim_galore --paired fefifo_8022_7__shotgun_1.fastq.gz fefifo_8022_7__shotgun_2.fastq.gz --illumina -j 36

We specified --illumina to indicate that the reads were generated using the Illumina sequencing platform.
Nonetheless, trim_galore can automatically detect adapters, providing flexibility for users who may know the specific sequencing platform.
Details of trimming are available in the trimming report file (Box 2).

Box 2: Example output of `trim_galore`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    In addition to the trimmed read files, `Trim_galore`` also generates a trimming report file.
    The trimming report contains details on read trimming, such as the number of trimmed reads.

.. code-block:: shell

    -rw-rw-r-- 1 jinfang jinfang  429M Oct 30 22:44 fefifo_8022_1__shotgun_1.fastq.gz
    -rw-rw-r-- 1 jinfang jinfang  4.1K Oct 31 05:15 fefifo_8022_1__shotgun_1.fastq.gz_trimming_report.txt
    -rw-rw-r-- 1 jinfang jinfang  390M Oct 31 05:16 fefifo_8022_1__shotgun_1_val_1.fq.gz
    -rw-rw-r-- 1 jinfang jinfang  540M Oct 30 22:44 fefifo_8022_1__shotgun_2.fastq.gz
    -rw-rw-r-- 1 jinfang jinfang  4.2K Oct 31 05:16 fefifo_8022_1__shotgun_2.fastq.gz_trimming_report.txt
    -rw-rw-r-- 1 jinfang jinfang  499M Oct 31 05:16 fefifo_8022_1__shotgun_2_val_2.fq.gz
    -rw-rw-r-- 1 jinfang jinfang  931M Oct 30 22:34 fefifo_8022_7__shotgun_1.fastq.gz
    -rw-rw-r-- 1 jinfang jinfang  4.2K Oct 31 05:17 fefifo_8022_7__shotgun_1.fastq.gz_trimming_report.txt
    -rw-rw-r-- 1 jinfang jinfang  861M Oct 31 05:20 fefifo_8022_7__shotgun_1_val_1.fq.gz
    -rw-rw-r-- 1 jinfang jinfang  1.1G Oct 30 22:34 fefifo_8022_7__shotgun_2.fastq.gz
    -rw-rw-r-- 1 jinfang jinfang  4.4K Oct 31 05:20 fefifo_8022_7__shotgun_2.fastq.gz_trimming_report.txt
    -rw-rw-r-- 1 jinfang jinfang 1003M Oct 31 05:20 fefifo_8022_7__shotgun_2_val_2.fq.gz

.. warning::

    During the trimming process, certain reads may be entirely removed due to low quality in its entirety.
    Using the ``--retain_unpaired`` parameter in ``trim_galore`` allows for the preservation of single-end reads.
    In this protocol, this option was not selected, so that both reads of a forward-revise pair were removed.

P3. Assemble reads into contigs (TIMING ~84min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Use Megahit for assembling reads into contigs:

.. code-block:: shell

    megahit -m 0.5 -t 32 -o megahit_fefifo_8022_1 -1 fefifo_8022_1__shotgun_1_val_1.fq.gz -2 fefifo_8022_1__shotgun_2_val_2.fq.gz --out-prefix fefifo_8022_1 --min-contig-len 1000
    megahit -m 0.5 -t 32 -o megahit_fefifo_8022_7 -1 fefifo_8022_7__shotgun_1_val_1.fq.gz -2 fefifo_8022_7__shotgun_2_val_2.fq.gz --out-prefix fefifo_8022_7 --min-contig-len 1000


`MEGAHIT` generates two output folders `megahit_fefifo_8022_1` and `megahit_fefifo_8022_7`.
Each contains five files and one sub-folder (Box 3). `fefifo_8022_1.contigs.fa` is the final contig sequence file.
We set `--min-contig-len 1000`, a common practice to retain all contigs longer than 1,000 base pairs.

Box 3: Example output of `MEGAHIT`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: shell

    -rw-rw-r--  1 jinfang jinfang  262 Oct 31 05:49 checkpoints.txt
    -rw-rw-r--  1 jinfang jinfang    0 Oct 31 05:49 done
    -rw-rw-r--  1 jinfang jinfang  97M Oct 31 05:49 fefifo_8022_1.contigs.fa
    -rw-rw-r--  1 jinfang jinfang 149K Oct 31 05:49 fefifo_8022_1.log
    drwxrwxr-x  2 jinfang jinfang 4.0K Oct 31 05:49 intermediate_contigs
    -rw-rw-r--  1 jinfang jinfang 1.1K Oct 31 05:27 options.json


P4. Predict genes by `Pyrodigal` (TIMING ~24min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    pyrodigal -q -p meta -i megahit_fefifo_8022_1/fefifo_8022_1.contigs.fa -d fefifo_8022_1.cds -a fefifo_8022_1.faa -f gff -o fefifo_8022_1.gff -j 36
    pyrodigal -q -p meta -i megahit_fefifo_8022_1/fefifo_8022_7.contigs.fa -d fefifo_8022_7.cds -a fefifo_8022_7.faa -f gff -o fefifo_8022_7.gff -j 36


Box 4: Example output of `Pyrodigal`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: shell

        -rw-rw-r--  1 jinfang jinfang  31M Oct 31 22:09 fefifo_8022_1.faa
        -rw-rw-r--  1 jinfang jinfang  98M Oct 31 21:50 fefifo_8022_1.cds
        -rw-rw-r--  1 jinfang jinfang 142M Oct 31 22:09 fefifo_8022_1.gff
        -rw-rw-r--  1 jinfang jinfang  31M Oct 31 22:09 fefifo_8022_7.faa
        -rw-rw-r--  1 jinfang jinfang  98M Oct 31 21:50 fefifo_8022_7.cds
        -rw-rw-r--  1 jinfang jinfang 142M Oct 31 22:09 fefifo_8022_7.gff

Hint: If you apply ``prodigal`` not ``pyrodigal``, use this to fix gene ID in gff files by dbcan_utils (TIMING ~1min)

.. code-block:: shell

    dbcan_utils gff_fix -i fefifo_8022_1.faa -g fefifo_8022_1.gff
    dbcan_utils gff_fix -i fefifo_8022_7.faa -g fefifo_8022_1.gff


Module 2. run_dbcan annotation (Fig. 3) to obtain CAZymes, CGCs, and substrates
```````````````````````````````````````````````````````````````````````````````

**CRITICAL STEP**

Users can skip P5 and P6, and directly run P7 (much slower though), if they want to predict not only CAZymes and CGCs, but also substrates.

P5. CAZyme annotation at the CAZyme family level (TIMING ~10min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    run_dbcan easy_CAZyme --input_raw_data fefifo_8022_1.faa --mode protein --output_dir fefifo_8022_1.CAZyme --db_dir db --input_format NCBI
    run_dbcan easy_CAZyme --input_raw_data fefifo_8022_7.faa --mode protein --output_dir fefifo_8022_7.CAZyme --db_dir db --input_format NCBI

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

    run_dbcan easy_CGC  --input_raw_data fefifo_8022_1.faa --mode protein --output_dir fefifo_8022_1.PUL --input_format NCBI --gff_type prodigal --input_gff fefifo_8022_1.gff --db_dir db
    run_dbcan easy_CGC  --input_raw_data fefifo_8022_7.faa --mode protein --output_dir fefifo_8022_7.PUL --input_format NCBI --gff_type prodigal --input_gff fefifo_8022_7.gff --db_dir db

.. warning::

    **Creating own gff file**
    If the users would like to create their own ``gff`` file (instead of using Prokka or Prodigal),
    it is important to make sure the value of ID attribute in the ``gff`` file matches the protein ID in the protein ``faa`` file.

    **[Troubleshooting]CGC not found**
    If no result is found in CGC output file, it is most likely because the sequence IDs in ``gff`` file and ``faa`` file do not match.
    Another less likely reason is that the contigs are too short and fragmented and not suitable for CGC prediction.


P7. Substrate prediction for CAZymes and CGCs (TIMING ~5h)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following commands will re-run run_dbcan to predict CAZymes, CGCs, and their substrates with the `--cgc_substrate` parameter.

.. code-block:: shell

    run_dbcan  easy_substrate --input_raw_data  fefifo_8022_1.faa --mode protein  --input_format NCBI --input_gff fefifo_8022_1.gff --gff_type prodigal   --output_dir fefifo_8022_1.dbCAN --db_dir db
    run_dbcan  easy_substrate --input_raw_data  fefifo_8022_7.faa --mode protein  --input_format NCBI --input_gff fefifo_8022_7.gff --gff_type prodigal   --output_dir fefifo_8022_7.dbCAN --db_dir db



In the `fefifo_8022_1.dbCAN <https://bcb.unl.edu/dbCAN_tutorial/dataset2-Wastyk2021/fefifo_8022_1.dbCAN/>_` directory, a total of 17 files and 1 folder are generated:



    Descriptions of Output Files:

    - ``PUL_blast.out``: BLAST results between CGCs and PULs.
    - ``CGC.faa``: CGC Fasta sequences.
    - ``cgc.gff``: reformatted from the user input gff file by marking CAZymes, TFs, TCs, and STPs.
    - ``cgc.out``: raw output of CGC predictions.

        1.	CGC_id: CGC1
        2.	type: CAZyme
        3.	contig_id: k141_32617
        4.	gene_id: fefifo_8022_1_00137
        5.	start: 1755
        6.	end: 3332
        7.	strand: -
        8.	annotation: GH13

    **Explanation**: Explanation: the gene fefifo_8022_1_00137 encodes a GH13 CAZyme in the CGC1 of the contig k141_32617. CGC1 also has other genes, which are provided in other rows. fefifo_8022_1_00137 is on the negative strand of k141_32617 from 1755 to 3332. The type can be one of the four signature gene types (CAZymes, TCs, TFs, STPs) or the null type (not annotated as one of the four signature genes).

    - ``cgc_standard.out.json``: JSON format of cgc_standard.out.
    - ``dbcan-sub.hmm.out``: HMMER search result against dbCAN-sub HMMdb, including a column with CAZyme substrates extracted from `fam-substrate-mapping-08012023.tsv`.
    - ``diamond.out``: DIAMOND search result against the CAZy annotated protein sequences (`CAZyDB.07262023.fa`).
    - ``dtemp.out``: temporary file.
    - ``hmmer.out``: HMMER search result against dbCAN HMMdb.
    - ``overview.txt``: summary of CAZyme annotation from three methods in TSV format. An example row has the following columns:

        1. ``Gene_ID``: fefifo_8022_1_00719
        2. ``EC#``: PL8_e13:2
        3. ``dbCAN``: PL8_2(368-612)
        4. ``dbCAN_sub``: PL8_e13
        5. ``DIAMOND``: PL8_2
        6. ``#ofTools``: 3

    **Explanation**: Explanation: the protein fefifo_8022_1_00719 is annotated by 3 tools to be a CAZyme: (1) PL8_2 (CAZy defined subfamily 2 of PL8) by HMMER vs dbCAN HMMdb with a domain range from aa position 368 to 612, (2) PL8_e13 (eCAMI defined subfamily e13; e indicates it is from eCAMI not CAZy) by HMMER vs dbCAN-sub HMMdb (derived from eCAMI subfamilies), and (3) PL8_2 by DIAMOND vs CAZy annotated protein sequences. The second column 4.2.2.20:2 is extracted from eCAMI, meaning that the eCAMI subfamily PL8_e13 contains two member proteins which have an EC 4.2.2.20 according to CAZy. In most cases, the 3 tools will have the same CAZyme family assignment. When they give different assignment. We recommend a preference order: dbCAN > eCAMI/dbCAN-sub > DIAMOND. See our dbCAN2 paper2, dbCAN3 paper3, and eCAMI4 for more details.

    - ``stp.out``: HMMER search result against the MiST5 compiled signal transduction protein HMMs from Pfam.
    - ``tf-1.out``: HMMER search result against the DBD6 compiled transcription factor HMMs from Pfam 7.
    - ``tf-2.out``: HMMER search result against the DBD compiled transcription factor HMMs from Superfamily 8.
    - ``tp.out``: DIAMOND search result against the TCDB 9 annotated protein sequences.
    - ``substrate.out``: summary of substrate prediction results for CGCs in TSV format from two approaches3 (dbCAN-PUL blast search and dbCAN-sub majority voting). An example row has the following columns:

        1. ``CGC_ID``: k141_31366|CGC2
        2. ``Best hit PUL_ID in dbCAN-PUL``: PUL0008
        3. ``Substrate of the hit PUL``: fructan
        4. ``Sum of bitscores for homologous gene pairs between CGC and PUL``: 6132.0
        5. ``Types of homologous gene pairs``: CAZyme-CAZyme;CAZyme-CAZyme;TC-TC;CAZyme-CAZyme;CAZyme-CAZyme;TC-TC
        6. ``Substrate predicted by majority voting of CAZymes in CGC``: fructan
        7. ``Voting score``: 2.0

    **Explanation**: The CGC1 of contig ``k141_31366`` has its best hit ``PUL0008`` (from ``PUL_blast.out``) with fructan as substrate (from ``dbCAN-PUL_12-12-2023.xlsx``). Six signature genes are matched between ``k141_31366|CGC2 and PUL0008 (from PUL_blast.out)``: four are CAZymes and the other two are TCs. The sum of blast bitscores of the six homologous pairs (``CAZyme-CAZyme, CAZyme-CAZyme, TC-TC, CAZyme-CAZyme, CAZyme-CAZyme and TC-TC``) is 6132.0. Hence, the substrate of ``k141_31366|CGC2`` is predicted to be fructan according to dbCAN-PUL blast search. The last two columns are based on the dbCAN-sub result (``dbcan-sub.hmm.out``), according to which two CAZymes in ``k141_31366|CGC2`` are predicted to have fructan substrate. The voting score is thus 2.0, so that according to the majority voting rule, ``k141_31366|CGC2`` is predicted to have a fructan substrate.

    **Note**: : for many CGCs, only one of the two approaches produces substrate prediction. In some cases, the two approaches produce different substrate assignments. We recommend a preference order: ``dbCAN-PUL blast search > dbCAN-sub`` majority voting. See our `dbCAN3 <https://academic.oup.com/nar/article/51/W1/W115/7147496>_` paper3 for more details.

    - ``synteny.pdf``: a folder with syntenic block alignment plots between all CGCs and PULs.
    - ``uniInput``: renamed Fasta file from input protein sequence file.


Module 3. Read mapping (Fig. 3) to calculate abundance for CAZyme families, subfamilies, CGCs, and substrates
``````````````````````````````````````````````````````````````````````````````````````````````````````````````

P8. Read mapping to all contigs of each sample (TIMING ~10 min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    mkdir samfiles 
    bwa index megahit_fefifo_8022_1/fefifo_8022_1.contigs.fa 
    bwa index megahit_fefifo_8022_7/fefifo_8022_7.contigs.fa 
    bwa mem -t 32 -o samfiles/fefifo_8022_1.sam megahit_fefifo_8022_1/fefifo_8022_1.contigs.fa fefifo_8022_1__shotgun_1_val_1.fq.gz fefifo_8022_1__shotgun_2_val_2.fq.gz 
    bwa mem -t 32 -o samfiles/fefifo_8022_7.sam megahit_fefifo_8022_7/fefifo_8022_7.contigs.fa fefifo_8022_7__shotgun_1_val_1.fq.gz fefifo_8022_7__shotgun_2_val_2.fq.gz 



P9. Sort SAM files by coordinates (TIMING ~6min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    cd samfiles
    samtools sort -@ 32 -o fefifo_8022_1.bam fefifo_8022_1.sam
    samtools sort -@ 32 -o fefifo_8022_7.bam fefifo_8022_7.sam
    samtools index fefifo_8022_1.bam
    samtools index fefifo_8022_7.bam
    rm -rf *sam
    cd ..


P10. Read count calculation for all proteins of each sample using dbcan_utils (TIMING ~1min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    mkdir fefifo_8022_1_abund && cd fefifo_8022_1_abund
    dbcan_utils cal_coverage -g ../fefifo_8022_1.fix.gff -i ../samfiles/fefifo_8022_1.bam -o fefifo_8022_1.depth.txt -t 6 --overlap_base_ratio 0.2 --mapping_quality 30 --identity 0.98  

    cd .. && mkdir fefifo_8022_7_abund && cd fefifo_8022_7_abund 
    dbcan_utils cal_coverage -g ../fefifo_8022_7.fix.gff -i ../samfiles/fefifo_8022_7.bam -o fefifo_8022_7.depth.txt -t 6 --overlap_base_ratio 0.2 --mapping_quality 30 --identity 0.98
    cd ..


P11. Read count calculation for a given region of contigs using Samtools (TIMING ~1min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    cd fefifo_8022_1_abund
    samtools index ../samfiles/fefifo_8022_1.bam
    samtools depth -r k141_2168:4235-19858 ../samfiles/fefifo_8022_1.bam > fefifo_8022_1.cgc.depth.txt
    cd ..


The parameter ``-r k141_2168:4235-19858`` specifies a region in a contig. For any CGC, its positional range can be found in the file ``cgc_standard.out`` produced by run_dbcan (Box 6). The ``depth.txt`` files contain the raw read counts for the specified region.


.. warning::

    The contig IDs are automatically generated by MEGAHIT. There is a small chance that a same contig ID appears in both samples. However, the two contigs in the two samples do not match each other even the ID is the same. For example, the contig ID ``k141_2168`` is most likely only found in the ``fefifo_8022_1`` sample. Even if there is a ``k141_2168`` in ``fefifo_8022_7``, the actual contigs in two samples are different.

P12. dbcan_utils to calculate the abundance of CAZyme families, subfamilies, CGCs, and substrates (TIMING ~1min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell
    
    cd fefifo_8022_1_abund
    dbcan_utils fam_abund -bt fefifo_8022_1.depth.txt -i ../fefifo_8022_1.dbCAN -a TPM
    dbcan_utils fam_substrate_abund -bt fefifo_8022_1.depth.txt -i ../fefifo_8022_1.dbCAN -a TPM
    dbcan_utils CGC_abund -bt fefifo_8022_1.depth.txt -i ../fefifo_8022_1.dbCAN -a TPM
    dbcan_utils CGC_substrate_abund -bt fefifo_8022_1.depth.txt -i ../fefifo_8022_1.dbCAN -a TPM

    cd .. && cd fefifo_8022_7_abund
    dbcan_utils fam_abund -bt fefifo_8022_7.depth.txt -i ../fefifo_8022_7.dbCAN -a TPM
    dbcan_utils fam_substrate_abund -bt fefifo_8022_7.depth.txt -i ../fefifo_8022_7.dbCAN -a TPM
    dbcan_utils CGC_abund -bt fefifo_8022_7.depth.txt -i ../fefifo_8022_7.dbCAN -a TPM
    dbcan_utils CGC_substrate_abund -bt fefifo_8022_7.depth.txt -i ../fefifo_8022_7.dbCAN -a TPM
    cd ..


We developed a set of Python scripts as ``dbcan_utils`` (included in the ``run_dbcan`` package) to take the raw read counts for all genes as input and output the normalized abundances (refer to Box 7) of CAZyme families, subfamilies, CGCs, and substrates (see Fig. 4). The parameter ``-a TPM`` can also be set to two other metrics: RPM, or RPKM61.

- **RPKM** is calculated as the number of mapped reads to a gene G divided by [(total number of mapped reads to all genes / 10^6) x (gene G length / 1000)].
- **RPM** is the number of mapped reads to a gene G divided by (total number of mapped reads to all genes / 10^6).
- **TPM** is calculated as [number of mapped reads to a gene G / (gene G length / 1000)] divided by the sum of [number of mapped reads to each gene / (the gene length / 1000)].


Box 7. Example output of dbcan_utils
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As an example, `fefifo_8022_1_abund <https://bcb.unl.edu/dbCAN_tutorial/dataset2-Wastyk2021/fefifo_8022_1_abund/>_` folder has 7 TSV files:

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

1. two abundance folders ``fefifo_8022_1_abund`` and ``fefifo_8022_7_abund``,
2. two CAZyme annotation ``folders fefifo_8022_1.dbCAN`` and ``fefifo_8022_7.dbCAN``, and
3. the ``dbCAN-PUL folder`` (under the db folder, released from ``dbCAN-PUL.tar.gz``).

P13. Heatmap for CAZyme substrate abundance across samples (Fig. S4B) (TIMING 1min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    dbcan_plot heatmap_plot --samples fefifo_8022_1,fefifo_8022_7 -i fefifo_8022_1_abund/ fam_substrate_abund.out,fefifo_8022_7_abund/fam_substrate_abund.out --show_abund --top 20


Here we plot the top 20 substrates in the two samples. The input files are the two CAZyme substrate abundance files calculated based on dbCAN-sub result. The default heatmap is ranked by substrate abundances. To rank the heatmap according to abundance profile using the function clustermap of seaborn package, users can invoke the ``--cluster_map`` parameter.

P14. Barplot for CAZyme family/subfamily/EC abundance across samples (Fig. S4C) (TIMING 1min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    dbcan_plot bar_plot --samples fefifo_8022_1,fefifo_8022_7 --vertical_bar --top 20 -i fefifo_8022_1_abund/fam_abund.out,fefifo_8022_7_abund/fam_abund.out --pdf fam.pdf --db_dir db

    dbcan_plot bar_plot --samples fefifo_8022_1,fefifo_8022_7 --vertical_bar --top 20 -i fefifo_8022_1_abund/subfam_abund.out,fefifo_8022_7_abund/subfam_abund.out --pdf subfam.pdf --db_dir db

    dbcan_plot bar_plot --samples fefifo_8022_1,fefifo_8022_7 --vertical_bar --top 20 -i fefifo_8022_1_abund/EC_abund.out,fefifo_8022_7_abund/EC_abund.out --pdf ec.pdf --db_dir db

Users can choose to generate a barplot instead of heatmap using the ``bar_plot`` method.

P15. Synteny plot between a CGC and its best PUL hit with read mapping coverage to CGC (Fig. S4A) (TIMING 1min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    dbcan_plot CGC_synteny_coverage_plot -i fefifo_8022_1.dbCAN --cgcid 'k141_2168|CGC1' --readscount fefifo_8022_1_abund/fefifo_8022_1.cgc.depth.txt


The ``fefifo_8022_1.dbCAN`` folder contains the ``PUL_blast.out`` file. Using this file, the ``cgc_standard.out`` file,
and the best PUL's ``gff`` file in ``dbCAN-PUL.tar.gz``, the CGC_synteny_plot method will create the ``CGC-PUL synteny plot``.
The ``-cgcid`` parameter is required to specify which CGC to be plotted (``'k141_2168|CGC1'`` in this example).
The ``fefifo_8022_1.cgc.depth.txt`` file is used to plot the read mapping coverage.

If users only want to plot the CGC structure:

.. code-block:: shell

    dbcan_plot CGC_plot -i fefifo_8022_1.dbCAN --cgcid 'k141_2168|CGC1'

If users only want to plot the CGC structure plus the read mapping coverage:

.. code-block:: shell

    dbcan_plot CGC_coverage_plot -i fefifo_8022_1.dbCAN --cgcid 'k141_2168|CGC1' --readscount fefifo_8022_1_abund/fefifo_8022_1.cgc.depth.txt

If users only want to plot the synteny between the CGC and PUL:

.. code-block:: shell

    dbcan_plot CGC_synteny_plot -i fefifo_8022_1.dbCAN --cgcid 'k141_2168|CGC1'


.. warning::

    The CGC IDs in different samples do not match each other. For example, specifying ``-i fefifo_8022_1.dbCAN`` is to plot
    the ``'k141_2168|CGC1'`` in the fefifo_8022_1 sample. The ``'k141_2168|CGC1'`` in the fefifo_8022_7 sample most likely does not exist,
    and even it does, the CGC has a different sequence even if the ID is the same.
