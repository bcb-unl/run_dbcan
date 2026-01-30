Run from Raw Reads(Cater 2023): Automated CAZyme and Glycan Substrate Annotation in Microbiomes: A Step-by-Step Protocol
========================================================================================================================

Introduction
------------

Overview
````````

In this tutorial, we present a comprehensive protocol to annotate CAZymes and glycan substrates in microbiome datasets. Using three real-world microbiome example datasets:

this guide will walk you through each step of the computational workflow for analyzing occurrence and abundance. The workflow, depicted in Fig. 1, is designed to be user-friendly and does not require extensive programming knowledge.


Workflow Steps
``````````````

1. **Pre-Processing of Raw Sequencing Reads:**
    Begin with the preprocessing of raw sequencing reads. This includes the removal of contaminants, adapter sequences, and trimming of low-quality reads. We'll use `trim_galore` and `Kraken2` for this purpose (Steps 1-2).

2. **Contig Assembly:**
    The clean reads from each sample are then assembled into contigs using `MEGAHIT` (Step 3).

3. **Gene Model Annotation:**
    These contigs are subsequently passed to `pyrodigal` for gene model annotation (Step 4).

4. **CAZyme and CGC Annotation:**
    The next phase involves annotating the contigs for CAZymes and CGCs. This is achieved by `run_dbcan`, utilizing the protein sequence (faa) and gene annotation (gff) files produced by `pyrodigal` (Step 5).

5. **Location Mapping and Substrate Prediction:**
    Step 6 involves mapping the location of annotated CAZymes and CGCs on the contigs. In Step 7, `run_dbcan`'s substrate prediction function infers glycan substrates for these CAZymes and CGCs.

6. **Abundance Calculation:**
    To quantify the abundance of CAZymes, substrates, and CGCs, clean reads from Step 2 are mapped to the nucleotide coding sequences (CDS) of proteins from Step 4 (Steps 8-14).

7. **Data Visualization:**
    Finally, this step focuses on visualizing the occurrence and abundance results. We provide Python scripts for creating publication-quality plots in PDF format.


.. image:: ../_static/img/Picture1.png
    :alt: workflow figure
    :width: 800px
    :align: center

.. |centered-text2| raw:: html

    <div style="text-align: center">Fig.1 <strong>Experimental design of CAZyme annotation in microbiomes</strong></div>

|centered-text2|

User Requirements
`````````````````

This protocol is designed for users who are comfortable with the Linux command-line interface and can execute Python scripts in the terminal. While extensive programming experience is not necessary, users should be familiar with editing Linux commands and plain-text scripts within a command-line environment.

Equipment
---------

Operating System
````````````````

All the modules of this protocol (Fig. 1) are designed to run on a command line (CLI) environment with a Linux OS (e.g., Ubuntu).
We recommend users install these modules and execute all commands on a high-performance Linux cluster or workstation with >32 CPUs
and 128GB of RAM instead of a laptop, as the assembly of raw reads has a high demand of CPU and RAM.

Once users finish the data visualization module (Fig. 1), the resulting image files (PDF format) can be copied
to a desktop or laptop with GUI for data visualization. In practice, users can choose not to use our read processing
module and read mapping module. They may instead use their preferred tools for preparing input data for run_dbcan module
and for calculating abundance for CAZymes and substrates. In that case, they can skip the installation of our read processing
module and read mapping module in this protocol.

Data Files
``````````

The example dataset ``Carter 2023`` is described above.
The raw read data, intermediate data from each analysis step, and final result
data and visualization files are organized in nested folders available on our
website (https://bcb.unl.edu/dbCAN_tutorial/dataset1-Carter2023/) and
https://dbcan.readthedocs.io. These websites also include data files and
protocols for four additional example datasets from ``Wastyk 2021`` , ``Priest 2023``, ``Amelia 2024``, and ``Emilson 2024``
which are not included in this protocol paper. We will use the independent sample
assembly route for ``Carter 2023`` in the main text to demonstrate all the commands.
Commands for the other routes are provided Supplementary Protocols.

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


Anaconda as the Software Management System
``````````````````````````````````````````
Anaconda will be used as the software package management system for this
protocol. Anaconda uses the ``conda`` command to create a virtual
environment to facilitate the easy installation of software packages
and running command line jobs. With the conda environment, users do
not need to worry about the potential issues of package dependencies
and version conflicts.

Like in all bioinformatics data analysis tasks, we recommend users organize
their data files by creating a dedicated folder for each data analysis
step.

.. _cater_2023: https://www.sciencedirect.com/science/article/pii/S0092867423005974

Example 1: Carter 2023 Dataset ``Carter 2023``
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

    The dbCAN.hmm and dbCAN_sub.hmm files are for HMMER search.

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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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


P3. Assemble reads into contigs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Use Megahit for assembling reads into contigs:


.. code-block:: shell

    megahit -m 0.5 -t 32 -o megahit_Wet2014 -1 Wet2014_1_val_1.fq.gz -2 Wet2014_2_val_2.fq.gz --out-prefix Wet2014 --min-contig-len 1000
    megahit -m 0.5 -t 32 -o megahit_Dry2014 -1 Dry2014_1_val_1.fq.gz -2 Dry2014_2_val_2.fq.gz --out-prefix Dry2014 --min-contig-len 1000


``MEGAHIT`` generates two output folders. Each contains five files and one sub-folder (Box 3).
``Wet2014.contigs.fa`` is the final contig sequence file. We set `--min-contig-len 1000`,
a common practice to retain all contigs longer than 1,000 base pairs.

Box 3: Example output of `MEGAHIT`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: shell

    -rw-rw-r--  1 jinfang jinfang  262 Dec 13 04:19 checkpoints.txt
    -rw-rw-r--  1 jinfang jinfang    0 Dec 13 04:19 done
    drwxrwxr-x  2 jinfang jinfang 4.0K Dec 13 04:19 intermediate_contigs
    -rw-rw-r--  1 jinfang jinfang 1.1K Dec 13 02:22 options.json
    -rw-rw-r--  1 jinfang jinfang 258M Dec 13 04:19 Wet2014.contigs.fa
    -rw-rw-r--  1 jinfang jinfang 208K Dec 13 04:19 Wet2014.log

.. warning::

    A common practice in metagenomics after assembly is to further bin contigs into metagenome-assembled genomes (MAGs).
    However, in this protocol, we chose not to generate MAGs because not all contigs can be binned into MAGs, and those un-binned
    contigs can also encode CAZymes.


P4. Predict genes by `Prodigal` (TIMING ~21min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    pyrodigal -p meta -i megahit_Wet2014/Wet2014.contigs.fa -d Wet2014.cds -a Wet2014.faa -f gff -o Wet2014.gff -j 36
    pyrodigal -p meta -i megahit_Dry2014/Dry2014.contigs.fa -d Dry2014.cds -a Dry2014.faa -f gff -o Dry2014.gff -j 36


The output files comprise of both protein and CDS sequences in Fasta format (e.g., Dry2014.faa and Wet2014.cds in Box 4).
The parameter -o and -f gff indicates the general feature format of gene annotation (e.g., Dry2014.gff).


Box 4: Example output of `Prodigal`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: shell

    -rw-rw-r--  1 jinfang jinfang 215M May 20 10:17 Dry2014.faa
    -rw-rw-r--  1 jinfang jinfang 510M May 20 10:17 Dry2014.cds
    -rw-rw-r--  1 jinfang jinfang 151M May 20 10:17 Dry2014.gff

Hint: If you apply ``prodigal`` not ``pyrodigal``, use this to fix gene ID in gff files by dbcan_utils (TIMING ~1min)

.. code-block:: shell

    dbcan_utils gff_fix -i Dry2014.faa -g Dry2014.gff
    dbcan_utils gff_fix -i Wet2014.faa -g Wet2014.gff

The value attribution CDS ID in gff file and protein ID in fasta file generated by prodigal is not match, which will lead to failure in CGC prediction.
The command dbcan_utils is designed to fix this issue. The output composes of two files: Dry2014.fix.gff and Wet2014.fix.gff.


Module 2. run_dbcan annotation (Fig. 2) to obtain CAZymes, CGCs, and substrates
```````````````````````````````````````````````````````````````````````````````

**CRITICAL STEP**

Users can skip P5 and P6, and directly run P7 (much slower though), if they want to predict not only CAZymes and CGCs, but also substrates.

P5. CAZyme annotation at the CAZyme family level (TIMING ~10min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    run_dbcan easy_CAZyme --input_raw_data Wet2014.faa --mode protein --output_dir Wet2014.CAZyme --db_dir db
    run_dbcan easy_CAZyme --input_raw_data Dry2014.faa --mode protein --output_dir Dry2014.CAZyme --db_dir db

Two arguments are required for ``run_dbcan``: the input sequence file (faa files) and the sequence type (protein).
By default, ``run_dbcan`` will use three methods (``pyHMMER`` vs ``dbCAN HMMdb``, ``DIAMOND`` vs ``CAZy``, ``pyHMMER`` vs ``dbCAN-sub HMMdb``) for
CAZyme annotation.


The sequence type can be `protein`, `prok`, `meta`. If the input sequence file contains metagenomic contig sequences (`fna` file),
the sequence type has to be `meta`, and `pyrodigal` will be called to predict genes.


P6. CGC prediction (TIMING ~15 min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following commands will re-run run_dbcan to not only predict CAZymes but also CGCs with protein `faa` and gene location `gff` files.

.. code-block:: shell

    run_dbcan easy_CGC --input_raw_data Wet2014.faa --mode protein --output_dir Wet2014.PUL --db_dir db  --input_format NCBI --input_gff Wet2014.gff --gff_type prodigal
    run_dbcan easy_CGC --input_raw_data Dry2014.faa --mode protein --output_dir Dry2014.PUL --db_dir db  --input_format NCBI --input_gff Wet2014.gff --gff_type prodigal


As mentioned above (see Fig. 1), CGC prediction is a featured function added into dbCAN2 in 2018.

.. warning::

    **Creating own gff file**
    If the users would like to create their own ``gff`` file (instead of using Prokka or Prodigal),
    it is important to make sure the value of ID attribute in the ``gff`` file matches the protein ID in the protein ``faa`` file.

    **[Troubleshooting]CGC not found**
    If no result is found in CGC output file, it is most likely because the sequence IDs in ``gff`` file and ``faa`` file do not match.
    Another less likely reason is that the contigs are too short and fragmented and not suitable for CGC prediction.


P7. Substrate prediction for CAZymes and CGCs (TIMING ~5h)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    run_dbcan easy_substrate --input_raw_data Wet2014.faa --mode protein --input_format NCBI --input_gff Wet2014.gff --gff_type prodigal --output_dir Wet2014.dbCAN
    run_dbcan easy_substrate --input_raw_data Dry2014.faa --mode protein --input_format NCBI --input_gff Wet2014.gff --gff_type prodigal --output_dir Dry2014.dbCAN



Box 6. Example output folder content of run_dbcan substrate prediction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    In the output directory (`Output Directory <https://bcb.unl.edu/dbCAN_tutorial/dataset1-Carter2023/individual_assembly/Wet2014.dbCAN/>`_),
    a total of 17 files and 1 folder are generated:




Descriptions of Output Files:
In the output directory, a total of 17 files and 1 folder are generated:

    - ``PUL_blast.out``: BLAST results between CGCs and PULs.
    - ``CGC.faa``: Protein Fasta sequences encoded in all CGCs.
    - ``cgc.gff``: Reformatted from the user input gff file by marking CAZymes, TFs, TCs, and STPs.
    - ``cgc.out``: Raw output of CGC predictions.
    - ``cgc_standard.out``: Simplified version of cgc.out for easy parsing in TSV format. Example columns include:
        1. ``CGC_id``: CGC1
        2. ``type``: CAZyme
        3. ``contig_id``: k141_272079
        4. ``gene_id``: k141_272079_6
        5. ``start``: 5827
        6. ``end``: 7257
        7. ``strand``: -
        8. ``annotation``: GH1

    **Explanation**: The gene Wet2014_00308 encodes a GH1 CAZyme in CGC1 of contig k141_272079. CGC1 also contains other genes, detailed in other rows. Wet2014_00308 is located on the negative strand of k141_272079 from positions 5827 to 7257. The type can be one of four signature gene types (CAZymes, TCs, TFs, STPs) or null type (not annotated as one of the signature genes).

    - ``dbCAN-sub.substrate.tsv``: HMMER search result against dbCAN-sub HMMdb, including a column with CAZyme substrates from fam-substrate-mapping.tsv.
    - ``diamond_results.tsv``: DIAMOND search result against the CAZy annotated protein sequences (CAZyDB.fa).
    - ``dbCAN_hmm_results.tsv``: HMMER search result against dbCAN HMMdb.
    - ``overview.tsv``: Summary of CAZyme annotation from three methods in TSV format. Example columns include:

        1. ``Gene_ID``: k141_355284_4
        2. ``EC#``: 3.2.1.99:3
        3. ``dbCAN``: GH43_4(42-453)
        4. ``dbCAN_sub``: GH43_e149
        5. ``DIAMOND``: GH43_4
        6. ``#ofTools``: 3

    **Explanation**: The protein Wet2014_000076 is annotated by three tools as a CAZyme: GH43_4 (CAZy defined subfamily 4 of GH43) by HMMER vs dbCAN HMMdb, GH43_e149 (eCAMI defined subfamily e149; 'e' indicates it is from eCAMI not CAZy) by HMMER vs dbCAN-sub HMMdb, and GH43_4 by DIAMOND vs CAZy annotated protein sequences. The EC number is extracted from eCAMI, indicating that the eCAMI subfamily GH43_e149 contains 3 member proteins with an EC 3.2.1.99 according to CAZy. The preference order for different assignments is dbCAN > eCAMI/dbCAN-sub > DIAMOND. Refer to dbCAN2 paper11, dbCAN3 paper12, and eCAMI41 for more details.

    - ``stp.out``: HMMER search result against the MiST5 compiled signal transduction protein HMMs from Pfam.
    - ``tf-1.out``: HMMER search result against the DBD6 compiled transcription factor HMMs from Pfam 7.
    - ``tf-2.out``: HMMER search result against the DBD compiled transcription factor HMMs from Superfamily 8.
    - ``tp.out``: DIAMOND search result against the TCDB 9 annotated protein sequences.
    - ``substrate.out``: summary of substrate prediction results for CGCs in TSV format from two approaches3 (dbCAN-PUL blast search and dbCAN-sub majority voting). An example row has the following columns:

        1. ``CGC_ID``: k141_227425|CGC1
        2. ``Best hit PUL_ID in dbCAN-PUL``: PUL0402
        3. ``Substrate of the hit PUL``: xylan
        4. ``Sum of bitscores for homologous gene pairs between CGC and PUL``: 2134.0
        5. ``Types of homologous gene pairs``: TC-TC;CAZyme-CAZyme
        6. ``Substrate predicted by majority voting of CAZymes in CGC``: xylan
        7. ``Voting score``: 2.0

    *Explanation*: The CGC1 of contig k141_227425 has its best hit PUL0402 (from PUL_blast.out) with xylan as substrate (from dbCAN-PUL_12-12-2023.xlsx). Two signature genes are matched between k141_227425|CGC1 and PUL0402: one is a CAZyme and the other is a TC. The sum of blast bit scores of the two homologous pairs (TC-TC and CAZyme-CAZyme) is 2134.0. Hence, the substrate of k141_227425|CGC1 is predicted to be xylan according to dbCAN-PUL blast search. The last two columns are based on the dbCAN-sub result (dbcan-sub.hmm.out), as the file indicates that two CAZymes in k141_227425|CGC1 are predicted to have xylan substrate. The voting score is 2.0, so according to the majority voting rule, k141_227425|CGC1 is predicted to have a xylan substrate.

    *Note*: For many CGCs, only one of the two approaches produces substrate prediction. In some cases, the two approaches produce different substrate assignments. The recommended preference order is dbCAN-PUL blast search > dbCAN-sub majority voting. Refer to dbCAN3 paper12 for more details.

    - ``synteny.pdf``: A folder with syntenic block alignment plots between all CGCs and PULs.
    - ``uniInput``: Renamed Fasta file from input protein sequence file.

Module 3. Read mapping (Fig. 2) to calculate abundance for CAZyme families, subfamilies, CGCs, and substrates
``````````````````````````````````````````````````````````````````````````````````````````````````````````````
P8. Read mapping to all CDS of each sample (TIMING ~20 min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    mkdir samfiles
    bwa index megahit_Wet2014/Wet2014.contigs.fa
    bwa index megahit_Dry2014/Dry2014.contigs.fa
    bwa mem -t 32 -o samfiles/Wet2014.sam megahit_Wet2014/Wet2014.contigs.fa Wet2014_1_val_1.fq.gz Wet2014_2_val_2.fq.gz
    bwa mem -t 32 -o samfiles/Dry2014.sam megahit_Dry2014/Dry2014.contigs.fa Dry2014_1_val_1.fq.gz Dry2014_2_val_2.fq.gz


Reads are mapped to the `contig` files from MEGAHIT.

P9. Sort SAM files by coordinates (TIMING ~4min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    cd samfiles
    samtools sort -@ 32 -o Wet2014.bam Wet2014.sam
    samtools sort -@ 32 -o Dry2014.bam Dry2014.sam
    samtools index Wet2014.bam
    samtools index Dry2014.bam
    rm -rf *sam
    cd ..


Reads are mapped to the `contig` files from MEGAHIT.

P10. Read count calculation for all proteins of each sample using dbcan_utils  (TIMING ~12min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell


    mkdir Wet2014_abund && cd Wet2014_abund
    dbcan_utils cal_coverage -g ../Wet2014.fix.gff -i ../samfiles/Wet2014.bam -o Wet2014.depth.txt -t 6 --overlap_base_ratio 0.2 --mapping_quality 30 --identity 0.98

    cd .. && mkdir Dry2014_abund && cd Dry2014_abund
    dbcan_utils cal_coverage -g ../Dry2014.fix.gff -i ../samfiles/Dry2014.bam -o Dry2014.depth.txt -t 6 --overlap_base_ratio 0.2 --mapping_quality 30 --identity 0.98
    cd ..


How to filter the mapped reads are the key step of the abundance calculation. Several parameters can be set to filter the reads. The command dbcan_utils provides --overlap_base_ratio, --mapping_quality and --identity to filter the mapped reads which will may produce more reliable abundance. The default values for these three parameters are 0.2, 30 and 0.98, respectively. Read counts are saved in depth.txt files of each sample.



P11. Read count calculation for a given region of contigs using Samtools (TIMING ~2min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    cd Wet2014_abund
    samtools index ../samfiles/Wet2014.bam
    samtools depth -r k141_41392:152403-165349 ../samfiles/Wet2014.bam > Wet2014.cgc.depth.txt




The parameter -r k141_41392:152403-165349 specifies a region in a contig. For any CGC, its positional range can be found in the file cgc_standard.out produced by run_dbcan (Box 6).

.. warning::
    The contig IDs are automatically generated by MEGAHIT. There is a small chance that a same contig ID appears in both samples. However, the two contigs in the two samples do not match each other even the ID is the same. For example, the contig ID k141_4139 is most likely only found in the Wet2014 sample. Even if there is a k141_41392 in Dry2014, the actual contigs in two samples are different.
    The depth.txt files contain the raw read counts for the specified region. It may show a warning but will not effect the result:
    [W::hts_idx_load3] The index file is older than the data file: ../samfiles/Wet2014.bam.bai


P12. dbcan_utils to calculate the abundance of CAZyme families, subfamilies, CGCs, and substrates (TIMING ~1min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    dbcan_utils fam_abund -bt Wet2014.depth.txt -i ../Wet2014.dbCAN -a TPM
    dbcan_utils fam_substrate_abund -bt Wet2014.depth.txt -i ../Wet2014.dbCAN -a TPM
    dbcan_utils CGC_abund -bt Wet2014.depth.txt -i ../Wet2014.dbCAN -a TPM
    dbcan_utils CGC_substrate_abund -bt Wet2014.depth.txt -i ../Wet2014.dbCAN -a TPM

    cd .. && cd Dry2014_abund
    dbcan_utils fam_abund -bt Dry2014.depth.txt -i ../Dry2014.dbCAN -a TPM
    dbcan_utils fam_substrate_abund -bt Dry2014.depth.txt -i ../Dry2014.dbCAN -a TPM
    dbcan_utils CGC_abund -bt Dry2014.depth.txt -i ../Dry2014.dbCAN -a TPM
    dbcan_utils CGC_substrate_abund -bt Dry2014.depth.txt -i ../Dry2014.dbCAN -a TPM
    cd ..


We developed a set of Python scripts as dbcan_utils (included in the run_dbcan package) to take the raw read counts for all genes as input and output the normalized abundances (Box 7) of CAZyme families, subfamilies, CGCs, and substrates (Fig. 4). The parameter -a TPM can also be two other metrics: RPM, or RPKM61.

- **RPKM** = # of mapped reads to a gene G / [(total # of mapped reads to all genes /106) x (gene G length/1000)]

- **RPM** = # of mapped reads to a gene G / (total # of mapped reads to all genes/106).

- **TPM** = [# of mapped reads to a gene G / (gene G length/1000)] / sum [# of mapped reads to each gene / (the gene length/1000)].


Box 7. Example output of dbcan_utils
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
As an example, the Wet2014_abund folder (https://bcb.unl.edu/dbCAN_tutorial/dataset1-Carter2023/individual_assembly/Wet2014_abund/) has 7 TSV files:

.. code-block:: shell

    -rw-rw-r--  1 jinfang jinfang 226K Jun 16 09:28 CGC_abund.out
    -rw-rw-r--  1 jinfang jinfang 2.3K Jun 16 09:28 CGC_substrate_majority_voting.out
    -rw-rw-r--  1 jinfang jinfang  20K Jun 16 09:28 CGC_substrate_PUL_homology.out
    -rw-rw-r--  1 jinfang jinfang 2.8K Jun 16 09:28 EC_abund.out
    -rw-rw-r--  1 jinfang jinfang 4.1K Jun 16 09:28 fam_abund.out
    -rw-rw-r--  1 jinfang jinfang  60K Jun 16 09:28 fam_substrate_abund.out
    -rw-rw-r--  1 jinfang jinfang  30K Jun 16 09:28 subfam_abund.out


Explanation of columns in these TSV files is as follows:

    - ``fam_abund.out``: CAZy family (from HMMER vs dbCAN HMMdb), sum of TPM, number of CAZymes in the family.
    - ``subfam_abund.out``: eCAMI subfamily (from HMMER vs dbCAN-sub HMMdb), sum of TPM, number of CAZymes in the subfamily.
    - ``EC_abund.out``: EC number (extracted from dbCAN-sub subfamily), sum of TPM, number of CAZymes with the EC.
    - ``fam_substrate_abund.out``: Substrate (from HMMER vs dbCAN-sub HMMdb), sum of TPM (all CAZymes in this substrate group), GeneID (all CAZyme IDs in this substrate group).
    - ``CGC_abund.out``: CGC_ID (e.g., k141_338400|CGC1), mean of TPM (all genes in the CGC), Seq_IDs (IDs of all genes in the CGC), TPM (of all genes in the CGC), Families (CAZyme family or other signature gene type of all genes in the CGC).
    - ``CGC_substrate_PUL_homology.out``: Substrate (from dbCAN-PUL blast search), sum of TPM, CGC_IDs (all CGCs predicted to have the substrate from dbCAN-PUL blast search), TPM (of CGCs in this substrate group).
    - ``CGC_substrate_majority_voting.out``: Substrate (from dbCAN-sub majority voting), sum of TPM, CGC_IDs (all CGCs predicted to have the substrate from dbCAN-sub majority voting), TPM (of CGCs in this substrate group).

.. warning::
    Proteins from multiple samples can be combined to generate a non-redundant set of proteins (Box 8). This may reduce the runtime for the run_dbcan step (step4), as only one faa file will be processed. However, this does not work for the CGC prediction, as contigs (fna files) from each sample will be needed. Therefore, this step is recommended if users only want the CAZyme annotation, and not recommended if CGCs are also to be predicted.


.. warning::
    As shown in Fig. 1 (step3), proteins from multiple samples can be combined to generate a non-redundant set of proteins (Box 8). This may reduce the runtime for the run_dbcan step (step4), as only one faa file will be processed. However, this does not work for the CGC prediction, as contigs (fna files) from each sample will be needed. Therefore, this step is recommended if users only want the CAZyme annotation, and not recommended if CGCs are also to be predicted.

Box 8. Combine proteins from multiple samples (optional TIMING ~3h)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This protein sequence clustering step will create a mapping table with sequence cluster ID and protein IDs from each sample.

.. code-block:: shell

    mkdir mmseqs_cluster && cd mmseqs_cluster
    ln -s ../db .
    cat ../Wet2014.faa ../Dry2014.faa > Dry_Wet.faa
    mmseqs easy-cluster --threads 32 -c 0.95 --min-seq-id 0.95 --cov-mode 2 Dry_Wet.faa Dry_Wet_cluster tmp
    mv Dry_Wet_cluster_cluster_rep.fasta Dry_Wet.cluster.faa

This Dry_Wet.cluster.faa file now contains the Fasta sequences of representative proteins of all mmseqs2 clusters, i.e., the non-redundant protein sequences from the two samples. The mapping table file Dry_Wet_cluster_cluster.tsv contains two columns, mmseqs2 cluster representative protein ID, and protein IDs in the cluster. Users can query this table to find the correspondence between CAZymes in Wet2014 and in Dry2014.




Module 4: dbcan_plot for data visualization (Fig. 2) of abundances of CAZymes, CGCs, and substrates (TIMING variable)
`````````````````````````````````````````````````````````````````````````````````````````````````````````````````````
**CRITICAL STEP**

To visualize the CAZyme annotation result, we provide a set of Python scripts as ``dbcan_plot`` to make publication-quality plots with the ``dbcan_utils`` results as the input. The ``dbcan_plot`` scripts are included in the ``run_dbcan`` package. Once the plots are made in PDF format, they can be transferred to users' Windows or Mac computers for visualization.

Five data folders will be needed as the input for ``dbcan_plot``:

1. Two abundance folders: ``Wet2014_abund`` and ``Dry2014_abund``.
2. Two CAZyme annotation folders: ``Wet2014.dbCAN`` and ``Dry2014.dbCAN``.
3. The ``dbCAN-PUL`` folder (located under the db folder, released from ``dbCAN-PUL.tar.gz``).



P13. Heatmap for CAZyme substrate abundance across samples (Fig. 6A) (TIMING 1min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    dbcan_plot heatmap_plot --samples Wet2014,Dry2014 -i Wet2014_abund/fam_substrate_abund.out,Dry2014_abund/fam_substrate_abund.out --show_abund --top 20

Here we plot the top 20 substrates in the two samples (Fig. 6A). The input files are the two CAZyme substrate abundance files calculated based on
dbCAN-sub result. The default heatmap is ranked by substrate abundances. To rank the heatmap according to abundance profile using
the clustermap function of the seaborn package (https://github.com/mwaskom/seaborn), users can invoke the ``--cluster_map`` parameter.

P14. Barplot for CAZyme family/subfamily/EC abundance across samples (Fig. 6B,C) (TIMING 1min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    dbcan_plot bar_plot --samples Wet2014,Dry2014 --vertical_bar --top 20 -i Wet2014_abund/fam_abund.out,Dry2014_abund/fam_abund.out --pdf fam_abund.pdf --db_dir db
    dbcan_plot bar_plot --samples Wet2014,Dry2014 --vertical_bar --top 20 -i Wet2014_abund/subfam_abund.out,Dry2014_abund/subfam_abund.out --pdf subfam_abund.pdf --db_dir db
    dbcan_plot bar_plot --samples Wet2014,Dry2014 --vertical_bar --top 20 -i Wet2014_abund/EC_abund.out,Dry2014_abund/EC_abund.out --pdf EC_abund.pdf --db_dir db



Users can choose to generate a barplot instead of heatmap using the bar_plot method.

P15. Synteny plot between a CGC and its best PUL hit with read mapping coverage to CGC (Fig. 6D) (TIMING 1min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    dbcan_plot CGC_synteny_coverage_plot -i Wet2014.dbCAN --cgcid 'k141_41392|CGC3' --readscount Wet2014_abund/Wet2014.cgc.depth.txt --db_dir db

The Wet2014.dbCAN folder contains the PUL.out file. Using this file, the cgc_standard.out file, and the best PUL’s gff file in dbCAN-PUL.tar.gz, the CGC_synteny_plot method will create the CGC-PUL synteny plot. The –cgcid parameter is required to specify which CGC to plot (k141_41392|CGC3 in this example). The Wet2014.cgc.depth.txt file is used to plot the read mapping coverage.
If users only want to plot the CGC structure:

.. code-block:: shell

    dbcan_plot CGC_plot -i Wet2014.dbCAN --cgcid 'k141_41392|CGC3'

If users only want to plot the CGC structure plus the read mapping coverage:

.. code-block:: shell

    dbcan_plot CGC_coverage_plot -i Wet2014.dbCAN --cgcid 'k141_41392|CGC3' --readscount Wet2014_abund/Wet2014.cgc.depth.txt

If users only want to plot the synteny between the CGC and PUL:

.. code-block:: shell

    dbcan_plot CGC_synteny_plot -i Wet2014.dbCAN --cgcid 'k141_41392|CGC3' --db_dir db

.. warning::

    The CGC IDs in different samples do not match each other. For example, specifying -i Wet2014.dbCAN is to plot the 'k141_41392|CGC3' in the Wet2014 sample. The 'k141_41392|CGC3' in the Dry2014 sample most likely does not exist, and even it does, the CGC has a different sequence even if the ID is the same.


Troubleshooting
---------------

We provide Table 3 to list possible issues and solutions. Users can also post issues on run_dbcan GitHub site.

Procedure Timing
----------------

The estimated time for completing each step of the protocol on the Carter2023 dataset is as follows:

- **Step P1. Contamination checking ~10min.**
- **Step P2. Raw reads processing ~20min.**
- **Step P3. Metagenome assembly ~4h20min.**
- **Step P4. Gene models prediction ~15min.**
- **Step P5. CAZyme annotation ~10min.**
- **Step P6. PUL prediction ~15min.**
- **Step P7. Substrate prediction both for CAZyme and PUL ~5h.**
- **Step P8-P11. Reads mapping ~36min.**
- **Step P12. Abundance estimation ~1min.**
- **Step P13-P15. Data visualization ~3min.**

Running this protocol on the Carter2023 dataset will take ~12h on a Linux computer with 40 CPUs and 128GB of RAM. The most time-consuming step is P7 (substrate prediction for CGCs and CAZymes). .  If users choose not to predict substates, this step will take ~15min. RAM usage was not specifically monitored during the execution. The step with the highest RAM usage is likely P3 (read assembly).
