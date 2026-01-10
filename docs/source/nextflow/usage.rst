.. _nextflow-usage:

Nextflow Pipeline: Usage
==========================

This document provides comprehensive guidance on using the Nextflow pipeline for CAZyme annotation in microbiome data. The pipeline supports multiple analysis modes and can be configured to suit various computational environments and data types.

Introduction
------------

This pipeline supports three analysis modes for CAZyme annotation in microbiome data:

- **Short reads** (``--type shortreads``): Assembly-based analysis for Illumina short-read data using MEGAHIT
- **Long reads** (``--type longreads``): Assembly-based analysis for PacBio/Nanopore long-read data using Flye
- **Assembly free** (``--type assemfree``): Direct annotation without assembly using seqtk and DIAMOND blastx

The assembly free mode is particularly useful for large datasets where assembly is computationally expensive or when you want to avoid potential assembly artifacts. It directly converts reads to FASTA format and uses DIAMOND blastx to search against the CAZyme database, followed by abundance calculation using dbcan_asmfree.

Samplesheet Input
-----------------

You will need to create a samplesheet file containing information about the samples you would like to analyze before running the pipeline. Use the ``--input`` parameter to specify the location of this file. The samplesheet must be a comma-separated file with 3 columns and a header row, as shown in the examples below.

.. code-block:: bash

   --input '[path to samplesheet file]'

Multiple Runs of the Same Sample
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``sample`` identifiers must be identical when you have re-sequenced the same sample more than once (e.g., to increase sequencing depth). The pipeline will automatically concatenate the raw reads from all runs before performing any downstream analysis. Below is an example for the same sample sequenced across 3 lanes:

.. code-block:: text

   sample,fastq_1,fastq_2
   CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
   CONTROL_REP1,AEG588A1_S1_L003_R1_001.fastq.gz,AEG588A1_S1_L003_R2_001.fastq.gz
   CONTROL_REP1,AEG588A1_S1_L004_R1_001.fastq.gz,AEG588A1_S1_L004_R2_001.fastq.gz

Full Samplesheet Format
~~~~~~~~~~~~~~~~~~~~~~~~

The pipeline will automatically detect whether a sample is single-end or paired-end using the information provided in the samplesheet. The samplesheet can have as many additional columns as desired; however, there is a strict requirement for the first 3 columns to match those defined in the table below.

A complete samplesheet file consisting of both single-end and paired-end data may look like the example below. This example includes 6 samples, where ``TREATMENT_REP3`` has been sequenced twice.

.. code-block:: text

   sample,fastq_1,fastq_2
   CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
   CONTROL_REP2,AEG588A2_S2_L002_R1_001.fastq.gz,AEG588A2_S2_L002_R2_001.fastq.gz
   CONTROL_REP3,AEG588A3_S3_L002_R1_001.fastq.gz,AEG588A3_S3_L002_R2_001.fastq.gz
   TREATMENT_REP1,AEG588A4_S4_L003_R1_001.fastq.gz,
   TREATMENT_REP2,AEG588A5_S5_L003_R1_001.fastq.gz,
   TREATMENT_REP3,AEG588A6_S6_L003_R1_001.fastq.gz,
   TREATMENT_REP3,AEG588A6_S6_L004_R1_001.fastq.gz,

.. list-table:: Samplesheet Column Descriptions
   :widths: 20 80
   :header-rows: 1

   * - Column
     - Description
   * - ``sample``
     - Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (``_``).
   * - ``fastq_1``
     - Full path to FastQ file for Illumina short reads 1. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".
   * - ``fastq_2``
     - Full path to FastQ file for Illumina short reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".

An example samplesheet file should be created following the format described above. Ensure that all file paths are absolute or relative to your working directory.

Running the Pipeline
---------------------

The typical command for running the pipeline is as follows:

.. code-block:: bash

   nextflow run nf-core/dbcanmicrobiome --input ./samplesheet.csv --outdir ./results -profile docker

This command launches the pipeline with the ``docker`` configuration profile, which ensures reproducible execution using containerized software. See below for more information about available profiles and configuration options.

Note that the pipeline will create the following files in your working directory:

.. code-block:: bash

   work                # Directory containing the nextflow working files
   <OUTDIR>            # Finished results in specified location (defined with --outdir)
   .nextflow_log       # Log file from Nextflow
   # Other nextflow hidden files, eg. history of pipeline runs and old logs.

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a ``yaml`` or ``json`` file via ``-params-file <file>``.

.. warning::
   Do not use ``-c <file>`` to specify parameters as this will result in errors. Custom config files specified with ``-c`` must only be used for `tuning process resource specifications <https://nf-co.re/docs/usage/configuration#tuning-workflow-resources>`_, other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

.. code-block:: bash

   nextflow run nf-core/dbcanmicrobiome -profile docker -params-file params.yaml

with:

.. code-block:: yaml
   :name: params.yaml

   input: './samplesheet.csv'
   outdir: './results/'
   <...>

You can also generate such ``YAML``/``JSON`` parameter files using the `nf-core launch tool <https://nf-co.re/launch>`_ or by manually creating them based on your requirements.

Updating the Pipeline
~~~~~~~~~~~~~~~~~~~~~~

When you run the pipeline command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline subsequently, it will always use the cached version if available, even if the pipeline has been updated since the last run. To ensure that you're running the latest version of the pipeline, regularly update the cached version using:

.. code-block:: bash

   nextflow pull nf-core/dbcanmicrobiome

Reproducibility
~~~~~~~~~~~~~~~

It is highly recommended to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used, providing reproducibility across different runs. If you consistently use the same version tag, you will always run the same version of the pipeline, even if there have been updates to the code since your initial run.

To ensure reproducibility, specify the pipeline version when running. First, check the `nf-core/dbcanmicrobiome releases page <https://github.com/nf-core/dbcanmicrobiome/releases>`_ to find the latest stable version (e.g., ``1.3.1``). Then specify this version when running the pipeline using the ``-r`` flag (single hyphen), for example: ``-r 1.3.1``. You can switch to any other version by changing the version number after the ``-r`` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and reuse `parameter files <#running-the-pipeline>`_ to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

.. tip::
   If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

Core Nextflow Arguments
------------------------

.. note::
   These options are part of Nextflow and use a *single* hyphen (pipeline parameters use a double-hyphen)

``-profile``
~~~~~~~~~~~~

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

.. important::
   We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline dynamically loads configurations from the `nf-core/configs repository <https://github.com/nf-core/configs>`_ at runtime, providing pre-configured profiles for various institutional clusters and compute environments. For more information and to check if your system is supported, please see the `nf-core/configs documentation <https://github.com/nf-core/configs#documentation>`_.

Note that multiple profiles can be loaded, for example: ``-profile test,docker`` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If ``-profile`` is not specified, the pipeline will run locally and expect all required software to be installed and available on the system ``PATH``. This approach is *not* recommended, as it can lead to inconsistent results across different machines due to variations in software versions and system configurations.

- ``test``
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- ``docker``
  - A generic configuration profile to be used with `Docker <https://docker.com/>`_
- ``singularity``
  - A generic configuration profile to be used with `Singularity <https://sylabs.io/docs/>`_
- ``podman``
  - A generic configuration profile to be used with `Podman <https://podman.io/>`_
- ``shifter``
  - A generic configuration profile to be used with `Shifter <https://nersc.gitlab.io/development/shifter/how-to-use/>`_
- ``charliecloud``
  - A generic configuration profile to be used with `Charliecloud <https://hpc.github.io/charliecloud/>`_
- ``apptainer``
  - A generic configuration profile to be used with `Apptainer <https://apptainer.org/>`_
- ``wave``
  - A generic configuration profile to enable `Wave <https://seqera.io/wave/>`_ containers. Use together with one of the above (requires Nextflow ``24.03.0-edge`` or later).
- ``conda``
  - A generic configuration profile to be used with `Conda <https://conda.io/docs/>`_. Please only use Conda as a last resort when containerization solutions (Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer) are not available or feasible.

Short Reads Analysis Modes: Subsampling and Co-assembly
--------------------------------------------------------

The following options are only available when using ``--type shortreads`` mode.

Subsampling Mode
~~~~~~~~~~~~~~~~

- **Purpose**: Downsample each sample before assembly to reduce computational requirements and enable quick pipeline validation.
- **Parameters**:
  - ``--subsample``: Enable subsampling mode.
  - ``--subsample_size``: Number of reads per file to retain (default: ``20000000`` as specified in the pipeline configuration).
- **Behavior**:
  - Subsampling is applied per-sample before MEGAHIT assembly using ``seqtk sample``.
  - This mode is mutually exclusive with ``--coassembly``.
- **Example**:

.. code-block:: bash

   nextflow run nf-core/dbcanmicrobiome \
     --type shortreads \
     --input samplesheet.csv \
     --outdir results_subsample \
     --subsample \
     --subsample_size 5000000 \
     -profile docker

Co-assembly Mode
~~~~~~~~~~~~~~~~

- **Purpose**: Co-assemble all short-read samples together to improve contig continuity and enhance detection of shared genomic features across samples.
- **Parameters**:
  - ``--coassembly``: Enable co-assembly mode across all samples.
- **Requirements & Behavior**:
  - Requires at least 2 samples; the pipeline will produce an error if fewer samples are provided.
  - Combines all reads (preserving paired-end or single-end structure) and performs a single MEGAHIT assembly.
  - This mode is mutually exclusive with ``--subsample``.
- **Example**:

.. code-block:: bash

   nextflow run nf-core/dbcanmicrobiome \
     --type shortreads \
     --input samplesheet.csv \
     --outdir results_coassembly \
     --coassembly \
     -profile docker

``-resume``
~~~~~~~~~~~

Specify this flag when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are identical, allowing the pipeline to continue from where it previously stopped. For inputs to be considered identical, both the file names and file contents must match exactly. For more information about this parameter and how Nextflow handles caching, see `this blog post <https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html>`_.

You can also supply a run name to resume a specific run: ``-resume [run-name]``. Use the ``nextflow log`` command to show previous run names.

``-c``
~~~~~~

Specify the path to a specific config file (this is a core Nextflow command). See the `nf-core website documentation <https://nf-co.re/usage/configuration>`_ for more information.

Custom Configuration
---------------------

Resource Requests
~~~~~~~~~~~~~~~~~~

While the default resource requirements set within the pipeline are designed to work for most users and input datasets, you may need to customize the compute resources requested by the pipeline. Each step in the pipeline has default requirements for the number of CPUs, memory, and time allocation. For most pipeline steps, if a job exits with certain error codes (typically related to resource exhaustion), it will automatically be resubmitted with increased resource requests (2x the original, then 3x the original). If the job still fails after the third attempt, the pipeline execution will be stopped.

To change the resource requests, please see the `max resources <https://nf-co.re/docs/usage/configuration#max-resources>`_ and `tuning workflow resources <https://nf-co.re/docs/usage/configuration#tuning-workflow-resources>`_ section of the nf-core website.

Custom Containers
~~~~~~~~~~~~~~~~~~

In some cases, you may wish to change the container or conda environment used by a particular pipeline step for a specific tool. By default, nf-core pipelines use containers and software from the `biocontainers <https://biocontainers.pro/>`_ or `bioconda <https://bioconda.github.io/>`_ projects. However, in some cases, the pipeline-specified version may be outdated or you may need a different version for compatibility reasons.

To use a different container from the default container or conda environment specified in a pipeline, please see the `updating tool versions <https://nf-co.re/docs/usage/configuration#updating-tool-versions>`_ section of the nf-core website.

Custom Tool Arguments
~~~~~~~~~~~~~~~~~~~~~~

A pipeline may not always support every possible argument or option of a particular tool used within it. Fortunately, nf-core pipelines provide flexibility for users to insert additional parameters that are not included by default in the pipeline configuration.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the `customising tool arguments <https://nf-co.re/docs/usage/configuration#customising-tool-arguments>`_ section of the nf-core website.

nf-core/configs
~~~~~~~~~~~~~~~

In most cases, you will only need to create a custom config file as a one-time setup. However, if you and others within your organization are likely to run nf-core pipelines regularly and need to use the same settings consistently, it may be beneficial to request that your custom config file be uploaded to the ``nf-core/configs`` git repository. Before submitting, please test that the config file works correctly with your pipeline of choice using the ``-c`` parameter. You can then create a pull request to the ``nf-core/configs`` repository including your config file, associated documentation (see examples in `nf-core/configs/docs <https://github.com/nf-core/configs/tree/master/docs>`_), and an amendment to `nfcore_custom.config <https://github.com/nf-core/configs/blob/master/nfcore_custom.config>`_ to include your custom profile.

See the main `Nextflow documentation <https://www.nextflow.io/docs/latest/config.html>`_ for more information about creating your own configuration files.

If you have any questions or issues regarding custom configurations, please send a message on `nf-core Slack <https://nf-co.re/join/slack>`_ in the `#configs channel <https://nfcore.slack.com/channels/configs>`_.

Running in the Background
--------------------------

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must remain active until the pipeline execution is complete.

The Nextflow ``-bg`` flag launches Nextflow in the background, detached from your terminal, so that the workflow continues running even if you log out of your session. When using this flag, logs are automatically saved to a file.

Alternatively, you can use terminal multiplexers such as ``screen`` or ``tmux`` to create a detached session that you can reconnect to at a later time. Some HPC setups also allow you to run Nextflow within a cluster job submitted through your job scheduler, from which it can submit additional jobs.

Nextflow Memory Requirements
-----------------------------

In some cases, the Nextflow Java virtual machine may request a large amount of memory, which can cause issues on systems with limited resources.

We recommend adding the following environment variable to limit Nextflow's memory usage. This should be added to your shell configuration file (typically ``~/.bashrc`` or ``~/.bash_profile``):

.. code-block:: bash

   NXF_OPTS='-Xms1g -Xmx4g'

