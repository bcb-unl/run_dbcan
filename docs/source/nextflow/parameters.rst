.. _nextflow-parameters:

Nextflow Pipeline: Parameters Reference
========================================

This document provides a comprehensive reference for all parameters available in the Nextflow pipeline for CAZyme annotation in microbiome data. Parameters are organized by functional category.

.. note::
   Before running the pipeline, make sure you have cloned the repository. See :ref:`nextflow-usage` for installation instructions. All examples assume you are in the ``dbcan-nf`` directory or use the full path to ``main.nf``.

Input/Output Options
--------------------

.. list-table:: Input/Output Parameters
   :widths: 20 10 10 60
   :header-rows: 1

   * - Parameter
     - Type
     - Required
     - Description
   * - ``--input``
     - string
     - Yes
     - Path to comma-separated file containing information about the samples in the experiment. Must have 3 columns (sample, fastq_1, fastq_2) and a header row. See :ref:`nextflow-usage` for samplesheet format details.
   * - ``--outdir``
     - string
     - Yes
     - The output directory where the results will be saved. Use absolute paths for Cloud infrastructure.
   * - ``--email``
     - string
     - No
     - Email address for completion summary. Set this to receive a summary email when the workflow exits. Can be set in ``~/.nextflow/config`` to avoid specifying on every run.
   * - ``--email_on_fail``
     - string
     - No
     - Email address for completion summary, only sent when pipeline fails.
   * - ``--plaintext_email``
     - boolean
     - No
     - Send plain-text email instead of HTML. Default: ``false``.

Analysis Mode Selection
------------------------

.. list-table:: Mode Selection Parameters
   :widths: 20 10 10 60
   :header-rows: 1

   * - Parameter
     - Type
     - Default
     - Description
   * - ``--type``
     - string
     - ``shortreads``
     - Analysis mode to use. Options: ``shortreads``, ``longreads``, ``assemfree``. See mode-specific documentation for details.

Quality Control Options
-------------------------

.. list-table:: Quality Control Parameters
   :widths: 20 10 10 60
   :header-rows: 1

   * - Parameter
     - Type
     - Default
     - Description
   * - ``--skip_fastqc``
     - boolean
     - ``false``
     - Skip FastQC quality control analysis. When enabled, FastQC steps are bypassed.
   * - ``--skip_trimming``
     - boolean
     - ``false``
     - Skip TrimGalore adapter trimming. When enabled, trimming steps are bypassed.

Kraken2 Options
---------------

.. list-table:: Kraken2 Parameters
   :widths: 20 10 10 60
   :header-rows: 1

   * - Parameter
     - Type
     - Default
     - Description
   * - ``--skip_kraken_extraction``
     - boolean
     - ``false``
     - Skip Kraken2 taxonomic classification and read extraction. When enabled, all reads are used without filtering.
   * - ``--kraken_db``
     - string
     - ``null``
     - Path to custom Kraken2 database directory. If not specified, the pipeline will build a standard database automatically.
   * - ``--kraken_tax``
     - string
     - ``9606``
     - NCBI taxonomy ID for taxonomic filtering. Default ``9606`` corresponds to *Homo sapiens*. Reads matching this taxon are extracted for downstream analysis.

Assembly Options (Short Reads Mode)
------------------------------------

These parameters are only applicable when ``--type shortreads`` is used.

.. list-table:: Short Reads Assembly Parameters
   :widths: 20 10 10 60
   :header-rows: 1

   * - Parameter
     - Type
     - Default
     - Description
   * - ``--subsample``
     - boolean
     - ``false``
     - Enable subsampling mode. Downsample each sample before assembly to reduce computational requirements. Mutually exclusive with ``--coassembly``. See :ref:`shortreads-subsample` for details.
   * - ``--subsample_size``
     - integer
     - ``20000000``
     - Number of reads per file to retain when subsampling is enabled. Only used when ``--subsample`` is set.
   * - ``--coassembly``
     - boolean
     - ``false``
     - Enable co-assembly mode. Combine all samples and perform a single MEGAHIT assembly. Requires at least 2 samples. Mutually exclusive with ``--subsample``. See :ref:`shortreads-coassembly` for details.

Long Reads Options
------------------

These parameters are only applicable when ``--type longreads`` is used.

.. list-table:: Long Reads Parameters
   :widths: 20 10 10 60
   :header-rows: 1

   * - Parameter
     - Type
     - Default
     - Description
   * - ``--flye_mode``
     - string
     - ``--pacbio-hifi``
     - Flye assembly mode. Options include ``--pacbio-hifi``, ``--pacbio-raw``, ``--nano-raw``, ``--nano-hq``. See Flye documentation for details.

Database Options
----------------

.. list-table:: Database Parameters
   :widths: 20 10 10 60
   :header-rows: 1

   * - Parameter
     - Type
     - Default
     - Description
   * - ``--dbcan_db``
     - string
     - ``null``
     - Path to custom dbCAN database directory. If not specified, the pipeline will download the database automatically.

MultiQC Options
---------------

.. list-table:: MultiQC Parameters
   :widths: 20 10 10 60
   :header-rows: 1

   * - Parameter
     - Type
     - Default
     - Description
   * - ``--multiqc_title``
     - string
     - ``null``
     - MultiQC report title. Printed as page header, used for filename if not otherwise specified.
   * - ``--multiqc_config``
     - string
     - ``null``
     - Custom config file to supply to MultiQC.
   * - ``--multiqc_logo``
     - string
     - ``null``
     - Custom logo file to supply to MultiQC. File name must also be set in the MultiQC config file.
   * - ``--multiqc_methods_description``
     - string
     - ``null``
     - Custom MultiQC YAML file containing HTML including a methods description.
   * - ``--max_multiqc_email_size``
     - string
     - ``25.MB``
     - File size limit when attaching MultiQC reports to summary emails.

Generic Options
---------------

These options are common to all nf-core pipelines and are typically set in a Nextflow config file.

.. list-table:: Generic Parameters
   :widths: 20 10 10 60
   :header-rows: 1

   * - Parameter
     - Type
     - Default
     - Description
   * - ``--publish_dir_mode``
     - string
     - ``copy``
     - Method used to save pipeline results to output directory. Options: ``symlink``, ``rellink``, ``link``, ``copy``, ``copyNoFollow``, ``move``.
   * - ``--monochrome_logs``
     - boolean
     - ``false``
     - Do not use coloured log outputs.
   * - ``--hook_url``
     - string
     - ``null``
     - Incoming hook URL for messaging service. Currently, MS Teams and Slack are supported.
   * - ``--validate_params``
     - boolean
     - ``true``
     - Boolean whether to validate parameters against the schema at runtime.

Parameter Usage Examples
-------------------------

Basic short reads analysis:

.. code-block:: bash

   nextflow run main.nf \
     --input samplesheet.csv \
     --outdir results \
     --type shortreads \
     -profile docker

Short reads with subsampling:

.. code-block:: bash

   nextflow run main.nf \
     --input samplesheet.csv \
     --outdir results \
     --type shortreads \
     --subsample \
     --subsample_size 5000000 \
     -profile docker

Long reads analysis:

.. code-block:: bash

   nextflow run main.nf \
     --input samplesheet.csv \
     --outdir results \
     --type longreads \
     --flye_mode --nano-raw \
     -profile docker

Assembly-free analysis:

.. code-block:: bash

   nextflow run main.nf \
     --input samplesheet.csv \
     --outdir results \
     --type assemfree \
     -profile docker

With custom databases:

.. code-block:: bash

   nextflow run main.nf \
     --input samplesheet.csv \
     --outdir results \
     --type shortreads \
     --dbcan_db /path/to/dbcan_db \
     --kraken_db /path/to/kraken_db \
     -profile docker

Skipping quality control steps:

.. code-block:: bash

   nextflow run main.nf \
     --input samplesheet.csv \
     --outdir results \
     --type shortreads \
     --skip_fastqc \
     --skip_trimming \
     -profile docker

Using parameter files:

.. code-block:: bash

   nextflow run main.nf \
     -profile docker \
     -params-file params.yaml

With ``params.yaml``:

.. code-block:: yaml

   input: './samplesheet.csv'
   outdir: './results/'
   type: 'shortreads'
   subsample: true
   subsample_size: 5000000
   skip_kraken_extraction: false
   kraken_tax: '9606'
