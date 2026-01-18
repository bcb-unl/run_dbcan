.. _shortreads-subsample:

Short Reads: Subsampling Mode
==============================

The subsampling mode allows you to downsample sequencing reads before assembly, reducing computational requirements and enabling quick pipeline validation on large datasets.

Overview
--------

Subsampling mode is useful when:

- You have very large datasets and want to reduce computational time and memory usage
- You want to quickly validate the pipeline workflow before running on full datasets
- You need to test different parameters on a subset of data
- Computational resources are limited

When subsampling is enabled, the pipeline uses ``seqtk sample`` to randomly downsample reads from each input file before assembly. This happens per-sample, so each sample is independently subsampled.

Activation
----------

To enable subsampling mode, use the ``--subsample`` flag along with ``--type shortreads``:

.. code-block:: bash

   nextflow run nf-core/dbcanmicrobiome \
     --input samplesheet.csv \
     --outdir results \
     --type shortreads \
     --subsample \
     -profile docker

Parameters
----------

.. list-table:: Subsampling Parameters
   :widths: 20 10 10 60
   :header-rows: 1

   * - Parameter
     - Type
     - Default
     - Description
   * - ``--subsample``
     - boolean
     - ``false``
     - Enable subsampling mode. Must be used with ``--type shortreads``.
   * - ``--subsample_size``
     - integer
     - ``20000000``
     - Number of reads per file to retain. Applied to each FASTQ file independently.

Usage Examples
--------------

Basic Subsampling
~~~~~~~~~~~~~~~~~

Subsample to 20 million reads per file (default):

.. code-block:: bash

   nextflow run nf-core/dbcanmicrobiome \
     --input samplesheet.csv \
     --outdir results \
     --type shortreads \
     --subsample \
     -profile docker

Custom Subsampling Size
~~~~~~~~~~~~~~~~~~~~~~~~

Subsample to 5 million reads per file:

.. code-block:: bash

   nextflow run nf-core/dbcanmicrobiome \
     --input samplesheet.csv \
     --outdir results \
     --type shortreads \
     --subsample \
     --subsample_size 5000000 \
     -profile docker

Subsampling with Other Options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Combine subsampling with other parameters:

.. code-block:: bash

   nextflow run nf-core/dbcanmicrobiome \
     --input samplesheet.csv \
     --outdir results \
     --type shortreads \
     --subsample \
     --subsample_size 10000000 \
     --skip_kraken_extraction \
     -profile docker

Behavior
--------

How Subsampling Works
~~~~~~~~~~~~~~~~~~~~~

1. **Per-File Subsampling**: Each FASTQ file (R1 and R2 for paired-end) is independently subsampled
2. **Random Sampling**: Uses ``seqtk sample`` with random seed for reproducible results
3. **Preserves Pairing**: For paired-end data, both R1 and R2 files are subsampled to maintain pairing
4. **Sample Naming**: Subsampled samples are renamed with ``_subsample`` suffix (e.g., ``sample1_dna_subsample``)

Limitations
~~~~~~~~~~~

- **Mutually Exclusive**: Cannot be used together with ``--coassembly``
- **No RNA-seq**: RNA-seq processing is automatically disabled when subsampling is enabled
- **Reduced Coverage**: Subsampling reduces sequencing depth, which may affect assembly quality and gene detection

Output Files
------------

Output files follow the same structure as standard short reads mode, but with modified sample names:

- Sample IDs are appended with ``_subsample`` suffix
- All downstream files (assembly, annotation, abundance) use the subsampled sample names
- Output directory structure remains the same as :ref:`shortreads-mode`

Example Output Structure
~~~~~~~~~~~~~~~~~~~~~~~~

::

   results/
   ├── megahit/
   │   └── sample1_dna_subsample_contigs.fa.gz
   ├── pyrodigal/
   │   ├── sample1_dna_subsample.faa.gz
   │   └── sample1_dna_subsample.gff.gz
   ├── rundbcan/
   │   └── sample1_dna_subsample_dbcan/
   └── ...

Best Practices
--------------

1. **Start Small**: Begin with a small subsample size (e.g., 5-10 million reads) to validate the pipeline
2. **Scale Up**: Gradually increase subsample size to find the optimal balance between quality and resources
3. **Compare Results**: Compare subsampled results with full dataset to assess impact on downstream analysis
4. **Resource Planning**: Use subsampling to estimate resource requirements for full dataset analysis

When to Use Subsampling
------------------------

**Recommended**:
- Initial pipeline validation
- Parameter optimization
- Resource-limited environments
- Very large datasets (>100M reads per sample)

**Not Recommended**:
- Final production analysis (use full dataset)
- Low-coverage samples (subsampling may further reduce coverage)
- When maximum sensitivity is required

See Also
--------

- :ref:`shortreads-mode` - Main short reads mode documentation
- :ref:`shortreads-coassembly` - Co-assembly mode (alternative to subsampling)
- :ref:`nextflow-parameters` - Complete parameter reference
