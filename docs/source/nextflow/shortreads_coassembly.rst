.. _shortreads-coassembly:

Short Reads: Co-assembly Mode
===============================

The co-assembly mode combines reads from all samples and performs a single joint assembly, improving contig continuity and enhancing detection of shared genomic features across samples.

Overview
--------

Co-assembly is particularly useful when:

- You want to improve assembly quality by combining sequencing depth from multiple samples
- You're analyzing samples from similar environments or conditions
- You want to detect shared CAZyme gene clusters across samples
- You need longer, more complete contigs for better gene prediction

In co-assembly mode, all reads from all samples are combined (preserving paired-end or single-end structure) and assembled together using MEGAHIT. The resulting assembly is then used for all downstream analysis, but read mapping and abundance calculation are performed per-sample.

Activation
----------

.. note::
   Before running the pipeline, make sure you have cloned the repository. See :ref:`nextflow-usage` for installation instructions.

To enable co-assembly mode, use the ``--coassembly`` flag along with ``--type shortreads`` (assuming you are in the ``dbcan-nf`` directory):

.. code-block:: bash

   nextflow run main.nf \
     --input samplesheet.csv \
     --outdir results \
     --type shortreads \
     --coassembly \
     -profile docker

Requirements
------------

- **Minimum Samples**: At least 2 samples are required. The pipeline will produce an error if fewer samples are provided.
- **Compatible Data**: All samples should be from similar sequencing runs or conditions for best results
- **Mutually Exclusive**: Cannot be used together with ``--subsample``

Parameters
----------

.. list-table:: Co-assembly Parameters
   :widths: 20 10 10 60
   :header-rows: 1

   * - Parameter
     - Type
     - Default
     - Description
   * - ``--coassembly``
     - boolean
     - ``false``
     - Enable co-assembly mode. Must be used with ``--type shortreads`` and requires at least 2 samples.

Usage Examples
--------------

Basic Co-assembly
~~~~~~~~~~~~~~~~~

Co-assemble all samples in the samplesheet:

.. code-block:: bash

   nextflow run main.nf \
     --input samplesheet.csv \
     --outdir results \
     --type shortreads \
     --coassembly \
     -profile docker

Co-assembly with Multiple Samples
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The samplesheet should contain at least 2 samples:

.. code-block:: text

   sample,fastq_1,fastq_2
   CONTROL_REP1,control1_R1.fastq.gz,control1_R2.fastq.gz
   CONTROL_REP2,control2_R1.fastq.gz,control2_R2.fastq.gz
   TREATMENT_REP1,treatment1_R1.fastq.gz,treatment1_R2.fastq.gz
   TREATMENT_REP2,treatment2_R1.fastq.gz,treatment2_R2.fastq.gz

Co-assembly with Other Options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Combine co-assembly with other parameters:

.. code-block:: bash

   nextflow run main.nf \
     --input samplesheet.csv \
     --outdir results \
     --type shortreads \
     --coassembly \
     --skip_kraken_extraction \
     -profile docker

Workflow Behavior
-----------------

Assembly Phase
~~~~~~~~~~~~~~

1. **Read Combination**: All reads from all samples are combined into a single set
2. **Single Assembly**: One MEGAHIT assembly is performed on the combined reads
3. **Assembly Naming**: The co-assembly is named ``coassembly`` in intermediate files

Annotation Phase
~~~~~~~~~~~~~~~~

1. **Single Annotation**: CAZyme annotation is performed once on the co-assembly
2. **Result Replication**: Annotation results are replicated to each original sample for downstream processing

Read Mapping Phase
~~~~~~~~~~~~~~~~~~

1. **Per-Sample Mapping**: Each sample's reads are mapped back to the co-assembly contigs
2. **Sample-Specific Coverage**: Coverage and abundance are calculated per-sample
3. **Preserved Sample Identity**: All output files maintain original sample names

Output Files
------------

Output Structure
~~~~~~~~~~~~~~~~

The co-assembly mode produces output files with a specific structure:

- **Co-assembly Files**: Assembly and annotation files use ``coassembly`` as the sample name
- **Sample-Specific Files**: Read mapping, coverage, and abundance files use original sample names
- **Replicated Annotations**: CAZyme annotation results are available for each sample

Example Output Structure
~~~~~~~~~~~~~~~~~~~~~~~~

::

   results/
   ├── megahit/
   │   └── coassembly_contigs.fa.gz          # Single co-assembly
   ├── pyrodigal/
   │   ├── coassembly.faa.gz                 # Genes from co-assembly
   │   └── coassembly.gff.gz
   ├── rundbcan/
   │   └── coassembly_dbcan/                 # Annotation from co-assembly
   ├── bwa_index_mem/
   │   ├── sample1_dna.bam                   # Per-sample mapping
   │   ├── sample1_dna.bam.bai
   │   ├── sample2_dna.bam
   │   └── sample2_dna.bam.bai
   ├── dbcan_utils_cal_abund/
   │   ├── sample1_dna_abund/                # Per-sample abundance
   │   └── sample2_dna_abund/
   └── ...

Key Differences from Standard Mode
-----------------------------------

1. **Single Assembly**: One assembly instead of per-sample assemblies
2. **Shared Contigs**: All samples share the same contig set
3. **Per-Sample Abundance**: Abundance is still calculated per-sample
4. **No RNA-seq**: RNA-seq processing is automatically disabled

Advantages
----------

- **Improved Contiguity**: Longer, more complete contigs due to increased sequencing depth
- **Better Gene Detection**: Shared genes are more likely to be detected and fully assembled
- **Reduced Fragmentation**: Fewer fragmented genes and gene clusters
- **Computational Efficiency**: Single assembly is more efficient than multiple per-sample assemblies

Considerations
--------------

- **Sample Compatibility**: Best results when samples are from similar environments
- **Heterogeneity**: Highly diverse samples may produce less optimal co-assemblies
- **Abundance Comparison**: Abundance values are comparable across samples as they use the same reference
- **Memory Requirements**: Co-assembly may require more memory than per-sample assembly

Best Practices
--------------

1. **Sample Selection**: Use samples from similar conditions or environments
2. **Quality Control**: Ensure all samples have similar quality before co-assembly
3. **Sample Size**: 2-10 samples typically work well; very large numbers may be computationally intensive
4. **Validation**: Compare co-assembly results with per-sample assemblies to assess improvement

When to Use Co-assembly
------------------------

**Recommended**:
- Samples from similar environments or conditions
- When improved contig continuity is important
- When detecting shared features across samples
- When computational resources allow single large assembly

**Not Recommended**:
- Highly diverse or unrelated samples
- When sample-specific assemblies are required
- When RNA-seq analysis is needed (automatically disabled)
- Single sample analysis (requires at least 2 samples)

Example Results
---------------

For example visualizations from co-assembly mode, see the :ref:`co-assembly results section <nextflow-results-examples>`.

See Also
--------

- :ref:`shortreads-mode` - Main short reads mode documentation
- :ref:`shortreads-subsample` - Subsampling mode (alternative to co-assembly)
- :ref:`nextflow-parameters` - Complete parameter reference
- :ref:`nextflow-results-examples` - Example results and visualizations