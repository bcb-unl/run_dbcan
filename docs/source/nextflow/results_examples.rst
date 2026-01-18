.. _nextflow-results-examples:

Nextflow Pipeline: Results Examples
====================================

This page showcases example results from the Nextflow pipeline for CAZyme annotation in microbiome data, demonstrating the output quality and visualization capabilities across different analysis modes.

Pipeline Overview
-----------------

The following figure illustrates the complete workflow of the Nextflow pipeline:

.. figure:: ../_static/img/nextflow/supplements/Figure-S1-pipeline.drawio.pdf
   :width: 100%
   :align: center
   :alt: Nextflow pipeline workflow diagram

   Figure S1: Complete workflow diagram of the Nextflow pipeline for CAZyme annotation in microbiome data. The pipeline supports three main modes: short reads assembly (MEGAHIT), long reads assembly (Flye), and assembly-free analysis (DIAMOND).

Performance Comparison
----------------------

Computational Performance
~~~~~~~~~~~~~~~~~~~~~~~~~

The following figure compares the computational requirements (CPU hours) across different analysis modes:

.. figure:: ../_static/img/nextflow/supplements/Figure-S2-CPU-hrs.pdf
   :width: 80%
   :align: center
   :alt: CPU hours comparison across analysis modes

   Figure S2: Computational performance comparison showing CPU hours required for different analysis modes. Assembly-free mode requires the least computational resources, while long reads assembly requires the most.

Results Statistics
~~~~~~~~~~~~~~~~~~

The following figure summarizes the results statistics across different modes:

.. figure:: ../_static/img/nextflow/supplements/Figure-S3-results-stats.pdf
   :width: 80%
   :align: center
   :alt: Results statistics comparison

   Figure S3: Summary statistics of CAZyme annotation results across different analysis modes, including number of CAZymes detected, families identified, and other key metrics.

Example Results by Analysis Mode
---------------------------------

Short Reads Mode - Standard Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following visualizations show example results from standard short reads analysis using two samples (Wet2014_dna and Dry2014_dna):

CAZyme Family Abundance Heatmap
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. figure:: ../_static/img/nextflow/supplements/Wet2014_dna_Dry2014_dna_pdf_main/heatmap.pdf
   :width: 100%
   :align: center
   :alt: CAZyme family abundance heatmap for standard short reads mode

   Heatmap showing CAZyme family abundance across samples in standard short reads mode.

CAZyme Family Distribution
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. figure:: ../_static/img/nextflow/supplements/Wet2014_dna_Dry2014_dna_pdf_main/fam.pdf
   :width: 100%
   :align: center
   :alt: CAZyme family distribution bar plot

   Bar plot showing the distribution of CAZyme families across samples.

CAZyme Subfamily Distribution
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. figure:: ../_static/img/nextflow/supplements/Wet2014_dna_Dry2014_dna_pdf_main/subfam.pdf
   :width: 100%
   :align: center
   :alt: CAZyme subfamily distribution bar plot

   Bar plot showing the distribution of CAZyme subfamilies across samples.

EC Number Distribution
^^^^^^^^^^^^^^^^^^^^^^

.. figure:: ../_static/img/nextflow/supplements/Wet2014_dna_Dry2014_dna_pdf_main/ec.pdf
   :width: 100%
   :align: center
   :alt: EC number distribution bar plot

   Bar plot showing the distribution of EC (Enzyme Commission) numbers across samples.

Short Reads Mode - Subsampling
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following visualizations show example results from subsampling mode:

.. figure:: ../_static/img/nextflow/supplements/Wet2014_dna_subsample_Dry2014_dna_subsample_pdf/heatmap.pdf
   :width: 100%
   :align: center
   :alt: CAZyme family abundance heatmap for subsampling mode

   Heatmap showing CAZyme family abundance in subsampling mode.

.. figure:: ../_static/img/nextflow/supplements/Wet2014_dna_subsample_Dry2014_dna_subsample_pdf/fam.pdf
   :width: 100%
   :align: center
   :alt: CAZyme family distribution for subsampling mode

   Family distribution in subsampling mode.

Short Reads Mode - Co-assembly
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following visualizations show example results from co-assembly mode:

.. figure:: ../_static/img/nextflow/supplements/coassembly_pdf/heatmap.pdf
   :width: 100%
   :align: center
   :alt: CAZyme family abundance heatmap for co-assembly mode

   Heatmap showing CAZyme family abundance in co-assembly mode.

.. figure:: ../_static/img/nextflow/supplements/coassembly_pdf/fam.pdf
   :width: 100%
   :align: center
   :alt: CAZyme family distribution for co-assembly mode

   Family distribution in co-assembly mode.

.. figure:: ../_static/img/nextflow/supplements/coassembly_pdf/subfam.pdf
   :width: 100%
   :align: center
   :alt: CAZyme subfamily distribution for co-assembly mode

   Subfamily distribution in co-assembly mode.

.. figure:: ../_static/img/nextflow/supplements/coassembly_pdf/ec.pdf
   :width: 100%
   :align: center
   :alt: EC number distribution for co-assembly mode

   EC number distribution in co-assembly mode.

Assembly-Free Mode
~~~~~~~~~~~~~~~~~~

The following visualizations show example results from assembly-free mode:

.. figure:: ../_static/img/nextflow/supplements/Dry2014_dna_Wet2014_dna_pdf_asmfree/heatmap.pdf
   :width: 100%
   :align: center
   :alt: CAZyme family abundance heatmap for assembly-free mode

   Heatmap showing CAZyme family abundance in assembly-free mode.

.. figure:: ../_static/img/nextflow/supplements/Dry2014_dna_Wet2014_dna_pdf_asmfree/fam.pdf
   :width: 100%
   :align: center
   :alt: CAZyme family distribution for assembly-free mode

   Family distribution in assembly-free mode.

.. figure:: ../_static/img/nextflow/supplements/Dry2014_dna_Wet2014_dna_pdf_asmfree/subfam.pdf
   :width: 100%
   :align: center
   :alt: CAZyme subfamily distribution for assembly-free mode

   Subfamily distribution in assembly-free mode.

.. figure:: ../_static/img/nextflow/supplements/Dry2014_dna_Wet2014_dna_pdf_asmfree/ec.pdf
   :width: 100%
   :align: center
   :alt: EC number distribution for assembly-free mode

   EC number distribution in assembly-free mode.

Interpreting the Results
-------------------------

The visualizations provided above demonstrate the comprehensive output of the Nextflow pipeline:

- **Heatmaps**: Show the abundance of CAZyme families across different samples, allowing for easy comparison and identification of differentially abundant CAZymes.

- **Family/Subfamily Bar Plots**: Display the distribution of CAZyme families and subfamilies, providing insights into the functional diversity of the microbiome.

- **EC Number Distribution**: Shows the distribution of enzyme commission numbers, indicating the functional categories of CAZymes present in the samples.

- **Performance Metrics**: The pipeline provides detailed statistics on computational requirements and result quality, helping users choose the most appropriate analysis mode for their datasets.

For more information about interpreting specific outputs, see the :ref:`nextflow-output` documentation.

See Also
--------

- :ref:`shortreads-mode` - Short reads mode documentation
- :ref:`shortreads-subsample` - Subsampling mode documentation
- :ref:`shortreads-coassembly` - Co-assembly mode documentation
- :ref:`assemfree-mode` - Assembly-free mode documentation
- :ref:`nextflow-output` - Output documentation
