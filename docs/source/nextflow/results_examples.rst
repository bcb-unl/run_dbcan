.. _nextflow-results-examples:

Nextflow Pipeline: Results Examples
====================================

This page showcases example results from the Nextflow pipeline for CAZyme annotation in microbiome data, demonstrating the output quality and visualization capabilities across different analysis modes.

Pipeline Overview
-----------------

The following figure illustrates the complete workflow of the Nextflow pipeline:

.. figure:: ../_static/img/nextflow/supplements/Figure-S1-pipeline.drawio.png
   :width: 100%
   :align: center
   :alt: Nextflow pipeline workflow diagram

   Figure S1: Complete workflow diagram of the Nextflow pipeline for CAZyme annotation in microbiome data. The pipeline supports three main modes: short reads assembly (MEGAHIT), long reads assembly (Flye), and assembly-free analysis (DIAMOND). `Download PDF version <../_static/img/nextflow/supplements/Figure-S1-pipeline.drawio.pdf>`_.

Performance Comparison
----------------------

Computational Performance
~~~~~~~~~~~~~~~~~~~~~~~~~

The following figure compares the computational requirements (CPU hours) across different analysis modes:

.. figure:: ../_static/img/nextflow/supplements/Figure-S2-CPU-hrs.png
   :width: 80%
   :align: center
   :alt: CPU hours comparison across analysis modes

   Figure S2: Computational performance comparison showing CPU hours required for different analysis modes. Assembly-free mode requires the least computational resources, while long reads assembly requires the most. `Download PDF version <../_static/img/nextflow/supplements/Figure-S2-CPU-hrs.pdf>`_.

Results Statistics
~~~~~~~~~~~~~~~~~~

The following figure summarizes the results statistics across different modes:

.. figure:: ../_static/img/nextflow/supplements/Figure-S3-results-stats.png
   :width: 80%
   :align: center
   :alt: Results statistics comparison

   Figure S3: Summary statistics of CAZyme annotation results across different analysis modes, including number of CAZymes detected, families identified, and other key metrics. `Download PDF version <../_static/img/nextflow/supplements/Figure-S3-results-stats.pdf>`_.

Example Results by Analysis Mode
---------------------------------

Short Reads Mode - Standard Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following visualizations show example results from standard short reads analysis using two samples (Wet2014_dna and Dry2014_dna):

CAZyme Family Abundance Heatmap
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. figure:: ../_static/img/nextflow/supplements/Wet2014_dna_Dry2014_dna_png_main/heatmap.png
   :width: 100%
   :align: center
   :alt: CAZyme family abundance heatmap for standard short reads mode

   Heatmap showing CAZyme family abundance across samples in standard short reads mode. `Download PDF version <../_static/img/nextflow/supplements/Wet2014_dna_Dry2014_dna_pdf_main/heatmap.pdf>`_.

CAZyme Family Distribution
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. figure:: ../_static/img/nextflow/supplements/Wet2014_dna_Dry2014_dna_png_main/fam.png
   :width: 100%
   :align: center
   :alt: CAZyme family distribution bar plot

   Bar plot showing the distribution of CAZyme families across samples. `Download PDF version <../_static/img/nextflow/supplements/Wet2014_dna_Dry2014_dna_pdf_main/fam.pdf>`_.

CAZyme Subfamily Distribution
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. figure:: ../_static/img/nextflow/supplements/Wet2014_dna_Dry2014_dna_png_main/subfam.png
   :width: 100%
   :align: center
   :alt: CAZyme subfamily distribution bar plot

   Bar plot showing the distribution of CAZyme subfamilies across samples. `Download PDF version <../_static/img/nextflow/supplements/Wet2014_dna_Dry2014_dna_pdf_main/subfam.pdf>`_.

EC Number Distribution
^^^^^^^^^^^^^^^^^^^^^^

.. figure:: ../_static/img/nextflow/supplements/Wet2014_dna_Dry2014_dna_png_main/ec.png
   :width: 100%
   :align: center
   :alt: EC number distribution bar plot

   Bar plot showing the distribution of EC (Enzyme Commission) numbers across samples. `Download PDF version <../_static/img/nextflow/supplements/Wet2014_dna_Dry2014_dna_pdf_main/ec.pdf>`_.

Short Reads Mode - Subsampling
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following visualizations show example results from subsampling mode:

.. figure:: ../_static/img/nextflow/supplements/Wet2014_dna_subsample_Dry2014_dna_subsample_png/heatmap.png
   :width: 100%
   :align: center
   :alt: CAZyme family abundance heatmap for subsampling mode

   Heatmap showing CAZyme family abundance in subsampling mode. `Download PDF version <../_static/img/nextflow/supplements/Wet2014_dna_subsample_Dry2014_dna_subsample_pdf/heatmap.pdf>`_.

.. figure:: ../_static/img/nextflow/supplements/Wet2014_dna_subsample_Dry2014_dna_subsample_png/fam.png
   :width: 100%
   :align: center
   :alt: CAZyme family distribution for subsampling mode

   Family distribution in subsampling mode. `Download PDF version <../_static/img/nextflow/supplements/Wet2014_dna_subsample_Dry2014_dna_subsample_pdf/fam.pdf>`_.

Short Reads Mode - Co-assembly
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following visualizations show example results from co-assembly mode:

.. figure:: ../_static/img/nextflow/supplements/coassembly_png/heatmap.png
   :width: 100%
   :align: center
   :alt: CAZyme family abundance heatmap for co-assembly mode

   Heatmap showing CAZyme family abundance in co-assembly mode. `Download PDF version <../_static/img/nextflow/supplements/coassembly_pdf/heatmap.pdf>`_.

.. figure:: ../_static/img/nextflow/supplements/coassembly_png/fam.png
   :width: 100%
   :align: center
   :alt: CAZyme family distribution for co-assembly mode

   Family distribution in co-assembly mode. `Download PDF version <../_static/img/nextflow/supplements/coassembly_pdf/fam.pdf>`_.

.. figure:: ../_static/img/nextflow/supplements/coassembly_png/subfam.png
   :width: 100%
   :align: center
   :alt: CAZyme subfamily distribution for co-assembly mode

   Subfamily distribution in co-assembly mode. `Download PDF version <../_static/img/nextflow/supplements/coassembly_pdf/subfam.pdf>`_.

.. figure:: ../_static/img/nextflow/supplements/coassembly_png/ec.png
   :width: 100%
   :align: center
   :alt: EC number distribution for co-assembly mode

   EC number distribution in co-assembly mode. `Download PDF version <../_static/img/nextflow/supplements/coassembly_pdf/ec.pdf>`_.

Assembly-Free Mode
~~~~~~~~~~~~~~~~~~

The following visualizations show example results from assembly-free mode:

.. figure:: ../_static/img/nextflow/supplements/Dry2014_dna_Wet2014_dna_png_asmfree/heatmap.png
   :width: 100%
   :align: center
   :alt: CAZyme family abundance heatmap for assembly-free mode

   Heatmap showing CAZyme family abundance in assembly-free mode. `Download PDF version <../_static/img/nextflow/supplements/Dry2014_dna_Wet2014_dna_pdf_asmfree/heatmap.pdf>`_.

.. figure:: ../_static/img/nextflow/supplements/Dry2014_dna_Wet2014_dna_png_asmfree/fam.png
   :width: 100%
   :align: center
   :alt: CAZyme family distribution for assembly-free mode

   Family distribution in assembly-free mode. `Download PDF version <../_static/img/nextflow/supplements/Dry2014_dna_Wet2014_dna_pdf_asmfree/fam.pdf>`_.

.. figure:: ../_static/img/nextflow/supplements/Dry2014_dna_Wet2014_dna_png_asmfree/subfam.png
   :width: 100%
   :align: center
   :alt: CAZyme subfamily distribution for assembly-free mode

   Subfamily distribution in assembly-free mode. `Download PDF version <../_static/img/nextflow/supplements/Dry2014_dna_Wet2014_dna_pdf_asmfree/subfam.pdf>`_.

.. figure:: ../_static/img/nextflow/supplements/Dry2014_dna_Wet2014_dna_png_asmfree/ec.png
   :width: 100%
   :align: center
   :alt: EC number distribution for assembly-free mode

   EC number distribution in assembly-free mode. `Download PDF version <../_static/img/nextflow/supplements/Dry2014_dna_Wet2014_dna_pdf_asmfree/ec.pdf>`_.

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
