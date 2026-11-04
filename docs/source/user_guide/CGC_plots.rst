.. _cgc_visualization:

CAZyme Gene Cluster Visualization
===================================

Introduction
--------------

.. hint::
   The CGC visualization module is currently in beta testing. We welcome feedback and suggestions for improvement.

After identifying CAZyme Gene Clusters (CGCs), visualizing their genomic distribution can provide valuable insights into their organization and potential functional relationships.
The ``cgc_circle_plot`` module generates circular genome plots that display currently:

* Position of CGCs around the genome
* Distribution of CAZymes

These visualizations help researchers identify patterns in CGC distribution and understand their genomic context.

Basic Usage
------------

The basic command to generate CGC circular plots is (run it after ``cgc_finder``; it reads ``cgc.gff`` and
``cgc_standard_out.tsv`` from the output directory, and adds substrate labels if ``substrate_prediction.tsv`` exists):

.. code-block:: shell

   run_dbcan cgc_circle_plot --output_dir <OUTPUT_DIRECTORY>

Examples for Different Genome Types
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: shell
   :caption: Prokaryotic genome analysis

   run_dbcan cgc_circle_plot --output_dir output_EscheriaColiK12MG1655_faa

.. code-block:: shell
   :caption: Eukaryotic genome analysis

   run_dbcan cgc_circle_plot --output_dir output_Xylona_heveae_TC161_faa


Parameters
-------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Parameter
     - Description
   * - ``--output_dir``
     - Directory containing CGC data and for output files

Output Files
--------------

The visualization module generates several output files:

* ``cgc_circos/`` The output of the CGC circle plot, including the SVG files. It will be generated in the output directory. For each contig/scaffold in the genome, a separate SVG file (``cgc_circos_<contig>.svg``) is created, and there's a whole genome overview SVG file (``cgc_circos_plot.svg``).


.. tip::
   The SVG format allows for editing in vector graphics programs like ``Inkscape`` or ``Adobe Illustrator``,
   which is useful for preparing publication-quality figures.
