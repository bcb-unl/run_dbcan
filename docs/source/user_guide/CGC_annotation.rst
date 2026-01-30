.. _cgc-finder:

CAZyme Gene Cluster (CGC) Identification
========================================

Introduction
---------------

After preparing the CGC annotation information and generating a ``cgc.gff`` file, the next step is to identify CAZyme Gene Clusters (CGCs) in your genome. The ``cgc_finder`` command analyzes the annotated genes to detect clusters involved in carbohydrate metabolism.

Basic Usage
--------------

Once the ``cgc.gff`` file is created (automatically saved in your output directory), you can run the CGC finder:

.. code-block:: shell
   :caption: Basic CGC Finder Usage
   :emphasize-lines: 1

   run_dbcan cgc_finder --output_dir <OUTPUT_DIRECTORY>

Examples for Different Genome Types
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: shell
   :caption: Prokaryotic protein analysis

   run_dbcan cgc_finder --output_dir output_EscheriaColiK12MG1655_faa

.. code-block:: shell
   :caption: Prokaryotic genome analysis

   run_dbcan cgc_finder --output_dir output_EscheriaColiK12MG1655_fna/

.. code-block:: shell
   :caption: Fungal protein analysis

   run_dbcan cgc_finder --output_dir output_Xylona_heveae_TC161_faa/

.. code-block:: shell
   :caption: JGI protein catalog analysis

   run_dbcan cgc_finder --output_dir output_Xylhe1_faa/

CGC Prediction Rules
------------------------

run_dbCAN supports two complementary rules for predicting CGCs:

1. **Null Gene Search** (Default)

   Forward and backward search with a defined number of non-significant genes. When a core/additional gene is found, the search extends to the next iteration.

2. **Distance-Based Search** (AntiSMASH-like)

   Uses base-pair distance (default: 15kb) to search forward and backward for core/additional genes. The distance is measured between consecutive significant genes.

.. note::

   You can use either rule individually or combine both for stricter CGC prediction criteria.

Advanced Usage
------------------

To customize CGC prediction parameters:

.. code-block:: shell
   :caption: Advanced CGC prediction with custom parameters
   :emphasize-lines: 1

   run_dbcan cgc_finder --output_dir output_dir --use_null_genes --num_null_gene 5 --use_distance --base_pair_distance 15000 --additional_genes TC --additional_genes TF --additional_genes STP

Key Parameters
~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Parameter
     - Description
   * - ``--output_dir``
     - Directory containing the ``cgc.gff`` file and for output files. It is generated from the ``CGC_information_generation`` step.
   * - ``--use_null_genes``
     - Enable null gene search strategy (true/false)
   * - ``--num_null_gene``
     - Maximum number of consecutive non-significant genes allowed
   * - ``--use_distance``
     - Enable distance-based search strategy (true/false)
   * - ``--base_pair_distance``
     - Maximum distance (bp) between significant genes
   * - ``--additional_genes``
     - Types of additional genes to include (TC: Transporter, TF: Transcription Factor, STP: Signal Transduction Protein)

.. hint::
   **CAZyme Gene Pairs**

   You can set ``--additional_genes`` with `CAZyme` and it will generate CGCs that include CAZyme pairs. This is useful for cases where you want to ensure that certain CAZymes are always included in the same cluster. For example, if you want to include a glycoside hydrolase (GH) and a glycosyltransferase (GT) together, you can specify `CAZyme` in the `--additional_genes` parameter.

   **Additional Gene Types**

   You can specify multiple additional gene types using the ``--additional_genes`` parameter:

   * Using multiple parameters: ``--additional_genes TC --additional_genes TF``

   **Custom Gene Types**

   Beyond the standard types (TC, TF, STP), you can include custom gene types such as peptidases or sulfatases. However, this requires:

   1. Annotating these functions in your genome.
   2. Manually updating the ``cgc.gff`` file to include these annotations.

   **Non-coding Elements**

   While tRNAs and other non-coding genes are included in the ``cgc.gff`` file, they are not considered formal components of CAZyme Gene Clusters (considered as null genes). You can include them using the ``--additional_genes`` parameter if needed for your analysis, but this deviates from standard CGC definitions.


Output Files
-------------

The CGC finder generates several output files:

* ``cgc_standard.tsv`` - Text file listing all identified CGCs and their components.
* ``cgc.gff`` - GFF-like format annotated with functional genes for CGCFinder and visualization.

.. admonition:: Next Steps

   After identifying CAZyme gene clusters, you can proceed to substrate prediction for those CGCs :doc:`Substrate Prediction <predict_CGC_substrate>`. This step will help you understand the potential substrates that these clusters can act upon.
   You can also proceed to :doc:`Visualizing CGCs <CGC_plots>` to create graphical visualization of your CGCs on the genome.

