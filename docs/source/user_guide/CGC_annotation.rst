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

``cgc_finder`` reads ``cgc.gff``. Genes annotated as ``CAZyme`` (proteins with ``#ofTools`` >= 2 in ``overview.tsv``)
are core signature genes; genes annotated with one of the ``--additional_genes`` types (default: ``TC``; the type is
the part before ``|`` of each ``+``-separated label, matched exactly) are additional signature genes. On each contig, consecutive signature genes are joined into one candidate cluster
as long as the following rules hold:

1. **Null Gene Search** (Default, ``--use_null_genes``)

   At most ``--num_null_gene`` (default: 2) non-signature genes may lie between two consecutive signature genes.

2. **Distance-Based Search** (AntiSMASH-like, off by default, ``--use_distance``)

   The distance between two consecutive signature genes must not exceed ``--base_pair_distance`` (default: 15000 bp).

A candidate is reported as a CGC only if it has at least ``--min_cluster_genes`` genes (default: 2), at least
``--min_core_cazyme`` CAZymes (default: 1), and the required additional gene types: with ``--additional_logic all``
(default) every listed type must be present; with ``any``, at least ``--additional_min_categories`` types.
With the defaults, a CGC therefore needs at least one CAZyme and one transporter (TC).

Genes with several labels (for example ``CAZyme|GH13_31+TC|8.A.9.1.1``) are handled as follows:

* **One gene, one role.** The core CAZyme(s) and each required additional type must come from *different* genes.
  A single CAZyme+TC gene cannot be both the CAZyme and the transporter of a CGC; it can be the transporter next to
  another CAZyme, or the CAZyme next to another transporter.
* **Same region.** A domain-level label (MEROPS peptidase unit, STP or fungal TF HMM domain) that overlaps a CAZyme
  domain of the same protein by more than 50% of the shorter interval is treated as the same domain classified by
  two databases (for example CE1 and peptidase family S09 on an alpha/beta-hydrolase) and does not count as an
  additional signature; the protein counts as a CAZyme only. The label stays in ``cgc.gff`` and
  ``cgc_standard_out.tsv``. TC, TF (DIAMOND) and sulfatase labels come from full-length protein databases and are
  never compared by coordinates. E-values of different databases are not compared.

.. note::

   You can use either rule individually or combine both for stricter CGC prediction criteria.
   ``--extend_mode bp`` or ``--extend_mode gene`` (with ``--extend_bp`` / ``--extend_gene_count``) adds flanking
   genes to each reported CGC after detection.

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
   * - ``--use_null_genes`` / ``--no-use_null_genes``
     - Enable (default) or disable the null gene rule
   * - ``--num_null_gene``
     - Maximum number of consecutive non-signature genes allowed (default: 2)
   * - ``--use_distance``
     - Enable the distance-based rule (flag, off by default)
   * - ``--base_pair_distance``
     - Maximum distance (bp) between signature genes (default: 15000; only used with ``--use_distance``)
   * - ``--additional_genes``
     - Types of additional genes required (default: ``TC``). Repeat the option for several types. Types written by
       ``gff_process``: ``TC`` (transporter), ``TF`` (transcription factor), ``STP`` (signal transduction protein),
       ``Sulfatase``, ``Peptidase``
   * - ``--additional_logic`` / ``--additional_min_categories``
     - ``all`` (default): every listed type is required; ``any``: at least ``--additional_min_categories`` types
   * - ``--min_core_cazyme`` / ``--min_cluster_genes``
     - Minimum number of CAZymes (default: 1) and of genes (default: 2) in a CGC

.. hint::
   **CAZyme Gene Pairs**

   You can set ``--additional_genes`` with `CAZyme` and it will generate CGCs that include CAZyme pairs. This is useful for cases where you want to ensure that certain CAZymes are always included in the same cluster. For example, if you want to include a glycoside hydrolase (GH) and a glycosyltransferase (GT) together, you can specify `CAZyme` in the `--additional_genes` parameter.

   **Additional Gene Types**

   You can specify multiple additional gene types using the ``--additional_genes`` parameter:

   * Using multiple parameters: ``--additional_genes TC --additional_genes TF``

   **Peptidases, Sulfatases and Custom Gene Types**

   ``gff_process`` also annotates peptidases (MEROPS) and sulfatases (SulfAtlas), so ``--additional_genes Peptidase`` and
   ``--additional_genes Sulfatase`` work without extra steps. Other gene types require:

   1. Annotating these functions in your genome.
   2. Manually updating the ``CGC_annotation`` attribute in ``cgc.gff`` to include these annotations.

   **Non-coding Elements**

   While tRNAs and other non-coding genes are included in the ``cgc.gff`` file, they are not considered formal components of CAZyme Gene Clusters (considered as null genes). You can include them using the ``--additional_genes`` parameter if needed for your analysis, but this deviates from standard CGC definitions.


Output Files
-------------

The CGC finder generates several output files:

* ``cgc_standard_out.tsv`` - One row per gene in each identified CGC (CGC ID, gene type, contig, protein ID, coordinates, strand, annotation).
* ``cgc_standard_out_summary.tsv`` - One row per CGC with gene counts per type and the CGC length.

``cgc_finder`` reads ``cgc.gff`` (written by ``gff_process``) and does not modify it.

.. note::

   ``--feature_type`` is accepted but currently has no effect: all lines of ``cgc.gff`` are used (known limitation).

Annotating null genes with Pfam (optional)
-------------------------------------------

``Pfam_null_cgc`` searches Pfam-A domains with pyHMMER and adds them to genes that have no signature annotation
(``null``) inside CGCs. It needs ``Pfam-A.hmm`` in ``--db_dir``; this file is **not** downloaded by
``run_dbcan database`` (see :doc:`prepare_the_database`).

.. code-block:: shell

   run_dbcan Pfam_null_cgc --output_dir output_dir --db_dir db --null_from_gff

With ``--null_from_gff``, the proteins of all ``null`` genes in ``cgc.gff`` are written to ``null_proteins.faa``
and only these are searched; without it, all proteins in ``uniInput.faa`` are searched, which is slower. The hits
(``Pfam_hmm_results.tsv``) are then added to the null genes in new copies of the CGC files:
``cgc_standard_out.pfam_annotated.tsv`` and ``cgc.pfam_annotated.gff``.

.. admonition:: Next Steps

   After identifying CAZyme gene clusters, you can proceed to substrate prediction for those CGCs :doc:`Substrate Prediction <predict_CGC_substrate>`. This step will help you understand the potential substrates that these clusters can act upon.
   You can also proceed to :doc:`Visualizing CGCs <CGC_plots>` to create graphical visualization of your CGCs on the genome.

