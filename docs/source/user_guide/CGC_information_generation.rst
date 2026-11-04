.. _cgc-information-generation:

CAZyme Gene Cluster (CGC) Annotation
====================================

Introduction
-------------

A CAZyme Gene Cluster (CGC) refers to a group of genes co-located on the genome that collectively participate in glycan metabolism. These clusters encode various enzymes and regulatory proteins that work together to process specific carbohydrate substrates. CGCs are particularly important in microbial genomes, where they enable efficient utilization of diverse carbohydrate sources in the environment.

Generating CGC Information
---------------------------

To identify and analyze CGCs, users must first prepare annotation information by converting their GFF file into a specialized CGC-ready format. The ``gff_process`` command must be run in the output directory of a finished ``CAZyme_annotation`` run (it reads ``uniInput.faa`` and ``overview.tsv``) and performs these steps:

1. Uses DIAMOND to annotate transporters (TC, against TCDB), prokaryotic transcription factors (TF, against ``TF.dmnd``; disable with ``--no-prokaryotic``), sulfatases and peptidases
2. Uses pyHMMER to identify signal transduction proteins (STP) and, only with ``--fungi``, fungal transcription factors (``TF.hmm``)
3. Combines and filters all results based on coverage and E-value thresholds into ``total_cgc_info.tsv``
4. Marks CAZymes using ``overview.tsv``: only proteins with ``#ofTools`` >= 2 count as CAZymes, annotated with their ``Recommend Results``
5. Integrates these annotations into the user-submitted GFF file to generate a ``cgc.gff`` file

``easy_CGC`` and ``easy_substrate`` run this step automatically.

Input GFF File Types
-----------------------

``cgc.gff`` is built by matching each gene in your GFF to a protein ID in ``uniInput.faa``. Which GFF attribute is
used depends on ``--gff_type``:

.. list-table::
   :widths: 18 82
   :header-rows: 1

   * - ``--gff_type``
     - How proteins are matched
   * - ``prodigal``
     - Every feature line is used. The protein ID is the ``ID`` attribute. This is the type of ``uniInput.gff``
       written by pyrodigal in ``prok``/``meta`` mode. Features are kept in file order (not re-sorted).
   * - ``NCBI_prok``
     - Only ``gene`` features are used. The protein ID is the ``protein_id`` attribute of the ``CDS`` child of the gene.
       Genes with a non-CDS child (tRNA, rRNA, ...) are written as ``protein_id=NA`` with ``CGC_annotation=Other|<type>``.
   * - ``NCBI_euk``
     - Only ``gene`` features are used. Each gene must have ``mRNA`` children; the protein ID is the ``protein_id``
       attribute of the first ``CDS`` under the first ``mRNA`` that has one. Genes with a child that is not ``mRNA``
       (tRNA, ncRNA, a ``CDS`` directly under the gene, ...) are written as ``protein_id=NA`` with
       ``CGC_annotation=Other|<type>``.
   * - ``JGI``
     - Only ``gene`` features are used. The protein ID is the ``proteinId`` attribute of the gene, compared with the
       third ``-``-separated field of the FASTA IDs. JGI headers such as ``jgi|Xylhe1|123456|name`` become
       ``jgi-Xylhe1-123456-name`` in protein mode, so ``proteinId=123456`` matches.

For ``NCBI_prok``, ``NCBI_euk`` and ``JGI`` the features of each contig are sorted by position, and ``cgc.gff``
contains one line per gene. Only ``prodigal``, ``NCBI_prok``, ``NCBI_euk`` and ``JGI`` are accepted.

.. warning::

   **Unmatched genes fail silently.** If the expected attribute is missing, the gene is written as
   ``protein_id=unknown``; if the ID does not occur in ``overview.tsv`` or ``total_cgc_info.tsv``, it is written as
   ``CGC_annotation=null``. No warning is printed. When no gene matches, ``cgc_finder`` finds no CGC and
   ``cgc_standard_out.tsv`` is empty. Check ``cgc.gff`` after this step:

   .. code-block:: shell

      grep -c "protein_id=unknown" output_dir/cgc.gff
      grep -c "CGC_annotation=CAZyme" output_dir/cgc.gff

   If most genes are ``unknown`` or no line contains ``CAZyme``, the GFF does not use the IDs of your protein FASTA.

Eukaryotic genomes (fungi, etc.)
---------------------------------

There is no eukaryotic genome mode. For fungi and other eukaryotes:

1. Predict genes with a eukaryotic gene finder (e.g. funannotate or BRAKER), or use an NCBI/JGI annotation.
2. Run ``CAZyme_annotation`` (or ``easy_CGC``) with ``--mode protein`` on the predicted proteins.
3. Provide the matching GFF with ``--input_gff`` and ``--gff_type NCBI_euk`` (or ``JGI`` for MycoCosm files).
   Add ``--fungi`` to search fungal transcription factors with ``TF.hmm``.

GFF3 files from funannotate, BRAKER and similar pipelines usually have ``gene`` → ``mRNA`` → ``CDS`` features, but
their ``CDS`` rows have no ``protein_id`` attribute, so every gene becomes ``protein_id=unknown``. Their protein
FASTA IDs are usually the transcript (``mRNA``) IDs. In that case, add ``protein_id=<FASTA ID>`` to each ``CDS``
row, taken from its ``Parent`` attribute, and use ``--gff_type NCBI_euk``:

.. code-block:: shell

   awk 'BEGIN { FS = OFS = "\t" }
   $3 == "CDS" && $9 !~ /protein_id=/ && match($9, /Parent=[^;]+/) {
       pid = substr($9, RSTART + 7, RLENGTH - 7)
       sub(/;$/, "", $9)
       $9 = $9 ";protein_id=" pid
   }
   { print }' annotation.gff3 > annotation.protein_id.gff3

   run_dbcan easy_CGC --mode protein --input_raw_data proteins.faa \
     --input_gff annotation.protein_id.gff3 --gff_type NCBI_euk \
     --output_dir output_fungus --db_dir db --fungi

Before running, compare a few FASTA IDs with the new ``protein_id`` values (they must be identical, at most
80 characters and unique; see :doc:`CAZyme_annotation`). If the file uses ``transcript`` instead of ``mRNA``
in column 3, rename those features to ``mRNA``. A GTF file must be converted to GFF3 with ``ID``/``Parent``
attributes first.

Command Examples
------------------

.. code-block:: shell

   run_dbcan gff_process --output_dir output_dir --db_dir db --input_gff EscheriaColiK12MG1655.gff --gff_type NCBI_prok

.. code-block:: shell

   run_dbcan gff_process --output_dir output_dir --db_dir db --input_gff output_dir/uniInput.gff --gff_type prodigal

.. code-block:: shell

   run_dbcan gff_process --output_dir output_dir --db_dir db --input_gff Xylhe1_GeneCatalog_proteins.gff --gff_type JGI

.. code-block:: shell

   run_dbcan gff_process --output_dir output_dir --db_dir db --input_gff eukaryotic_genome.gff --gff_type NCBI_euk --fungi

.. note::

   ``gff_process`` has no ``--mode`` option. If ``--input_gff`` is omitted, it uses ``<output_dir>/uniInput.gff`` with
   ``--gff_type prodigal``, which only exists after a ``prok``/``meta`` run. For protein input always pass
   ``--input_gff`` and ``--gff_type``. If the GFF file is missing, an error is logged but no ``cgc.gff`` is written.

.. note::

   Eukaryotic GFF processing is currently in beta. While functional, we're still validating and improving this feature.

Key Parameters
~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - Parameter
     - Description
   * - ``--output_dir``
     - Directory for both input and output files, it is generated from the ``CAZyme_annotation`` step
   * - ``--db_dir``
     - Directory containing database files, it is generated from the ``database`` step
   * - ``--input_gff``
     - Path to your GFF annotation file (default: ``<output_dir>/uniInput.gff``)
   * - ``--gff_type``
     - Format of the GFF file: ``NCBI_prok``, ``prodigal``, ``JGI``, or ``NCBI_euk`` (default: ``prodigal``)
   * - ``--fungi`` / ``--no-fungi``
     - Run the fungal TF HMM search (default: off)
   * - ``--prokaryotic`` / ``--no-prokaryotic``
     - Run the prokaryotic TF DIAMOND search (default: on)

E-value and coverage thresholds for each search (``--e_value_threshold_tc``, ``--coverage_threshold_tc``, ...)
are listed by ``run_dbcan gff_process --help``. DIAMOND coverage options are percentages of the query length (0-100);
pyHMMER coverage options are fractions of the HMM length (0-1).

Output Files
-------------

This step generates several output files:

* ``diamond.out.tc`` - DIAMOND output for Transporters (TCs) via TCDB database (https://www.tcdb.org/)
* ``diamond.out.tf`` - DIAMOND output for prokaryotic Transcription Factors (TFs)
* ``TF_hmm_results.tsv`` - pyHMMER output for fungal Transcription Factors (only with ``--fungi``)
* ``STP_hmm_results.tsv`` - pyHMMER output for Signal Transduction Proteins (STPs) using the STP-HMM database
* ``diamond.out.sulfatase`` and ``diamond.out.peptidase`` - DIAMOND output for sulfatases and peptidases
* ``total_cgc_info.tsv`` - Combined table of all non-CAZyme signature genes (TC, TF, STP, Sulfatase, Peptidase)
* ``cgc.gff`` - Simplified GFF with ``protein_id`` and ``CGC_annotation`` attributes, the input of ``cgc_finder``

.. admonition:: Next Steps

   After generating the CGC annotation information, proceed to :doc:`CGC identification <CGC_annotation>` to identify and analyze CAZyme gene clusters in your genome.
