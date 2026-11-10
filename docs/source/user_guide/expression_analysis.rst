Expression and DEG Analysis
============================

The ``run_dbcan expression`` command automates read mapping (optional), gene-level
counting, CAZyme/CGC abundance estimation, PyDESeq2 differential expression, and
CGC expression visualization.

Requirements
------------

Install the optional expression dependencies::

    pip install 'dbcan[expression]'

This pins ``numpy>=2.0,<2.4`` for compatibility with ``anndata`` (used by PyDESeq2).
If you see ``numpy.dtypes has no attribute 'StringDType'``, downgrade NumPy::

    pip install 'numpy>=2.0,<2.4'

External tools (conda environment ``CAZyme_annotation``):

- ``bwa``, ``samtools`` (when aligning FASTQ)
- Existing ``run_dbcan`` outputs on a **single reference assembly**

Samplesheet
-----------

Tab-separated file with header::

    sample_id	condition	bam	r1	r2	replicate_group
    SRF_MT	control	/path/SRF.bam		
    BML_MT	treatment		BML_1.fq.gz	BML_2.fq.gz

Each row must provide ``bam`` **or** ``r1`` (and optionally ``r2``).

**DESeq2 replicates:** PyDESeq2 needs at least **two samples per condition** (e.g. two
``control`` and two ``treatment`` rows). One sample per condition label cannot estimate
dispersion. To compare four libraries as control vs treatment, assign the same
``condition`` to replicates instead of unique labels (``test1``, ``test2``, …) per sample.
Use ``--skip-deseq2`` if you only need abundance tables and CGC plots.

GFF input
---------

``--gff`` is optional. If omitted, the pipeline selects (in order):

1. ``{input_dir}/cgc.gff`` — **recommended** after ``run_dbcan``; uses ``protein_id=`` aligned with ``cgc_standard_out.tsv``
2. ``{input_dir}/*.fix.gff`` — from ``dbcan_utils gff_fix`` (Prodigal ``ID=``)
3. ``{input_dir}/uniInput.gff``

Both ``cgc.gff`` and ``*.fix.gff`` are supported for read counting.

Reads are counted per gene with the ``dbcan_utils cal_coverage`` rules and defaults (see :doc:`abundance`):
the two mates of a pair count once (fragments), counting is unstranded, secondary/supplementary/QC-fail/duplicate
alignments are skipped, genes are keyed by ``protein_id`` (else ``ID``, ``Name``, ``locus_tag``) and rows that
share an ID are one gene. ``--overlap-base-ratio``, ``--mapping-quality``, ``--identity`` and ``--hifi`` are
passed on.

**Eukaryotes (multi-exon genes):** ``cgc.gff`` from NCBI eukaryotic annotation stores **one row per gene**
(full locus from transcription start to end, not one row per CDS exon). Read counts and CGC plots therefore
represent the **whole gene**; introns are inside the arrow span. This matches ``cgc_standard_out.tsv`` and
is the intended behavior for CGC-level expression figures (one point per gene, not per exon).

Example
-------

.. code-block:: shell

    run_dbcan expression \
      -i MAG_sample.dbCAN \
      --reference-fasta MAG_sample.contigs.fa \
      --samplesheet samples.tsv \
      --output-dir MAG_sample.expression \
      --run-plots

``cgc.gff`` under ``MAG_sample.dbCAN/`` is used automatically when ``--gff`` is omitted.

Outputs
-------

- ``count_matrix.tsv`` — raw gene counts per sample
- ``gene_deg.tsv``, ``DEG.tsv`` — gene-level DESeq2 results
- ``cgc_deg.tsv``, ``cgc_gene_deg_map.tsv`` — CGC-level summaries
- ``cgc_gene_counts.tsv`` — long table for plotting
- ``gene_norm_counts.tsv`` — DESeq2 normalized counts (after DE analysis)
- ``plots/CGC_expression/*.expression.pdf`` — SeMa-Trap style dual-panel figures (step expression + gene map)
- ``plots/CGC_expression/*.gene_heatmap.pdf`` — per-CGC gene × sample/condition heatmap

Plotting
--------

CGC expression figures use a **gene-level expression track** (SeMa-Trap style) in the top
panel: for each gene, a **narrow rectangle** spans its genomic start–end (same x as arrows
below); the rectangle's **vertical position** is the expression value and its **height is
fixed**. Adjacent genes are joined by thin lines between rectangle edge centers, forming one
colored **trajectory** per ``condition``. The bottom panel shows gene arrows with functional
labels and DEG markers. Default metric is ``log2(TPM+1)``.

Plot all CGCs (one PDF + heatmap per cluster)::

    dbcan_plot CGC_expression_plot \
      -i MAG_sample.dbCAN \
      --expression-dir MAG_sample.expression \
      --metric log2_tpm \
      --heatmap-rows sample \
      --deg-marker both

The other expression figures::

    dbcan_plot DEG_volcano -i MAG_sample.dbCAN --expression-dir MAG_sample.expression
    dbcan_plot DEG_CGC_bar -i MAG_sample.dbCAN --expression-dir MAG_sample.expression
    dbcan_plot CGC_expression_heatmap -i MAG_sample.dbCAN --expression-dir MAG_sample.expression --top 50

They write ``plots/DEG_volcano.pdf``, ``plots/DEG_CGC_bar.pdf`` and ``plots/CGC_expression_heatmap.pdf`` under
``--expression-dir`` (``--pdf`` sets another path) and exit with status 1 when their input table is missing or empty.

Plot a single CGC::

    dbcan_plot CGC_expression_plot \
      -i MAG_sample.dbCAN \
      --expression-dir MAG_sample.expression \
      --cgcid 'scaffold_1|CGC1' \
      --metric log2fc

``--plot-metric`` / ``dbcan_plot --metric`` choices:

- **Recommended for expression tracks (log-transformed):** ``log2_tpm`` (default), ``log2_fpkm``,
  ``log2_rpkm``, ``log2_rpm``, ``log2_cpm``, ``log2_norm``
- **Differential / comparative:** ``log2fc`` (needs DESeq2 + replicates), ``gene_zscore`` (within-CGC
  z-score per condition)
- **Linear scale (optional):** ``tpm``, ``fpkm``, ``rpkm``, ``rpm``, ``counts``

``cgc_gene_counts.tsv`` stores ``count``, ``tpm``, ``fpkm``, ``rpm``, and ``rpkm`` per gene for plotting.
Re-run quantification after upgrading if older outputs lack the new columns.
``--heatmap-rows``: ``sample`` (default) or ``condition``.
