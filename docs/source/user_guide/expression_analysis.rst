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
- ``plots/CGC_expression/*.expression.pdf`` — dual-panel CGC figures

Plotting
--------

Plot all CGCs (one PDF per cluster)::

    dbcan_plot CGC_expression_plot \
      -i MAG_sample.dbCAN \
      --expression-dir MAG_sample.expression

Plot a single CGC::

    dbcan_plot CGC_expression_plot \
      -i MAG_sample.dbCAN \
      --expression-dir MAG_sample.expression \
      --cgcid 'scaffold_1|CGC1'
