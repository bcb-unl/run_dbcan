Unreleased
~~~~~~~~~~

.. rubric:: Changed defaults and command lines

#. ``--mode`` accepts only ``prok``, ``meta`` or ``protein``. pyrodigal (``prok``/``meta``) is for bacteria, archaea
   and metagenomes; for fungi and other eukaryotes predict genes with funannotate, BRAKER or similar and run
   ``--mode protein`` with their GFF. A warning is logged when pyrodigal proteins look like a split eukaryotic genome.
#. ``--mode protein`` (and any user ``--input_gff``) requires ``--gff_type``. New types: ``funannotate``, ``BRAKER``,
   ``GFF3``, ``GTF`` (case-insensitive). The gene/protein match is checked before step 1 and the run stops when fewer
   than ``--gff_min_match`` (default 0.9) of the protein-coding genes match; the report is ``gff_resolution.tsv``.
#. ``--prokaryotic`` / ``--fungi``: TF annotation uses TF DIAMOND for prokaryotes (default) and the TF HMMs with
   ``--fungi``.
#. ``substrate_prediction``: ``-cpn`` now takes effect; ``-uqcgn``/``-upghn`` compare the right side; mismatched
   ``-ept``/``-eptn`` exit with status 2; ``-iden`` is a percent identity (0-100); ``-odbcanpul`` defaults to
   ``false``; ``--pul``, ``-o``, ``-w``, ``-rerun`` and ``-env`` are hidden and ignored with a warning; the legacy
   ``--mode``/``--input_raw_data`` are accepted silently. New ``--synteny_revcomp/--no-synteny_revcomp`` (default on).
#. ``dbcan_utils cal_coverage`` counts **fragments** by default (``--count-mode fragments|reads``): paired-end raw
   counts are about half of 5.2.9; RPM/TPM/FPKM are unaffected. Secondary, supplementary, QC-fail and duplicate
   alignments are skipped. New ``--id-attr``, ``--feature-type``, ``--strandedness``, ``--include-duplicates``.
   ``*_abund`` write into ``-o``; ``-a RPKM`` is accepted as FPKM.
#. ``--run_pfam``, ``--force_topology`` and ``--feature_type`` are hidden and ignored with a warning.
#. Failed steps exit with a non-zero status (missing/mismatched GFF, CGCFinder, overview, SignalP/DeepTMHMM,
   plots, synteny figures, abundance ID mismatches).

.. rubric:: Changed outputs

#. ``overview.tsv``: gene IDs are always read as text (numeric or ``NA`` IDs are no longer lost); the ``DIAMOND``
   column keeps only CAZy family tokens, once each; an EC without a count is written ``EC:1``; the ``Substrate``
   column looks the EC up across the eHMM's composition families, falls back to the family-level substrate, and uses
   ``dbCAN-sub-substrate.tsv`` from the database directory when present (same result as the dbCAN-sub web server).
   Only the tools selected with ``--methods`` are merged.
#. dbCAN-sub subfamily names are no longer rewritten by the ``GT2_*`` rename (it applies to ``dbCAN.hmm`` only).
#. CGC signatures: TC, TF (DIAMOND), Sulfatase and Peptidase are full-length DIAMOND searches and keep the best hit
   per protein without a coordinate overlap test against dbCAN domains. Peptidase hits are filtered by subject
   (MEROPS unit) coverage, and MEROPS inhibitor units are dropped. Domain-level sources (STP, fungal TF HMM,
   Peptidase) that cover the same region as a CAZyme domain are recorded as additional roles. ``total_cgc_info.tsv``
   has a ``Source`` column. A CGC is accepted only when distinct genes fill its required roles (one gene, one role).
#. ``--extend_mode`` extensions stop before the core of a neighbouring CGC.
#. ``CGC.faa`` headers (``<contig>|<CGC#>|<protein>|<Gene Type>|<family>``): the family is now that of the Gene
   Type role of a multi-role gene (CAZyme > TC > TF > STP > Sulfatase > Peptidase), not the second ``|`` field of
   the whole annotation. For example ``STP|LacI+TF|CCPA_PRIMG+STP|...`` with Gene Type ``TF`` gives ``CCPA_PRIMG``,
   and ``TF|LACI_ECOLI+STP|LacI+TC|3.A.1.2.1+...`` with Gene Type ``TC`` gives ``3.A.1``. On *E. coli* K-12 57
   headers change; ``substrate_prediction.tsv`` is identical after sorting (scoring uses the Gene Type only). Scripts
   that parse the query IDs of ``PUL_blast.out`` see the new family values.
#. ``substrate_prediction.tsv`` rows are sorted by contig, then by CGC number. Contig and protein IDs containing
   ``|`` are resolved against ``cgc_standard_out.tsv``. New optional tables ``substrate_prediction_dbcan_sub.tsv``
   (``-odbcan_sub true``) and ``substrate_prediction_dbcan_pul.tsv`` (``-odbcanpul true``).
#. GFF3 without gene rows (gffread style): transcripts are grouped into one gene by ``geneID``/``gene_id``/``gene``
   only when they share contig and strand and overlap, so repeated gene symbols (two ``tnpA`` copies, several
   ``trnL`` tRNAs) stay separate genes. Pseudogenes are non-coding rows.
#. ``dbcan_utils``: EC keys are single EC numbers; substrates are split on ``;`` and ``,``; all dbCAN-sub rows of a
   protein are used; ``NA`` placeholder rows no longer shift the ``CGC_abund`` columns; ``Fams`` shows the family of
   each member's Gene Type role.

.. rubric:: Figures

#. One palette for every CGC figure (the original dbCAN colours). A gene with several roles takes the colour of its
   highest-priority role (CAZyme > TC > TF > STP > Sulfatase > Peptidase); the other roles appear only in its label.
#. ``dbcan_plot CGC_plot``: scaled arrows, stacked short labels, kb ruler and scale bar, no flanking genes and no
   caption lines. ``CGC_coverage_plot`` uses the same gene track and accepts several ``--reads-count`` files.
#. Synteny figures (``substrate_prediction`` and ``dbcan_plot CGC_synteny_plot``) share one renderer: the dbCAN-PUL is
   aligned to the CGC and reverse-complemented when it lies on the opposite strand (``--synteny-revcomp``, default
   on); only hits that pass the substrate-prediction cut-offs are drawn, each labelled with its % identity. Ribbons are narrower than the genes and a gene with several hits gives
   each its own slot, so multiple matches stay visible. ``CGC_synteny_coverage_plot`` writes
   ``<cgc>_<PUL>-syntenic-cov.<ext>``.
#. ``cgc_circle_plot --layout auto|circos|overview``: circos for one contig, an overview of the contigs with CGCs
   otherwise; PNG and PDF instead of one SVG per contig. With ``--layout circos`` on several contigs each sector is
   named inside the innermost ring and the ticks are labelled on the largest sector only.
#. ``dbcan_plot``: ``-o`` is honoured by every plot; ``--vertical-bar`` draws vertical bars; ``seaborn-v0_8*`` styles
   are accepted.

.. rubric:: Fixes

#. ``--run_signalp`` without ``signalp6`` on ``PATH`` and ``--run_deeptmhmm`` without ``--deeptmhmm_dir`` stop with
   a usage error before the CAZyme searches start; ``--run_deeptmhmm`` works without ``--run_signalp``.
#. ``--threads`` reaches pyrodigal.
#. Protein FASTA in ``--mode protein``: a trailing ``*``/``.`` is removed and internal ones become ``X``.
#. Assembly-free DIAMOND counting: CAZy family tokens are parsed at any position, JGI sequence IDs keep the protein
   number, and ``.fq/.fastq/.fa/.fasta/.fna`` (optionally ``.gz``) read files are accepted.
#. The GFF match report counts a childless gene with a non-coding biotype (e.g. *E. coli* ``ralA``) as non-coding.
