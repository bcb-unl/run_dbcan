Unreleased
~~~~~~~~~~

.. rubric:: Changed defaults and command lines

#. ``--mode`` accepts only ``prok``, ``meta`` or ``protein``. pyrodigal (``prok``/``meta``) is for bacteria, archaea
   and metagenomes; for fungi and other eukaryotes predict genes with funannotate, BRAKER or similar and run
   ``--mode protein`` with their GFF. A warning is logged when pyrodigal proteins look like a split eukaryotic genome.
#. CGC prediction from proteins (``easy_CGC``/``easy_substrate`` with ``--mode protein``) requires ``--input_gff``
   and ``--gff_type``, and ``gff_process`` requires ``--gff_type`` with any user ``--input_gff``;
   ``CAZyme_annotation`` still needs only the FASTA. New types: ``funannotate``, ``BRAKER``, ``GFF3``, ``GTF``
   (case-insensitive). The gene/protein match is checked first (before step 1 in ``easy_*``) and the run stops when
   fewer than ``--gff_min_match`` (default 0.9) of the protein-coding genes match; report: ``gff_resolution.tsv``.
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

#. ``cgc.gff``: on a CGC CAZyme (``#ofTools`` >= 2) the TC, TF (DIAMOND) and sulfatase labels are no longer used
   (full-length database hits cannot be placed relative to the CAZyme domain; CAZyme wins; the hits stay in
   ``total_cgc_info.tsv``). Domain-level labels (STP, fungal TF HMM, peptidase) are kept. *E. coli* K-12: 24 genes
   lose a label and CGC44 (``GT121 + TC 9.B.67`` on WaaL) is gone (53 -> 52 CGCs); *Xylona*: 27 genes, CGC111
   (``GT48 + TC 9.B.119`` on FKS1) is gone (111 -> 110). Substrate predictions are unchanged.
#. ``cgc.gff``: a non-CAZyme protein keeps only **one** full-length label (TC, TF DIAMOND or sulfatase): the hit
   with the highest DIAMOND bit score (ties: higher min(query, subject) coverage, lower E-value, TC > TF >
   Sulfatase; the E-value depends on database size). Domain-level labels are not affected; all hits stay in
   ``total_cgc_info.tsv``. The DIAMOND searches report ``BITSCORE`` (after ``Database`` in ``diamond.out.*``),
   ``total_cgc_info.tsv`` has a last column ``Bitscore`` (an older file without it keeps all labels), and each search
   keeps its best hit per protein by bit score. *E. coli* K-12: of 21 proteins with a TF and a TC label, 20
   regulators (two-component, LacI, NtrC families; TC hit 50-230 bits vs 147-1353 for TF) keep TF and their Gene
   Type changes from ``TC`` to ``TF``; RbsB (``NP_418207.1``) keeps ``TC|3.A.1.2.1`` (551 vs 54 bits). With the
   default ``--additional_genes TC`` the TF genes are no longer signatures: CGC7 (lac) loses ``lacI``, CGC9 five
   genes (bridged only by ``phoB``), CGC32 ``fhlA`` and one null gene. Still 52 CGCs with the same substrates; the
   PUL0751 score of CGC7 drops from 3968 to 3296. *Xylona* (``--fungi``), PUL0008 and PUL0230: unchanged.
#. Sulfatase hits must also cover >= 70% of the SulfAtlas protein (``--subject_coverage_threshold_sulfatase``):
   SulfAtlas entries are often sulfatases fused to GH domains, and a GH32 or GH2 protein that matched only the GH
   part was labelled ``S1_7`` / ``S1_19``. *E. coli* K-12: 113 -> 6 hits, *Xylona*: 267 -> 10; CGCs and substrates
   unchanged.
#. ``overview.tsv``: gene IDs are always read as text (numeric or ``NA`` IDs are no longer lost); the ``DIAMOND``
   column keeps only CAZy family tokens, once each; an EC without a count is written ``EC:1``; the ``Substrate``
   column looks the EC up across the eHMM's composition families, falls back to the family-level substrate, and uses
   ``dbCAN-sub-substrate.tsv`` from the database directory when present (same result as the dbCAN-sub web server).
   Only the tools selected with ``--methods`` are merged; a result table with rows left by an earlier run of an
   unselected tool is ignored with a warning (the empty placeholder table is ignored silently).
#. dbCAN-sub subfamily names are no longer rewritten by the ``GT2_*`` rename (it applies to ``dbCAN.hmm`` only).
#. CGC signatures: TC, TF (DIAMOND), Sulfatase and Peptidase are full-length DIAMOND searches and keep the best hit
   per protein without a coordinate overlap test against dbCAN domains. Peptidase hits are filtered by subject
   (MEROPS unit) coverage, and MEROPS inhibitor units are dropped. ``total_cgc_info.tsv`` has a ``Source`` column.
#. ``cgc_finder``: each gene counts once, by its final Gene Type (CAZyme > TC > TF > STP > Sulfatase > Peptidase);
   its other labels stay in ``Gene Annotation`` but fill no role. A ``CAZyme|GH13+STP|...`` gene next to a GH31 is
   no longer a CGC with ``--additional_genes STP``, and a ``TC|...+STP|...`` gene is only a TC. Defaults (``TC``):
   no change on *E. coli* K-12, *Xylona*, PUL0008 and PUL0230. With TC, TF, STP and Peptidase all required
   (``--additional_logic all``): *E. coli* 19 -> 14 CGCs, *Xylona* 23 -> 15; each lost CGC had its STP or peptidase
   only as a second label of a TC, TF or CAZyme gene. With ``--additional_logic any``: no change. A gene whose
   Gene Type is not listed is now also a null gene for joining, so CGCs can be split or trimmed: with
   ``--additional_genes STP`` alone *E. coli* goes from 46 to 28 CGCs. ``--additional_genes`` values are matched
   case-insensitively (``stp``, ``SULFATLAS``).
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

.. rubric:: Databases

#. New maintainer script ``python -m dbcan.utils.build_peptidase_db`` rebuilds ``peptidase_db.dmnd`` from MEROPS
   ``pepunit.lib``: identifiers are read from the bracket before the family field (the currently distributed
   database has 322 entries named after a bracket in the protein name or organism, e.g. ``misleading`` or
   ``Clostridium``), malformed identifiers and exact duplicates are left out, inhibitor units are kept as decoys
   (run_dbcan never reports them as Peptidase; ``--drop-inhibitors`` leaves them out), and ``diamond makedb`` is
   run. See :ref:`rebuild-peptidase-db`. The distributed database is not changed by this.
#. New ``fam-substrate-mapping.tsv`` (version line ``2026-11-10``, md5 ``0129c4b557ff3d5a9bc2ba316a66c2a7``,
   1,095 rows), built from the ``db_current`` table and checked against the CAZy family pages: 99 corrected rows
   (wrong class, EC or name), 64 activities listed by CAZy but missing from the table, rows for the families new in
   ``dbCAN.hmm`` V15 (CBM107-CBM115, CE25, CE26, GT139, GT140), family-level rows that keep a family's default class,
   one ``Substrate_curated`` style ("substrate (context)", ``;`` between substrates), and the new classes ``cutin``
   and ``mycothiol S-conjugate``. Over all 53,411 dbCAN-sub profiles, 744 change their substrate relative to the S3
   release ``db_v5-2-9_11-3-2026``; the release has to replace the table in both download locations.

.. rubric:: Figures

#. One palette for every CGC figure (the original dbCAN colours). A gene with several roles takes the colour of its
   highest-priority role (CAZyme > TC > TF > STP > Sulfatase > Peptidase); the other roles appear only in its label.
#. ``dbcan_plot CGC_plot``: scaled arrows, stacked short labels, kb ruler and scale bar, no flanking genes and no
   caption lines. ``CGC_coverage_plot`` uses the same gene track and accepts several ``--reads-count`` files.
#. Synteny figures (``substrate_prediction`` and ``dbcan_plot CGC_synteny_plot``) share one renderer: the dbCAN-PUL is
   aligned to the CGC and reverse-complemented when it lies on the opposite strand (``--synteny-revcomp``, default
   on); only hits that pass the substrate-prediction cut-offs are drawn, each labelled with its % identity in a tag
   on the ribbon (one row for the one-to-one pairs; tags of multiple matches above it, near the CGC end of their
   ribbon). Ribbons are narrower than the genes and a gene with several hits gives each its
   own slot, so multiple matches stay visible. ``CGC_synteny_coverage_plot`` writes
   ``<cgc>_<PUL>-syntenic-cov.<ext>``.
#. ``cgc_circle_plot --layout auto|circos|overview``: circos for one contig, an overview of the contigs with CGCs
   otherwise; PNG and PDF instead of one SVG per contig. With ``--layout circos`` on several contigs each sector is
   named inside the innermost ring and the ticks are labelled on the largest sector only.
#. ``dbcan_plot``: ``-o`` is honoured by every plot; ``--vertical-bar`` draws vertical bars; ``seaborn-v0_8*`` styles
   are accepted.
#. Documentation: new :ref:`figure_gallery` page with one example of every figure and the command that made it
   (real runs of *E. coli*, *Xylona heveae* and PUL0008); a pipeline overview figure (index, quick start) and a CGC
   rule figure (:ref:`cgc-rules-figure`); the user guide is rewritten to lead with the command and its outputs.

.. rubric:: Fixes

#. ``--run_signalp`` without ``signalp6`` on ``PATH`` and ``--run_deeptmhmm`` without ``--deeptmhmm_dir`` stop with
   a usage error before the CAZyme searches start; ``--run_deeptmhmm`` works without ``--run_signalp``.
#. ``--threads`` reaches pyrodigal.
#. Protein FASTA in ``--mode protein``: a trailing ``*``/``.`` is removed and internal ones become ``X``.
#. Assembly-free DIAMOND counting: CAZy family tokens are parsed at any position, JGI sequence IDs keep the protein
   number, and ``.fq/.fastq/.fa/.fasta/.fna`` (optionally ``.gz``) read files are accepted.
#. The GFF match report counts a childless gene with a non-coding biotype (e.g. *E. coli* ``ralA``) as non-coding.
#. ``dbcan_asmfree`` substrate abundance reads ``fam-substrate-mapping.tsv`` with leading ``#`` version lines and
   splits a class cell on ``;`` as well as ``,``.
#. SignalP6 or DeepTMHMM that exits 0 but writes no (or an empty) result file is a failed step (non-zero exit; the
   other tool's column is still written); an ``overview.tsv`` without genes gets the ``-`` columns of the requested tools.
#. CAZyme overview: a result table of a method selected with ``--methods`` that cannot be read, or that has rows but
   lacks required columns, fails the step (non-zero exit) instead of being skipped; empty tables are still accepted.

.. rubric:: Known issues / to test

#. ``dbCAN.hmm`` V15 merges the ``GT2_*`` models into one ``GT2`` profile (881 profiles). This is kept as is for
   now; its effect on GT2 sensitivity (fewer GT2 calls than with the split models were seen on test genomes) is to
   be tested before the next release. The fungal taxon profile on the dev branch is locked to the previous
   ``dbCAN.hmm`` and has to be recalibrated for V15.
#. ``fam-substrate-mapping.tsv``: for an incomplete EC (``a.b.c.-``) with several substrates under the same family
   (25 keys in the 2026-11-10 table, e.g. GH3 + 3.2.1.-) the last substrate alphabetically is used, as in 5.2.9. A
   per-subfamily substrate table (``dbCAN-sub-substrate.tsv``) in the database directory replaces this lookup. Two
   activities CAZy lists were left out of the table because of this rule (GH172 alpha-D-fructofuranosidase and
   PL38 exo-glucuronan lyase, both 3.2.1.-/4.2.2.-); their text is kept on the sibling row.
#. ``dbCAN-sub`` profiles of the old GT families (GT1-GT118) have no substrate rows (only GT119-GT140 do).
#. dbCAN-PUL has to be rebuilt with this version before the next database release: the PUL gene annotations
   (``dbCAN-PUL.xlsx`` / ``PUL.dmnd`` and the per-PUL GFF/CGC files) were made with 5.2.9 rules, and this version
   changes CGC membership (final Gene Type per gene, CAZyme wins over full-length labels, best full-length label by
   bitscore, sulfatase subject coverage), the TC family in ``PUL_blast.out`` query IDs (cut to three levels) and the
   substrate mapping. Re-run CAZyme annotation, ``gff_process`` and ``cgc_finder`` on every PUL, rebuild
   ``PUL.dmnd`` and the PUL table, and compare the PUL-to-substrate assignments with the current release.
#. ``dbcan_asmfree``: ``CAZyID_subfam_mapping.tsv`` and ``subfam_EC_mapping.tsv`` are not downloaded by
   ``run_dbcan database`` and are keyed by full CAZy FASTA headers, of which only about 56% match the headers of the
   current ``CAZy.dmnd``. Both tables have to be regenerated from the current CAZy data (TODO for the next database
   update).
