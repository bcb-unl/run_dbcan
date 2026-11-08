Unreleased
~~~~~~~~~~

.. rubric:: Changed outputs

#. ``CGC.faa`` headers (``<contig>|<CGC#>|<protein>|<Gene Type>|<family>``): the family is now that of the Gene
   Type role of a multi-role gene (CAZyme > TC > TF > STP > Sulfatase > Peptidase), not the second ``|`` field of
   the whole annotation. For example ``STP|LacI+TF|CCPA_PRIMG+STP|...`` with Gene Type ``TF`` gives ``CCPA_PRIMG``,
   and ``TF|LACI_ECOLI+STP|LacI+TC|3.A.1.2.1+...`` with Gene Type ``TC`` gives ``3.A.1``. On *E. coli* K-12 57
   headers change; ``substrate_prediction.tsv`` is identical after sorting (scoring uses the Gene Type only). Scripts
   that parse the query IDs of ``PUL_blast.out`` see the new family values.
#. ``substrate_prediction.tsv`` rows are sorted by contig, then by CGC number.
#. GFF3 without gene rows (gffread style): transcripts are grouped into one gene by ``geneID``/``gene_id``/``gene``
   only when they share contig and strand and overlap, so repeated gene symbols (two ``tnpA`` copies, several
   ``trnL`` tRNAs) stay separate genes.
#. ``cgc_circle_plot --layout circos`` on several contigs names each sector inside the innermost ring and labels
   the ticks on the largest sector only.

.. rubric:: Fixes

#. ``--run_signalp`` without ``signalp6`` on ``PATH`` and ``--run_deeptmhmm`` without ``--deeptmhmm_dir`` stop with
   a usage error before the CAZyme searches start.
