.. _cgc_visualization:

CAZyme Gene Cluster Visualization
===================================

run_dbcan draws four kinds of CGC figures from the files of a run directory
(``cgc.gff``, ``cgc_standard_out.tsv`` and, when present, ``substrate_prediction.tsv``,
``PUL_blast.out`` and ``DEG.tsv``):

.. list-table::
   :widths: 28 32 40
   :header-rows: 1

   * - Figure
     - Command
     - Shows
   * - Genome view (circos or overview)
     - ``run_dbcan cgc_circle_plot``
     - Where the CGCs are on the genome / assembly and which have a predicted substrate
   * - CGC map
     - ``dbcan_plot CGC_plot``
     - The genes of one CGC
   * - CGC map with read depth
     - ``dbcan_plot CGC_coverage_plot``
     - The same map under one read-depth panel per sample
   * - Synteny with dbCAN-PUL
     - ``run_dbcan substrate_prediction`` (all CGCs), ``dbcan_plot CGC_synteny_plot`` /
       ``CGC_synteny_coverage_plot`` (one CGC)
     - One CGC above its best (or a chosen) dbCAN-PUL, homologous genes linked

Colours and labels (all figures)
---------------------------------

Every figure uses the same gene-type colours:

.. list-table::
   :widths: 30 30
   :header-rows: 1

   * - Gene type
     - Colour
   * - CAZyme
     - orange ``#E67E22``
   * - TC (transporter)
     - green ``#2ECC71``
   * - TF (transcription factor)
     - purple ``#9B59B6``
   * - STP (signal transduction)
     - yellow ``#F1C40F``
   * - Peptidase
     - teal ``#16A085``
   * - Sulfatase
     - dark blue ``#010E1B``
   * - Other (no signature hit, RNA)
     - grey ``#95A5A6``

A gene with several signature roles (for example ``CAZyme|GH13_31+TC|8.A.9.1.1`` in ``cgc_standard_out.tsv``) gets
**one** colour, chosen with the same priority as the ``Gene Type`` column of CGCFinder:
CAZyme > TC > TF > STP > Sulfatase > Peptidase. Its other roles appear only in its text label
(``GH13_31 + TC 8.A.9.1.1``). The ``SULFATLAS`` and ``PEPTIDASE`` spellings of the dbCAN-PUL reference files are
treated as Sulfatase and Peptidase. Legends list only the types present in the figure. When a dbCAN-PUL
reference gene carries both a TF number and a TF family name (``TF|0036286+TF|HTH_AraC``), the label shows the
family name (``HTH_AraC``).

The **CGC extents** of the genome views use two colours that are not gene-type colours, so that they cannot be
mistaken for genes: slate ``#34495E`` for a CGC with a substrate prediction and cool grey ``#ABB2B9`` for a CGC
without one.

Genome view: ``cgc_circle_plot``
---------------------------------

Run it after ``cgc_finder`` (and after ``substrate_prediction`` if you want substrate labels):

.. code-block:: shell

   run_dbcan cgc_circle_plot --output_dir output_EscheriaColiK12MG1655_fna
   run_dbcan cgc_circle_plot --output_dir output_Xylona_heveae_TC161_faa            # 23 scaffolds -> overview
   run_dbcan cgc_circle_plot --output_dir output_Xylona_heveae_TC161_faa --layout circos

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Parameter
     - Description
   * - ``--output_dir``
     - run_dbcan output directory; the figures are written to its ``cgc_circos/`` folder
   * - ``--layout auto|circos|overview``
     - ``auto`` (default): a circos plot when ``cgc.gff`` has exactly one contig (a complete circular genome), the
       overview otherwise. ``circos`` / ``overview`` force one layout.

**Circos layout** (``cgc_circos/cgc_circos_plot.png`` and ``.pdf``). Rings, from the outside in: coordinate axis;
CGC extents (slate = substrate predicted, cool grey = no prediction) labelled ``CGCn substrate``; the genes of the CGCs
in their gene-type colour; all CAZyme genes of the genome. A gene is placed on the CGC ring when its contig, start
and end match a row of ``cgc_standard_out.tsv``. With more than 60 CGCs only the CGCs with a predicted substrate are
labelled. When ``DEG.tsv`` (protein ID and log2 fold change, tab-separated, no header) exists in the output directory,
two inner rings show up-/down-regulated genes (red / blue) and the log2FC profile; genes without a numeric value
(``NA``) are not drawn as regulated. The figure is sized to the circle and the room its labels take, so a small genome with few
labels has no empty band around the rings. With ``--layout circos`` on several contigs, each contig is one sector
named inside the rings (sectors narrower than 2.5 degrees stay unnamed); every sector starts at 0, so the tick marks
are drawn on all sectors but labelled on the largest one only, and the subtitle gives the tick interval.

**Overview layout** (``cgc_circos/cgc_overview.png`` and ``.pdf``). One row per contig that carries a CGC (longest
first, up to 60 rows; all rows share one bp scale), the contig name on the left and its length at the end of the
bar. CGCs are blocks on the bar, slate with a ``CGCn substrate`` label when a substrate was predicted and cool
grey otherwise (when no substrate prediction exists the CGC ids are shown where they fit). Contigs without CGCs are not
drawn.

Contig lengths are taken from ``##sequence-region`` lines of ``cgc.gff`` when present, otherwise from the last
annotated feature of each contig.

CGC map: ``dbcan_plot CGC_plot``
---------------------------------

.. code-block:: shell

   dbcan_plot CGC_plot -i output_EscheriaColiK12MG1655_fna --cgcid 'NC_000913.3|CGC7' -o figures --format pdf,png

The title gives the CGC and its location (``CGC7  NC_000913.3:352,706–367,510``) and one line the predicted
substrate(s) (names only, from ``substrate_prediction.tsv``). Genes are arrows drawn to scale and pointing in the
direction of their strand, with a short family label above each gene. When labels would overlap they are stacked
in tiers and joined to their gene by a thin leader line; labels and leader lines never cross each other's text (a
lifted label may start or end at its leader instead of being centred on it).
A kb ruler and a scale bar are drawn below the genes and the gene-type legend at the bottom. Only the genes of the
CGC are drawn (no flanking genes). ``--show-protein-ids`` also prints the protein IDs under the genes.
Output: ``<contig>_<CGC>.cgc.<format>``.

CGC map with read depth: ``dbcan_plot CGC_coverage_plot``
-----------------------------------------------------------

.. code-block:: shell

   samtools depth -a -r 'NC_000913.3:352706-367510' ctrl1.bam > ctrl1.CGC7.depth.txt
   samtools depth -a -r 'NC_000913.3:352706-367510' treat1.bam > treat1.CGC7.depth.txt
   dbcan_plot CGC_coverage_plot -i output_dir --cgcid 'NC_000913.3|CGC7' \
       --reads-count ctrl1.CGC7.depth.txt --reads-count treat1.CGC7.depth.txt --samples ctrl1,treat1 -o figures

``--reads-count`` takes a ``samtools depth`` file (contig, position, depth; tab- or space-separated) and can be
repeated; each file becomes one panel above the gene map (first file on top). Only the rows of the CGC's contig
inside the CGC window are used, so a depth file of the whole assembly works too. The contig names must be the ones
of ``cgc.gff``: a depth file without any row for the CGC's contig (for example ``NC_000913`` instead of
``NC_000913.3``) is an error (exit status 1), not an empty panel. A contig that is present but has no row inside the
CGC window only gives a warning and a zero-depth panel. All panels share one y axis (``read depth``) so samples can be
compared; ``--samples`` names the panels (default: file names). Output: ``<contig>_<CGC>.cgc-coverage.<format>``.

Synteny with dbCAN-PUL
-----------------------

``run_dbcan substrate_prediction`` (and ``easy_substrate``) writes one PDF per CGC that has a dbCAN-PUL hit to
``synteny_pdf/<contig>_<CGC>_<PUL>-syntenic.pdf``. ``dbcan_plot`` draws the same figure for one CGC:

.. code-block:: shell

   dbcan_plot CGC_synteny_plot -i output_dir --cgcid 'NC_000913.3|CGC7' --db-dir db -o figures
   dbcan_plot CGC_synteny_plot -i output_dir --cgcid 'NC_000913.3|CGC7' --db-dir db --pul PUL0751 --no-synteny-revcomp
   dbcan_plot CGC_synteny_coverage_plot -i output_dir --cgcid 'NC_000913.3|CGC7' --db-dir db --reads-count ctrl1.CGC7.depth.txt

The CGC is drawn above the dbCAN-PUL (read from ``<db_dir>/dbCAN-PUL/<PUL>.out/cgc.gff``). The left labels name the
two tracks and the title gives ``CGC vs PUL`` and the PUL substrate. The reference is shifted so homologous genes sit
under each other and, with ``--synteny-revcomp`` (default on), it is drawn reverse-complemented when most homologous
gene pairs lie on opposite strands. Homologous genes are joined by ribbons whose grey level is the protein identity
(colour bar ``identity (%)``). Only DIAMOND hits from ``PUL_blast.out`` that pass the per-hit cut-offs of
substrate prediction are drawn: ``substrate_prediction`` uses its own ``-evalue``/``-iden``/``-cov`` values and
``dbcan_plot`` takes ``--evalue-cutoff`` (default 0.01), ``--identity-cutoff`` (percent, default 0) and
``--coverage-cutoff`` (fraction, default 0). ``--pul`` compares the CGC with another dbCAN-PUL than its best hit.
``CGC_synteny_coverage_plot`` adds the read-depth panels of ``CGC_coverage_plot`` above the CGC.

``dbcan_plot`` options
-----------------------

.. list-table::
   :widths: 32 68
   :header-rows: 1

   * - Option
     - Description
   * - ``-i`` / ``--input``
     - run_dbcan output directory (CGC figures) or comma-separated abundance files (``bar_plot``, ``heatmap_plot``)
   * - ``-o`` / ``--output-dir``
     - Directory for the figures (default: current directory; created if needed)
   * - ``--cgcid``
     - CGC id ``contig|CGCn`` as in ``substrate_prediction.tsv`` (contig names may contain ``|``)
   * - ``--format``
     - CGC figures: ``pdf`` (default), ``png``, ``svg`` or a comma-separated list; ``--dpi`` sets the PNG resolution
       (default 200)
   * - ``--db-dir``
     - Database directory holding ``dbCAN-PUL/`` (synteny figures)
   * - ``--pul``
     - dbCAN-PUL id for the synteny figures (default: best hit of the CGC)
   * - ``--reads-count``
     - samtools depth file; repeat for several samples
   * - ``--samples``
     - Comma-separated names: one per ``--reads-count`` file, or one per abundance file
   * - ``--synteny-revcomp`` / ``--no-synteny-revcomp``
     - Reverse-complement the dbCAN-PUL when it is on the opposite strand (default: on)
   * - ``--evalue-cutoff``, ``--identity-cutoff``, ``--coverage-cutoff``
     - Per-hit filters of the synteny links (same meaning and defaults as in ``substrate_prediction``)
   * - ``--show-protein-ids``
     - Print protein IDs under the genes of ``CGC_plot`` / ``CGC_coverage_plot``
   * - ``--top``, ``--col``, ``--value``
     - Abundance charts: number of rows, or the rows whose column ``--col`` has one of the ``--value`` values
   * - ``--vertical-bar``
     - ``bar_plot``: vertical bars (default: horizontal bars)
   * - ``--plot-style``
     - Matplotlib style of the abundance charts: ``ggplot`` (default), ``seaborn-v0_8``, ``seaborn-v0_8-poster``,
       ``seaborn-v0_8-whitegrid``, ``default``, ``classic``, ``bmh``, ``fivethirtyeight``; the old names ``seaborn``
       and ``seaborn-poster`` are mapped to the ``seaborn-v0_8`` ones
   * - ``--pdf``
     - File name of ``bar_plot`` (default ``bar_plot.pdf``) or ``heatmap_plot`` (default ``heatmap.pdf`` /
       ``heatmap_cluster.pdf``); a bare name is written into ``--output-dir``
   * - ``--show-abund``, ``--cluster-map``, ``--palette``, ``--show-fig``
     - ``heatmap_plot``: print values in the cells, cluster rows, colormap name. ``--show-fig`` also opens the
       heatmap in a window after writing the PDF; without a display (no interactive matplotlib backend) the PDF is
       still written but ``dbcan_plot`` exits with status 1

``dbcan_plot`` exits with status 1 and an error message when a figure cannot be drawn (unknown ``--cgcid``, missing
depth file, depth file without the CGC's contig, no dbCAN-PUL hit, empty abundance file, ...).

.. tip::
   PDF and SVG output keep the text editable in ``Inkscape`` or ``Adobe Illustrator``.
