"""``run_dbcan cgc_circle_plot``: genome-wide view of the CGCs of one run.

``--layout auto`` (default) draws a circos plot when cgc.gff has exactly one contig (a
complete circular genome) and the overview ideogram (:mod:`dbcan.plot.cgc_overview`)
otherwise; ``--layout circos`` / ``--layout overview`` force one of them.

Circos rings, outside in: coordinate axis; CGC extents (labelled ``CGCn substrate``);
genes of the CGCs, coloured by Gene Type; all CAZyme genes; and, when ``DEG.tsv`` is
present, DEG markers and a log2FC line.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List, Optional

import matplotlib.pyplot as plt
from matplotlib.patches import Patch

import dbcan.constants.plot_cgc_circle_constants as C
import dbcan.constants.plots_constants as P
from dbcan.plot import cgc_viz as V
from dbcan.plot.cgc_overview import plot_overview

logger = logging.getLogger(__name__)


def read_deg(path: Path) -> Dict[str, float]:
    """DEG.tsv (protein_id <TAB> log2FC, no header) -> {protein_id: log2FC}; NA rows dropped."""
    out: Dict[str, float] = {}
    with Path(path).open() as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 2 or not f[0].strip():
                continue
            try:
                value = float(f[1])
            except ValueError:
                continue            # header, NA, empty
            if value == value:      # not NaN
                out[f[0].strip()] = value
    return out


class CGCCircosPlot:
    def __init__(self, config):
        self.config = config
        self.input_dir = Path(str(config.output_dir).strip())
        self.output_dir = self.input_dir / C.CGC_CIRCOS_DIR
        self.layout = (getattr(config, "layout", None) or C.CGC_DEFAULT_LAYOUT).lower()
        if self.layout not in C.CGC_LAYOUT_CHOICES:
            raise V.PlotInputError(f"Unknown layout '{self.layout}' (use {', '.join(C.CGC_LAYOUT_CHOICES)})")
        gff = self.input_dir / C.CGC_GFF_FILE
        tsv = self.input_dir / C.CGC_RESULT_FILE
        for f in (gff, tsv):
            if not f.exists():
                raise V.PlotInputError(f"Required input not found: {f}")
        self.genes_by_contig = V.read_gff_genes(gff)
        self.sizes = V.contig_sizes(gff)
        self.cgcs = V.read_cgc_table(tsv)
        self.calls = V.read_substrates(self.input_dir / C.CGC_SUB_PREDICTION_FILE)
        deg_file = self.input_dir / C.DEG_FILE
        self.deg: Optional[Dict[str, float]] = read_deg(deg_file) if deg_file.exists() else None
        if self.deg is not None and not self.deg:
            logger.warning(f"{deg_file} has no numeric log2FC values; DEG rings are not drawn.")
            self.deg = None

    def resolved_layout(self) -> str:
        if self.layout != "auto":
            return self.layout
        return "circos" if len(self.sizes) == 1 else "overview"

    def plot(self) -> List[Path]:
        if not self.sizes:
            raise V.PlotInputError(f"No features in {self.input_dir / C.CGC_GFF_FILE}")
        if not self.cgcs:
            logger.warning("No CGCs in cgc_standard_out.tsv; the figure shows no CGC.")
        layout = self.resolved_layout()
        logger.info(f"[CGC plot] layout={layout} ({len(self.sizes)} contig(s), {len(self.cgcs)} CGCs)")
        if layout == "overview":
            return plot_overview(self.cgcs, dict(self.sizes), self.calls, self.output_dir, formats=C.CGC_PLOT_FORMATS)
        return self.plot_circos()

    # ------------------------------------------------------------------ circos
    def _cgc_member_keys(self):
        """(contig, start, end) of every gene listed in cgc_standard_out.tsv."""
        return {(g.contig, g.start, g.end) for genes in self.cgcs.values() for g in genes}

    def cgc_ring_genes(self, contig: str) -> List[V.Gene]:
        """GFF genes of ``contig`` that belong to a CGC row (matched by contig, start, end).

        Matching by protein_id painted every RNA (protein_id=NA) onto the CGC ring as soon
        as one CGC contained an RNA.
        """
        members = self._cgc_member_keys()
        return [g for g in self.genes_by_contig.get(contig, []) if (g.contig, g.start, g.end) in members]

    def _cgc_label(self, cid: str) -> str:
        parts = V.split_cgc_key(cid)
        name = parts[1] if parts else cid
        call = self.calls.get(cid)
        names = call.names() if call else []
        return f"{name} {', '.join(names)}" if names else name

    def plot_circos(self) -> List[Path]:
        from pycirclize import Circos, config as circos_config

        V.setup_matplotlib()
        sizes = {c: max(1, s) for c, s in self.sizes.items()}
        space = 0 if len(sizes) == 1 else min(2, max(0, int(360 // len(sizes)) - 1))
        circos = Circos(sectors=sizes, space=space)
        cgc_genes = []
        total = sum(sizes.values())
        major = V.nice_step(total, 10 if len(sizes) == 1 else 24)    # one tick interval for all contigs
        label_all = len(self.cgcs) <= C.CGC_CIRCOS_MAX_ID_LABELS
        for sector in circos.sectors:
            genes = self.genes_by_contig.get(sector.name, [])
            size = sizes[sector.name]

            axis = sector.add_track(C.CGC_AXIS_TRACK_RANGE)
            axis.axis(fc="black", ec="none")
            if size > major:
                multi = len(sizes) > 1            # several contigs: no '0' at every sector start
                axis.xticks_by_interval(major, outer=False, tick_length=1.5, label_size=C.CGC_TICK_LABEL_SIZE,
                                        label_orientation="horizontal",
                                        label_formatter=lambda v, s=major, m=multi: (
                                            "" if (m and v == 0) else V.fmt_axis_tick(v, s, total)))

            rng = sector.add_track(C.CGC_RANGE_TRACK_RANGE)
            rng.axis(fc=C.CGC_TRACK_BG_COLOR, ec="none")
            feat = sector.add_track(C.CGC_FEATURE_TRACK_RANGE)
            feat.axis(fc=C.CGC_TRACK_BG_COLOR, ec="none")
            caz = sector.add_track(C.CGC_CAZYME_TRACK_RANGE)
            caz.axis(fc=C.CGC_TRACK_BG_COLOR, ec="none")

            for g in genes:
                lo, hi = max(0, g.start - 1), min(size, g.end)
                if hi <= lo:
                    continue
                if g.gtype == "CAZyme":
                    caz.rect(lo, hi, fc=V.gene_color("CAZyme"), ec="none")
            for g in self.cgc_ring_genes(sector.name):
                lo, hi = max(0, g.start - 1), min(size, g.end)
                if hi > lo:
                    feat.rect(lo, hi, fc=V.gene_color(g.gtype), ec="none")
                    cgc_genes.append(g)

            for cid, cg in self.cgcs.items():
                if cg[0].contig != sector.name:
                    continue
                lo, hi = min(g.start for g in cg) - 1, max(g.end for g in cg)
                if lo >= size:
                    logger.warning(f"{cid} ({lo}-{hi}) lies outside contig {sector.name} ({size} bp); skipped.")
                    continue
                hi = min(hi, size)
                predicted = bool(self.calls.get(cid) and self.calls[cid].names())
                rng.rect(lo, hi, fc=C.CGC_RANGE_PREDICTED_COLOR if predicted else C.CGC_RANGE_COLOR, ec="none")
                if not (predicted or label_all):
                    continue                      # many CGCs: only those with a substrate are labelled
                rng.annotate((lo + hi) / 2, self._cgc_label(cid), min_r=C.CGC_LABEL_MIN_R,
                             max_r=C.CGC_LABEL_MAX_R, label_size=C.CGC_LABEL_SIZE, shorten=None,
                             text_kws={"color": P.PLOT_INK if predicted else P.PLOT_INK_MUTED})

            if self.deg is not None:
                self._deg_tracks(sector, genes, size)

        circos_config.ann_adjust.enable = True
        W, H = C.CGC_FIGURE_SIZE, C.CGC_FIGURE_SIZE - 0.3
        fig = plt.figure(figsize=(W, H))
        side = 0.72                      # circle diameter as a fraction of the width
        ax = fig.add_axes([(1 - side) / 2, 0.095, side, side * W / H], projection="polar")
        circos.plotfig(ax=ax)

        name = next(iter(sizes)) if len(sizes) == 1 else f"{len(sizes)} contigs"
        fig.text(0.04, 0.985, f"CGCs of {name}", ha="left", va="top", fontsize=P.PLOT_TITLE_FONT_SIZE + 1,
                 fontweight="bold", color=P.PLOT_INK)
        n_pred = sum(1 for cid in self.cgcs if self.calls.get(cid) and self.calls[cid].names())
        fig.text(0.04, 0.958, f"{V.plural(len(self.cgcs), 'CGC')} · {n_pred} with substrate prediction · "
                 f"{V.fmt_bp(sum(sizes.values()))}", ha="left", va="top", fontsize=P.PLOT_SUBTITLE_FONT_SIZE,
                 color=P.PLOT_INK_MUTED)
        self._legend(fig, cgc_genes)
        return V.save_figure(fig, self.output_dir, C.CGC_CIRCOS_FILE_TEMPLATE, {}, C.CGC_PLOT_FORMATS,
                             P.PLOT_DEFAULT_DPI)

    def _deg_tracks(self, sector, genes, size):
        marks = sector.add_track(C.DEG_TRACK_RANGE)
        marks.axis(fc=C.CGC_TRACK_BG_COLOR, ec="none")
        line = sector.add_track(C.DEG_LOG2FC_RANGE)
        line.axis(fc="none", ec="#BDBDBD", lw=0.4)
        pts = []
        for g in genes:
            fc = self.deg.get(g.pid)
            if fc is None:
                continue          # no value (NA): neither up nor down
            lo, hi = max(0, g.start - 1), min(size, g.end)
            if fc != 0 and hi > lo:
                marks.rect(lo, hi, fc=C.DEG_UP_COLOR if fc > 0 else C.DEG_DOWN_COLOR, ec="none")
            pts.append((min(size, int(g.mid)), fc))
        if len(pts) > 1:
            pts.sort()
            vmax = max(1.0, max(abs(v) for _, v in pts))
            line.line([0, size], [0, 0], vmin=-vmax, vmax=vmax, color="#9A9A9A", lw=0.5, ls="dotted")
            line.line([p for p, _ in pts], [v for _, v in pts], vmin=-vmax, vmax=vmax, color=C.DEG_LINE_COLOR, lw=0.7)

    def _legend(self, fig, cgc_genes):
        hs, ls = V.type_legend_handles(V.types_present(cgc_genes))
        ls = [f"{t} gene" for t in ls]
        predicted = [bool(self.calls.get(cid) and self.calls[cid].names()) for cid in self.cgcs]
        if any(not p for p in predicted):
            hs, ls = [Patch(facecolor=C.CGC_RANGE_COLOR)] + hs, ["CGC"] + ls
        if any(predicted):
            hs, ls = [Patch(facecolor=C.CGC_RANGE_PREDICTED_COLOR)] + hs, ["CGC with substrate prediction"] + ls
        if self.deg is not None:
            hs += [Patch(facecolor=C.DEG_UP_COLOR), Patch(facecolor=C.DEG_DOWN_COLOR)]
            ls += ["up-regulated", "down-regulated"]
        fig.legend(hs, ls, loc="lower center", bbox_to_anchor=(0.5, 0.01), ncol=min(5, len(hs)), frameon=False,
                   fontsize=P.PLOT_LABEL_FONT_SIZE, handlelength=1.4, columnspacing=1.4)
