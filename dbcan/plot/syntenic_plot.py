"""CGC vs dbCAN-PUL synteny figures.

One renderer (:func:`plot_synteny`) is shared by ``run_dbcan substrate_prediction`` /
``easy_substrate`` (one PDF per CGC with a dbCAN-PUL hit, in ``synteny_pdf/``) and by
``dbcan_plot CGC_synteny_plot`` / ``CGC_synteny_coverage_plot``.

* The reference PUL is shifted so homologous genes line up with the CGC and, when
  ``revcomp`` is on (default) and most homologous pairs lie on opposite strands, it is
  drawn reverse-complemented.
* Links are the DIAMOND hits (PUL_blast.out) that pass the substrate-prediction per-hit
  cut-offs (E-value, identity, query coverage); their grey level is the % identity, which is also
  printed beside each link. Ribbons are narrower than the genes and a gene with several links gives
  each its own slot, so crossing links and multiple matches of one gene stay visible.
* Text is minimal: the two track names on the left, a title (CGC vs PUL) and the
  substrate. Gene colours use the shared dbCAN palette.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import matplotlib.colors as mcolors
from matplotlib import patheffects
import matplotlib.pyplot as plt
from matplotlib.path import Path as MPath
from matplotlib.patches import PathPatch

import dbcan.constants.cgc_substrate_prediction_constants as CSP
import dbcan.constants.plots_constants as P
from dbcan.IO.cgc_ids import cgc_query_id, parse_query_id, split_cgc_key
from dbcan.plot import cgc_viz as V
from dbcan.plot.cgc_map import draw_depth_panels, load_depth_tracks, substrate_text, PANEL_IN, PANEL_GAP_IN

logger = logging.getLogger(__name__)

U = P.SYNTENY_UNIT_IN
GENE_H = P.SYNTENY_GENE_H
Y_CGC = P.SYNTENY_CGC_TRACK_Y
TITLE_IN = P.SYNTENY_TITLE_IN
LEGEND_IN = P.SYNTENY_LEGEND_IN
R_IN = P.SYNTENY_RIGHT_IN
IDENTITY_CMAP = mcolors.LinearSegmentedColormap.from_list("dbcan_identity", list(P.SYNTENY_IDENTITY_RAMP))


@dataclass
class HomologyCutoffs:
    """Per-hit filters of substrate_prediction (same comparisons as its DIAMOND reader)."""
    evalue: float = CSP.DEFAULT_EVALUE_CUTOFF
    identity: float = CSP.DEFAULT_IDENTITY_CUTOFF
    coverage: float = CSP.DEFAULT_COVERAGE_CUTOFF

    @classmethod
    def from_config(cls, cfg) -> "HomologyCutoffs":
        def pick(name, default):
            v = getattr(cfg, name, None)
            return default if v is None else float(v)
        return cls(pick("evalue_cutoff", cls.evalue), pick("identity_cutoff", cls.identity),
                   pick("coverage_cutoff", cls.coverage))

    def passes(self, pident: float, evalue: float, qstart: int, qend: int, qlen: float) -> bool:
        if evalue > self.evalue or pident < self.identity:
            return False
        if qlen and (qend - qstart + 1) / qlen < self.coverage:
            return False
        return True


@dataclass
class Link:
    q: V.Gene
    s: V.Gene
    pident: float
    bitscore: float


def read_pul_links(blast_path: Path, cgcid: str, pul: str, cgc_genes: Sequence[V.Gene],
                   pul_genes: Sequence[V.Gene], cutoffs: HomologyCutoffs) -> List[Link]:
    """Best hit per (CGC gene, PUL gene) pair passing the cut-offs.

    Query IDs are read with the parser of substrate_prediction (dbcan.IO.cgc_ids): first the
    exact IDs that substrate_prediction writes to CGC.faa for the genes of this CGC, then the
    ``<contig>|<CGC#>|<protein>|<type>|<family>`` pattern; contigs and proteins may contain '|'.
    """
    blast_path = Path(blast_path)
    if not blast_path.exists():
        return []
    qby = {g.pid: g for g in cgc_genes}
    sby = {g.pid: g for g in pul_genes}
    cgc_parts = split_cgc_key(cgcid)
    known = {}
    if cgc_parts:
        for g in cgc_genes:
            q = cgc_query_id(g.contig, cgc_parts[1], g.pid, g.gene_type or "null", g.annotation)
            known[q.format()] = q
    best: Dict[Tuple[str, str], Link] = {}
    with blast_path.open() as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            query = parse_query_id(f[0], known)
            if query is None or query.cgc_id != cgcid:
                continue
            sp = f[1].split(":")
            if len(sp) < 3 or sp[1] != pul:
                continue
            q, s = qby.get(query.protein), sby.get(sp[2])
            if q is None or s is None:
                continue
            try:
                pident, evalue, bits = float(f[2]), float(f[10]), float(f[11])
                qstart, qend = int(f[6]), int(f[7])
                qlen = float(f[12]) if len(f) > 12 and f[12] else 0.0
            except ValueError:
                continue
            if not cutoffs.passes(pident, evalue, qstart, qend, qlen):
                continue
            key = (q.pid, s.pid)
            if key not in best or bits > best[key].bitscore:
                best[key] = Link(q, s, pident, bits)
    return list(best.values())


def load_pul_genes(db_dir: Path, pul: str) -> List[V.Gene]:
    gff = Path(db_dir) / P.DBCAN_PUL_DIR / f"{pul}.out" / P.DBCAN_PUL_GFF
    if not gff.exists():
        raise V.PlotInputError(f"dbCAN-PUL gff for {pul} not found: {gff}")
    genes = [g for gs in V.read_gff_genes(gff, feature_types=("CDS", "gene")).values() for g in gs]
    if not genes:
        raise V.PlotInputError(f"No genes in {gff}")
    return sorted(genes, key=lambda g: g.start)


def _ribbon(ax, qx0, qx1, sx0, sx1, ytop, ybot, color, zorder):
    ym = (ytop + ybot) / 2
    verts = [(qx0, ytop), (qx0, ym), (sx0, ym), (sx0, ybot), (sx1, ybot), (sx1, ym), (qx1, ym),
             (qx1, ytop), (qx0, ytop)]
    codes = [MPath.MOVETO, MPath.CURVE4, MPath.CURVE4, MPath.CURVE4, MPath.LINETO, MPath.CURVE4,
             MPath.CURVE4, MPath.CURVE4, MPath.CLOSEPOLY]
    ax.add_patch(PathPatch(MPath(verts, codes), facecolor=color, edgecolor=P.SYNTENY_LINK_EDGE_COLOR,
                           lw=P.SYNTENY_LINK_EDGE_LW, zorder=zorder))


def _ribbon_x_at(qx: float, sx: float, ytop: float, ybot: float, y: float) -> float:
    """x at height y of a ribbon edge running from qx (CGC track) to sx (PUL track)."""
    ym = (ytop + ybot) / 2
    best_t, best_d = 0.5, float("inf")
    for i in range(201):
        t = i / 200
        u = 1 - t
        yt = ytop * u ** 3 + 3 * ym * u * u * t + 3 * ym * u * t * t + ybot * t ** 3
        if abs(yt - y) < best_d:
            best_t, best_d = t, abs(yt - y)
    t, u = best_t, 1 - best_t
    return qx * u ** 3 + 3 * qx * u * u * t + 3 * sx * u * t * t + sx * t ** 3


def link_ports(links, pmap) -> Dict[int, Tuple[float, float, float, float]]:
    """Ends of each link ribbon: (qx0, qx1, sx0, sx1).

    A ribbon is narrower than its genes (P.SYNTENY_LINK_WIDTH_FRACTION of the gene width), so
    crossing ribbons never hide each other. A gene with several links (multiple matches) gets one
    slot per link, ordered by the partner's position, so the ribbons leave the gene side by side."""
    def slots(items, lo, hi):
        k = len(items)
        slot = (hi - lo) / k
        w = min(P.SYNTENY_LINK_WIDTH_FRACTION * (hi - lo), P.SYNTENY_LINK_SLOT_FILL * slot)
        return {id(l): (lo + (i + 0.5) * slot - w / 2, lo + (i + 0.5) * slot + w / 2) for i, l in enumerate(items)}

    smid = {id(l): sum(pmap(l.s)[:2]) / 2 for l in links}
    by_q: Dict[tuple, List[Link]] = {}
    by_s: Dict[tuple, List[Link]] = {}
    for l in links:
        by_q.setdefault((l.q.contig, l.q.start, l.q.end), []).append(l)
        by_s.setdefault((l.s.contig, l.s.start, l.s.end), []).append(l)
    q_end, s_end = {}, {}
    for ls in by_q.values():
        g = ls[0].q
        q_end.update(slots(sorted(ls, key=lambda l: smid[id(l)]), min(g.start, g.end), max(g.start, g.end)))
    for ls in by_s.values():
        a, b, _ = pmap(ls[0].s)
        s_end.update(slots(sorted(ls, key=lambda l: l.q.mid), min(a, b), max(a, b)))
    return {id(l): (*q_end[id(l)], *s_end[id(l)]) for l in links}


def draw_link_labels(ax, links, ports, ytop: float, ybot: float, x_per_in: float) -> int:
    """Write the % identity of every link beside its ribbon, half-way between the tracks (right
    side first, then left, then the same at the heights in P.SYNTENY_LINK_LABEL_LEVELS), so a
    number never covers a ribbon or another number. A thin white halo keeps it readable over
    crossing ribbons. Returns the number of labels written."""
    ym = (ytop + ybot) / 2
    fs = P.SYNTENY_LINK_LABEL_FONT_SIZE
    gap = P.SYNTENY_LINK_LABEL_GAP_IN * x_per_in
    placed = {lvl: [] for lvl in P.SYNTENY_LINK_LABEL_LEVELS}
    halo = [patheffects.withStroke(linewidth=2.2, foreground="white")]
    n = 0
    for l in sorted(links, key=lambda l: (l.q.mid, -l.pident)):
        qx0, qx1, sx0, sx1 = ports[id(l)]
        text = P.SYNTENY_LINK_LABEL_FORMAT.format(l.pident)
        w = V.text_width_data_in(text, fs) * x_per_in
        choice = None
        for lvl in P.SYNTENY_LINK_LABEL_LEVELS:
            y = ym + lvl
            right = max(_ribbon_x_at(qx1, sx1, ytop, ybot, y), _ribbon_x_at(qx0, sx0, ytop, ybot, y)) + gap
            left = min(_ribbon_x_at(qx1, sx1, ytop, ybot, y), _ribbon_x_at(qx0, sx0, ytop, ybot, y)) - gap - w
            for x0 in (right, left):
                if all(x0 + w <= p0 or x0 >= p1 for p0, p1 in placed[lvl]):
                    choice = (lvl, x0)
                    break
            if choice:
                break
        if choice is None:
            choice = (P.SYNTENY_LINK_LABEL_LEVELS[0], (qx1 + sx1) / 2 + gap)
        lvl, x0 = choice
        placed[lvl].append((x0, x0 + w))
        ax.text(x0, ym + lvl, text, ha="left", va="center", fontsize=fs, color=P.PLOT_INK, zorder=1.6,
                path_effects=halo)
        n += 1
    return n


def draw_links(ax, links, pmap, norm, x_per_in) -> List[Link]:
    """Draw every homology link as a narrow ribbon (grey by identity) and label it with its
    identity. Ribbons are narrower than the genes and a gene with several links gives each its own
    slot, so multiple matches (paralogues) and crossing links all stay visible. Returns the links."""
    ytop, ybot = Y_CGC - GENE_H / 2, GENE_H / 2
    ports = link_ports(links, pmap)
    for l in sorted(links, key=lambda l: (l.pident, l.bitscore)):
        _ribbon(ax, *ports[id(l)], ytop, ybot,
                IDENTITY_CMAP(norm(min(max(l.pident, P.SYNTENY_IDENTITY_MIN), P.SYNTENY_IDENTITY_MAX))),
                zorder=1.2)
    draw_link_labels(ax, links, ports, ytop, ybot, x_per_in)
    return list(links)


def reference_transform(cgc_genes, pul_genes, links, revcomp: bool):
    """Return (map(gene)->(x0, x1, strand), flipped) placing the PUL under the CGC."""
    opp = sum(l.bitscore for l in links if l.q.strand != l.s.strand)
    same = sum(l.bitscore for l in links if l.q.strand == l.s.strand)
    flip = bool(revcomp and opp > same)
    p_lo, p_hi = min(g.start for g in pul_genes), max(g.end for g in pul_genes)

    def oriented(g):
        if not flip:
            return g.start, g.end, g.strand
        return p_lo + p_hi - g.end, p_lo + p_hi - g.start, {"+": "-", "-": "+"}.get(g.strand, g.strand)

    diffs: List[float] = []
    for l in links:
        a, b, _ = oriented(l.s)
        diffs += [l.q.mid - (a + b) / 2] * max(1, int(l.bitscore // 100))
    diffs.sort()
    offset = diffs[len(diffs) // 2] if diffs else min(g.start for g in cgc_genes) - p_lo

    def pmap(g):
        a, b, st = oriented(g)
        return a + offset, b + offset, st

    return pmap, flip


def plot_synteny(cgc_genes: Sequence[V.Gene], cgcid: str, pul_genes: Sequence[V.Gene], pul: str,
                 links: Sequence[Link], out_dir: Path, *, substrates: Optional[Sequence[str]] = None,
                 revcomp: bool = True, depth_files=None, formats=("pdf",), dpi: int = P.PLOT_DEFAULT_DPI,
                 template: str = P.SYNTENY_FILE_TEMPLATE) -> List[Path]:
    if not cgc_genes or not pul_genes:
        raise V.PlotInputError(f"Nothing to draw for {cgcid} vs {pul}")
    V.setup_matplotlib()
    cgc_genes = sorted(cgc_genes, key=lambda g: g.start)
    pmap, _flip = reference_transform(cgc_genes, pul_genes, links, revcomp)
    placed = [pmap(g) for g in pul_genes]
    c_lo, c_hi = min(g.start for g in cgc_genes), max(g.end for g in cgc_genes)
    p_lo, p_hi = min(a for a, _, _ in placed), max(b for _, b, _ in placed)
    xmin, xmax = min(c_lo, p_lo), max(c_hi, p_hi)
    span = max(1, xmax - xmin)
    xlo, xhi = xmin - 0.012 * span, xmax + 0.012 * span
    head = 0.014 * span

    parsed = V.split_cgc_key(cgcid)
    cgc_name = parsed[1] if parsed else cgcid
    fs = P.PLOT_LABEL_FONT_SIZE
    name_w = max(V.text_width_data_in(t, P.PLOT_SUBTITLE_FONT_SIZE, fontweight="bold") for t in (cgc_name, pul))
    L_IN = 0.3 + name_w + 0.25
    W = max(9.0, min(17.0, L_IN + 3.5 + 0.55 * max(len(cgc_genes), len(pul_genes))))
    ax_w = W - L_IN - R_IN

    c_x = [g.mid for g in cgc_genes]
    p_x = [(a + b) / 2 for a, b, _ in placed]
    lay_c = V.layout_labels([g.label for g in cgc_genes], c_x, xlo, xhi, W, ax_w / W, fs,
                            extents=[(g.start, g.end) for g in cgc_genes])
    lay_p = V.layout_labels([g.label for g in pul_genes], p_x, xlo, xhi, W, ax_w / W, fs,
                            extents=[(a, b) for a, b, _ in placed])
    tier_h = (fs * 1.25 / 72) / U + 0.08
    y_top = Y_CGC + GENE_H / 2 + 0.3 + max(1, lay_c.ntiers) * tier_h + 0.1
    y_scale = -GENE_H / 2 - 0.3 - max(1, lay_p.ntiers) * tier_h - 0.75
    y_bot = y_scale - 0.75
    ax_h = (y_top - y_bot) * U

    tracks = []
    if depth_files:
        tracks = load_depth_tracks(depth_files, cgc_genes[0].contig, int(xlo), int(xhi) + 1)
    depth_in = len(tracks) * (PANEL_IN + PANEL_GAP_IN)
    H = TITLE_IN + depth_in + ax_h + LEGEND_IN
    fig = plt.figure(figsize=(W, H))
    ax = fig.add_axes([L_IN / W, LEGEND_IN / H, ax_w / W, ax_h / H])
    ax.set_xlim(xlo, xhi)
    ax.set_ylim(y_bot, y_top)
    ax.axis("off")

    ax.plot([c_lo, c_hi], [Y_CGC, Y_CGC], color=P.PLOT_BACKBONE_COLOR, lw=1.0, zorder=1)
    ax.plot([p_lo, p_hi], [0, 0], color=P.PLOT_BACKBONE_COLOR, lw=1.0, zorder=1)
    norm = mcolors.Normalize(P.SYNTENY_IDENTITY_MIN, P.SYNTENY_IDENTITY_MAX)
    draw_links(ax, links, pmap, norm, (xhi - xlo) / ax_w)
    for g in cgc_genes:
        V.draw_gene(ax, g.start, g.end, g.strand, g.gtype, Y_CGC, GENE_H, head)
    for g, (a, b, st) in zip(pul_genes, placed):
        V.draw_gene(ax, a, b, st, g.gtype, 0, GENE_H, head)
    V.draw_labels(ax, lay_c, Y_CGC + GENE_H / 2 + 0.3, tier_h, +1, fs, Y_CGC + GENE_H / 2 + 0.03, span)
    V.draw_labels(ax, lay_p, -GENE_H / 2 - 0.3, tier_h, -1, fs, -GENE_H / 2 - 0.03, span)

    # track names (left margin), vertically centred on each track
    for name, y in ((cgc_name, Y_CGC), (pul, 0.0)):
        fy = (LEGEND_IN + (y - y_bot) * U) / H
        fig.text(0.3 / W, fy, name, ha="left", va="center", fontsize=P.PLOT_SUBTITLE_FONT_SIZE,
                 fontweight="bold", color=P.PLOT_INK)

    V.draw_scale_bar(ax, xmax, y_scale, span)
    cb_w_in = 1.4
    cax = fig.add_axes([L_IN / W, (LEGEND_IN + (y_scale - y_bot) * U - 0.04) / H, cb_w_in / W, 0.08 / H])
    sm = plt.cm.ScalarMappable(cmap=IDENTITY_CMAP, norm=norm)
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal", ticks=[20, 60, 100])
    cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=P.PLOT_SMALL_FONT_SIZE - 0.5, length=2, colors=P.PLOT_INK_MUTED, pad=1)
    cax.text(1.04, 0.5, "identity (%)", transform=cax.transAxes, ha="left", va="center",
             fontsize=P.PLOT_SMALL_FONT_SIZE, color=P.PLOT_INK_MUTED)

    if tracks:
        draw_depth_panels(fig, tracks, cgc_genes, xlo, xhi, W, H, LEGEND_IN + ax_h + PANEL_GAP_IN,
                          left_in=L_IN, right_in=R_IN)

    fig.text(0.3 / W, 1 - 0.18 / H, f"{cgcid}  vs  {pul}", ha="left", va="top",
             fontsize=P.PLOT_TITLE_FONT_SIZE, fontweight="bold", color=P.PLOT_INK)
    fig.text(0.3 / W, 1 - 0.42 / H, substrate_text(substrates), ha="left", va="top",
             fontsize=P.PLOT_SUBTITLE_FONT_SIZE, color=P.PLOT_INK_MUTED)
    V.add_type_legend(fig, V.types_present(list(cgc_genes) + list(pul_genes)), y=0.02 / H)
    return V.save_figure(fig, out_dir, template, {"cgc": cgcid, "pul": pul}, formats, dpi)


def pul_substrate_from_db(db_dir: Path, pul: str) -> List[str]:
    """substrate_final of a PUL in dbCAN-PUL.xlsx (used when the PUL is not the CGC's best hit)."""
    xlsx = Path(db_dir) / CSP.PUL_EXCEL_FILE
    if not xlsx.exists():
        return []
    try:
        import pandas as pd
        table = pd.read_excel(xlsx)
        row = table[table["ID"].astype(str) == pul]
        if len(row):
            return V.SubstrateCall(pul, pul, str(row.iloc[0].get("substrate_final", "") or "")).names()
    except Exception as e:  # metadata is optional for the figure
        logger.warning(f"Could not read substrate of {pul} from {xlsx}: {e}")
    return []


def draw_cgc_vs_pul(run_dir: Path, db_dir: Path, cgcid: str, out_dir: Path, *, pul: Optional[str] = None,
                    cutoffs: Optional[HomologyCutoffs] = None, revcomp: bool = True, depth_files=None,
                    formats=("pdf",), dpi: int = P.PLOT_DEFAULT_DPI, template: str = P.SYNTENY_FILE_TEMPLATE,
                    cgc_table=None, substrate_calls=None) -> List[Path]:
    """Read a run directory and draw ``cgcid`` against ``pul`` (default: its best dbCAN-PUL hit)."""
    run_dir = Path(run_dir)
    cgcs = cgc_table if cgc_table is not None else V.read_cgc_table(run_dir / P.CGC_RESULT_FILE)
    if cgcid not in cgcs:
        raise V.PlotInputError(f"CGC '{cgcid}' not found in {run_dir / P.CGC_RESULT_FILE}")
    calls = substrate_calls if substrate_calls is not None else V.read_substrates(run_dir / P.CGC_SUB_PREDICTION_FILE)
    call = calls.get(cgcid)
    pul = pul or (call.pul if call else "")
    if not pul:
        raise V.PlotInputError(f"No dbCAN-PUL hit for {cgcid} in {run_dir / P.CGC_SUB_PREDICTION_FILE}; "
                               "run substrate_prediction first or pass --pul.")
    if call and call.pul == pul:
        substrates = V.SubstrateCall(cgcid, pul, call.pul_substrate).names()
    else:
        substrates = pul_substrate_from_db(db_dir, pul)
    pul_genes = load_pul_genes(db_dir, pul)
    links = read_pul_links(run_dir / P.PUL_DIAMOND_FILE, cgcid, pul, cgcs[cgcid], pul_genes,
                           cutoffs or HomologyCutoffs())
    logger.debug(f"Synteny {cgcid} vs {pul}: {len(links)} homology link(s) drawn")
    if not links:
        logger.warning(f"Synteny {cgcid} vs {pul}: no homology link passes the cut-offs in {P.PUL_DIAMOND_FILE}")
    return plot_synteny(cgcs[cgcid], cgcid, pul_genes, pul, links, out_dir, substrates=substrates,
                        revcomp=revcomp, depth_files=depth_files, formats=formats, dpi=dpi, template=template)


class SyntenicPlot:
    """``substrate_prediction`` step: one synteny PDF per CGC that has a dbCAN-PUL hit."""

    def __init__(self, config):
        self.config = config
        self.output_dir = Path(config.output_dir).resolve()
        self.db_dir = Path(config.db_dir).resolve()
        self.pdf_dir = self.output_dir / P.SYNTENY_PDF_DIR
        self.revcomp = bool(getattr(config, "synteny_revcomp", True))
        self.cutoffs = HomologyCutoffs.from_config(config)

    def syntenic_plot_allpairs(self) -> Tuple[int, int]:
        """Draw every CGC/PUL pair; returns (written, failed)."""
        sub_file = self.output_dir / P.CGC_SUB_PREDICTION_FILE
        self.pdf_dir.mkdir(parents=True, exist_ok=True)
        calls = V.read_substrates(sub_file)
        pairs = [c for c in calls.values() if c.pul]
        if not pairs:
            logger.info(f"No CGC with a dbCAN-PUL hit in {sub_file}; no synteny plots.")
            return 0, 0
        cgcs = V.read_cgc_table(self.output_dir / P.CGC_RESULT_FILE)
        ok = failed = 0
        for call in pairs:
            try:
                draw_cgc_vs_pul(self.output_dir, self.db_dir, call.cgcid, self.pdf_dir, pul=call.pul,
                                cutoffs=self.cutoffs, revcomp=self.revcomp, cgc_table=cgcs,
                                substrate_calls=calls)
                ok += 1
            except Exception as e:
                failed += 1
                logger.error(f"Synteny plot {call.cgcid} vs {call.pul} failed: {e}")
        logger.info(f"Generated {ok} synteny plot(s) in {self.pdf_dir}" + (f"; {failed} failed" if failed else ""))
        return ok, failed
