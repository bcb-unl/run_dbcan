from pathlib import Path
import logging
import warnings
import matplotlib.pyplot as plt
from matplotlib.patches import PathPatch, Polygon
from matplotlib.lines import Line2D
from typing import Dict, Tuple, List
import matplotlib.path as mpath

from dbcan.configs.cgc_substrate_config import SynPlotConfig

import dbcan.constants.plots_constants as plots_constants

warnings.filterwarnings("ignore", category=UserWarning)
logger = logging.getLogger(__name__)

class SyntenicPlot:
    """Syntenic plots between CGCs and PULs"""

    def __init__(self, config: SynPlotConfig):
        self.config = config
        self.output_dir = Path(config.output_dir).resolve()
        self.db_dir = Path(config.db_dir).resolve()

        self.input_sub_out = self.output_dir / plots_constants.CGC_SUB_PREDICTION_FILE
        self.blastp = self.output_dir / plots_constants.PUL_DIAMOND_FILE
        self.cgc = self.output_dir / plots_constants.CGC_RESULT_FILE

        fallback_sub = self.output_dir / "substrate_prediction.tsv"
        if (not self.input_sub_out.exists()) and fallback_sub.exists():
            logging.warning(f"[substrate] {self.input_sub_out.name} not found, using fallback file: {fallback_sub.name}")
            self.input_sub_out = fallback_sub

        self.pdf_dir = self.output_dir / "synteny_pdf"
        self.pdf_dir.mkdir(parents=True, exist_ok=True)

        logging.info(f"[synteny init] output_dir={self.output_dir}")
        logging.info(f"[synteny init] substrate_file={self.input_sub_out} (exists={self.input_sub_out.exists()})")
        logging.info(f"[synteny init] blast_file={self.blastp} (exists={self.blastp.exists()})")
        logging.info(f"[synteny init] cgc_file={self.cgc} (exists={self.cgc.exists()})")
        logging.info(f"[synteny init] dbCAN-PUL dir={self.db_dir / 'dbCAN-PUL'} (exists={(self.db_dir / 'dbCAN-PUL').is_dir()})")

    def syntenic_plot_allpairs(self):
        """overall function to plot all CGC-PUL pairs"""
        if not self.blastp.exists() or self.blastp.stat().st_size == 0:
            logger.warning(f"[skip] BLAST result not found or empty: {self.blastp}")
            return
        cgcpul_blastp = read_blast_result_cgc(str(self.blastp))
        logging.info(f"[load] BLAST pair groups: {len(cgcpul_blastp)}")

        if not self.cgc.exists():
            logger.error(f"[abort] CGC standard out not found: {self.cgc}")
            return
        cgc_proteinid2gene, cgcid2gene, cgcid2geneid = read_UHGG_CGC_stanrdard_out(str(self.cgc))
        logging.info(f"[load] CGC IDs: {len(cgcid2gene)} (proteins={len(cgc_proteinid2gene)})")

        PULid_proteinid2gene, PULid2gene, PULid2geneid = self.read_PUL_cgcgff()
        logging.info(f"[load] PUL IDs: {len(PULid2gene)} (proteins={len(PULid_proteinid2gene)})")

        if not self.input_sub_out.exists() or self.input_sub_out.stat().st_size == 0:
            logger.warning(f"[skip] Substrate prediction file not found or empty: {self.input_sub_out}")
            return

        plot_count = 0
        candidate_pairs = 0
        no_blast_pairs = 0
        missing_cgc = 0
        missing_pul = 0

        with self.input_sub_out.open() as fh:
            header = next(fh, "")
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 2:
                    continue
                cgc, pul = parts[0].strip(), parts[1].strip()
                substrate = parts[2].strip() if len(parts) > 2 else ""
                if not pul:
                    continue
                candidate_pairs += 1

                if cgc not in cgcid2gene:
                    missing_cgc += 1
                    logger.debug(f"[miss CGC] {cgc} not in cgc_standard_out")
                    continue
                if pul not in PULid2gene:
                    missing_pul += 1
                    logger.debug(f"[miss PUL] {pul} not in PUL gff db")
                    continue

                cgcpul_key = f"{cgc}:{pul}"

                bed_cgc = cgcid2gene[cgc]
                bed_pul = PULid2gene[pul]

                starts1, ends1, strands1, types1 = Get_parameters_for_plot(bed_cgc)
                starts2, ends2, strands2, types2 = Get_parameters_for_plot(bed_pul)
                genes1 = cgcid2geneid[cgc]
                genes2 = PULid2geneid[pul]

                blocks = []
                if cgcpul_key in cgcpul_blastp:
                    for rec in cgcpul_blastp[cgcpul_key]:
                        qseqid = rec.qseqid
                        sseqid = rec.sseqid
                        try:
                            cgc_proteinid = qseqid.split("|")[2]
                        except Exception:
                            continue
                        _, pul_proteinid = parse_pul_ids(sseqid)
                        try:
                            idx1 = genes1.index(cgc_proteinid)
                            idx2 = genes2.index(pul_proteinid)
                            blocks.append(f"{idx1}-{idx2}-{rec.pident}")
                        except ValueError:
                            continue
                else:
                    no_blast_pairs += 1
                    logger.debug(f"[no BLAST] {cgcpul_key} no alignment blocks; plotting empty synteny")

                self._create_syntenic_plot(
                    starts1, starts2, ends1, ends2,
                    strands1, strands2, types1, types2,
                    blocks, cgc, pul, substrate
                )
                plot_count += 1

        if plot_count == 0:
            logger.warning(
                "[diagnostic]error candidate=%d, blast_groups=%d, CGC=%d, PUL=%d. "
                "missing_cgc=%d, missing_pul=%d, no_blast=%d",
                candidate_pairs, len(cgcpul_blastp), len(cgcid2gene), len(PULid2gene),
                missing_cgc, missing_pul, no_blast_pairs
            )
        else:
            logger.info(
                "Generated %d syntenic plots (candidate=%d, missing_cgc=%d, missing_pul=%d, no_blast=%d)",
                plot_count, candidate_pairs, missing_cgc, missing_pul, no_blast_pairs
            )

    def read_PUL_cgcgff(self):
        """
        Read PUL cgc.gff files to get gene annotations.
        """
        PULidgeneid2gene: Dict[str, CGC_stanrdard] = {}
        pul_dir = self.db_dir / "dbCAN-PUL"
        if not pul_dir.is_dir():
            logger.warning(f"dbCAN-PUL directory not found: {pul_dir}")
            return {}, {}, {}

        for entry in pul_dir.iterdir():
            if entry.is_dir() and entry.name.startswith("PUL") and entry.name.endswith(".out"):
                gff_path = entry / "cgc.gff"
                read_cgcgff(gff_path, PULidgeneid2gene)

        cgcid2gene: Dict[str, List[CGC_stanrdard]] = {}
        cgcid2geneid: Dict[str, List[str]] = {}
        for PULidgeneid, gene in PULidgeneid2gene.items():
            cgcid2gene.setdefault(gene.CGCID, []).append(gene)
            cgcid2geneid.setdefault(gene.CGCID, []).append(gene.Protein_ID)

        return PULidgeneid2gene, cgcid2gene, cgcid2geneid

    def _create_syntenic_plot(self, starts, starts1, ends, ends1, strands, strands1, types, types1, blocks, cgcid, pulid, substrate):
        plot_config = {'output_dir': self.output_dir}
        syntenic_plot(starts, starts1, ends, ends1, strands, strands1,
                      types, types1, blocks, cgcid, pulid, substrate, plot_config)


class blastp_hit(object):
    def __init__(self, lines):
        # Expecting outfmt 6: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore [qlen slen]
        self.qseqid = lines[0]
        self.sseqid = lines[1]
        self.pident = float(lines[2]) if len(lines) > 2 else 0.0
        self.length = int(lines[3]) if len(lines) > 3 else 0
        self.mismatch = int(lines[4]) if len(lines) > 4 else 0
        self.gapopen = int(lines[5]) if len(lines) > 5 else 0
        self.qstart = int(lines[6]) if len(lines) > 6 else 0
        self.qend = int(lines[7]) if len(lines) > 7 else 0
        self.sstart = int(lines[8]) if len(lines) > 8 else 0
        self.send = int(lines[9]) if len(lines) > 9 else 0
        self.evalue = float(lines[10]) if len(lines) > 10 else 1.0
        self.bitscore = float(lines[11]) if len(lines) > 11 else 0.0
        self.qlen = int(lines[12]) if len(lines) > 12 else 0
        self.slen = int(lines[13]) if len(lines) > 13 else 0

    def __repr__(self):
        return "\t".join([str(self.__dict__[attr]) for attr in self.__dict__])


def identity_map(seqsim: float) -> str:
    if 80 <= seqsim <= 100:
        return "red"
    if 60 <= seqsim < 80:
        return "blue"
    if 40 <= seqsim < 60:
        return "green"
    if 20 <= seqsim < 40:
        return "cyan"
    return "gray"


def Get_Position(starts, ends, strands, maxbp, yshift=0, up=1):
    Width = 1000
    Height = 160
    poly_heigth = 5
    Triangle_length = 4

    # Copy to avoid modifying the caller's list
    starts = list(starts)
    ends = list(ends)
    strands = list(strands)

    plot_start_y = Height / 2 - poly_heigth - yshift
    polygens = []
    # Shift positions
    shfit_pos = min(starts) if starts else 0
    for i in range(len(starts)):
        starts[i] -= shfit_pos
        ends[i] -= shfit_pos

    pixeachbp = Width / maxbp if maxbp > 0 else 1.0

    blocks = []
    lines = []
    for i in range(len(starts)):
        if strands[i] == "+":
            positions_str = f"{starts[i] * pixeachbp} {plot_start_y} "
            positions_str += f"{ends[i] * pixeachbp - Triangle_length} {plot_start_y} "
            if up == 1:
                blocks.append(positions_str)
            positions_str += f"{ends[i] * pixeachbp} {plot_start_y + poly_heigth} "
            positions_str += f"{ends[i] * pixeachbp - Triangle_length} {plot_start_y + 2*poly_heigth} "
            positions_str += f"{starts[i] * pixeachbp} {plot_start_y + 2*poly_heigth}"

            positions_str1 = f"{starts[i] * pixeachbp} {plot_start_y + 2*poly_heigth} "
            positions_str1 += f"{ends[i] * pixeachbp - Triangle_length} {plot_start_y + 2*poly_heigth} "
            if up == 2:
                blocks.append(positions_str1)
        else:
            positions_str = f"{starts[i] * pixeachbp} {plot_start_y + poly_heigth} "
            positions_str += f"{starts[i] * pixeachbp + Triangle_length} {plot_start_y} "
            positions_str += f"{ends[i] * pixeachbp} {plot_start_y} "
            positions_str1 = f"{ends[i] * pixeachbp} {plot_start_y} "
            positions_str1 += f"{starts[i] * pixeachbp + Triangle_length} {plot_start_y} "
            if up == 1:
                blocks.append(positions_str1)
            positions_str += f"{ends[i] * pixeachbp} {plot_start_y + 2* poly_heigth} "
            positions_str += f"{starts[i]* pixeachbp + Triangle_length} {plot_start_y + 2* poly_heigth}"
            positions_str1 = f"{ends[i] * pixeachbp} {plot_start_y + 2* poly_heigth} "
            positions_str1 += f"{starts[i]* pixeachbp + Triangle_length} {plot_start_y + 2* poly_heigth}"
            if up == 2:
                blocks.append(positions_str1)

        polygens.append(positions_str)
        if i < len(starts) - 1:
            positions_str = f"{ends[i] * pixeachbp} {plot_start_y + poly_heigth}  "
            positions_str += f"{starts[i+1]*pixeachbp} {plot_start_y + poly_heigth}"
            lines.append(positions_str)

    return polygens, blocks, lines, [], []


def plot_Polygon_homologous(polygens1, polygens2, types1, types2, size, ax):
    # colors_map = {
    #     "CAZYME": "#E67E22",
    #     "TC": "#2ECC71",
    #     "TF": "#9B59B6",
    #     "STP": "#F1C40F",
    #     "PEPTIDASE": "#16A085",
    #     "SULFATASE": "#010E1B",
    #     "OTHER": "#95A5A6"
    # }
    colors_map = plots_constants.CGC_FEATURE_COLORS
    default_color = plots_constants.CGC_FEATURE_COLORS.get("Other")
    
    for j in range(len(polygens1)):
        polygen = polygens1[j].split()
        pts = [[float(polygen[2*i]), float(polygen[2*i+1])] for i in range(len(polygen)//2)]
        t = (types1[j] or "Other")
        color = colors_map.get(t, default_color)
        ax.add_patch(Polygon(pts, color=color, alpha=0.5, lw=0))

    for j in range(len(polygens2)):
        polygen = polygens2[j].split()
        pts = [[float(polygen[2*i]), float(polygen[2*i+1])] for i in range(len(polygen)//2)]
        t = (types2[j] or "Other")
        color = colors_map.get(t, default_color)
        ax.add_patch(Polygon(pts, color=color, alpha=0.5, lw=0))


def decode_block(block: str) -> Tuple[int, int, float]:
    o1, o2, sim = block.split("-")
    return int(o1), int(o2), float(sim)


def points2(coord: str) -> Tuple[float, float, float, float]:
    x1, y1, x2, y2 = coord.split()
    return float(x1), float(y1), float(x2), float(y2)


def Shade_curve(x11,x12,y11,y12,x21,x22,y21,y22,xmid,ymid,color):
    # Use matplotlib.path.Path (avoid confusion with pathlib.Path)
    PathCls = mpath.Path
    M, C4, L, CP = PathCls.MOVETO, PathCls.CURVE4, PathCls.LINETO, PathCls.CLOSEPOLY
    pathdata = [
        (M,  (x11, y11)),
        (C4, (x11, ymid)),
        (C4, (x21, ymid)),
        (C4, (x21, y21)),
        (L,  (x22, y22)),
        (C4, (x22, ymid)),
        (C4, (x12, ymid)),
        (C4, (x12, y12)),
        (CP, (x11, y11))
    ]
    codes = [c for c, _ in pathdata]
    verts = [p for _, p in pathdata]
    path = PathCls(verts, codes)
    return PathPatch(path, color=color, alpha=0.2, lw=0)


def plot_syntenic_block(blocks, blocks1_coor, blocks2_coor, ax):
    for block in blocks:
        order1, order2, sim = decode_block(block)
        coord1 = blocks1_coor[order1]
        coord2 = blocks2_coor[order2]
        x11, y11, x12, y12 = points2(coord1)
        x21, y21, x22, y22 = points2(coord2)

        color = identity_map(sim)
        xmid = (x11 + x22 + x21 + x22) / 4
        ymid = (y11 + y22) / 2
        ax.add_patch(Shade_curve(x11, x12, y11, y12, x21, x22, y21, y22, xmid, ymid, color))


def plot_genome_line(lines_coor1, lines_coor2, ax):
    for line in lines_coor1:
        x1, y1, x2, y2 = points2(line)
        ax.add_patch(Polygon([(x1, y1), (x2, y2)], color="gray", lw=2))
    for line in lines_coor2:
        x1, y1, x2, y2 = points2(line)
        ax.add_patch(Polygon([(x1, y1), (x2, y2)], color="gray", lw=2))


def syntenic_plot(starts, starts1, ends, ends1, strands, strands1, types, types1,
                  blocks, cgcid, pulid, substrate, config):
    custom_lines = [
        Line2D([0], [0], color="red", lw=4, alpha=0.5),
        Line2D([0], [0], color="blue", lw=4, alpha=0.5),
        Line2D([0], [0], color="green", lw=4, alpha=0.5),
        Line2D([0], [0], color="cyan", lw=4, alpha=0.5),
        Line2D([0], [0], color="gray", lw=4, alpha=0.5)
    ]
    labels = ["80-100", "60-80", "40-60", "20-40", "0-20"]

    # genelabelcolor = ["#E67E22", "#2ECC71", "#9B59B6", "#F1C40F", "#16A085", "#34495E", "#95A5A6"]
    # geneslabels = ["CAZyme", "TC", "TF", "STP", "PEPTIDASE", "SULFATASE", "Other"]
    genelabelcolor= plots_constants.GENE_LABEL_COLOR
    geneslabels = plots_constants.GENE_LABELS
    genecustom_lines = [Line2D([0], [0], color=c, lw=6, alpha=0.5) for c in genelabelcolor]

    px = 1/plt.rcParams['figure.dpi']
    width = 1600
    height = 320*2
    fig = plt.figure(figsize=((width*px), (height*px*2/4)))
    ax = fig.add_subplot(111)

    maxbp = max([max(ends) - min(starts), max(ends1) - min(starts1)]) if starts and starts1 else 1

    polygens, blocks_coor, lines_coor, _, _ = Get_Position(starts, ends, strands, maxbp, yshift=0, up=1)
    polygens1, blocks1_coor, lines_coor1, _, _ = Get_Position(starts1, ends1, strands1, maxbp, yshift=40, up=2)

    plot_Polygon_homologous(polygens, polygens1, types, types1, 2, ax)
    plot_syntenic_block(blocks, blocks_coor, blocks1_coor, ax)
    plot_genome_line(lines_coor, lines_coor1, ax)

    legend1 = ax.legend(custom_lines, labels, frameon=False, loc='upper left', title="Identity")
    ax.add_artist(legend1)
    legend2 = ax.legend(genecustom_lines, geneslabels, frameon=False, loc='lower left', title="Gene")
    ax.add_artist(legend2)

    ax.text(500, 90, cgcid, fontsize=18, ha='center')
    ax.text(500, 0, pulid, fontsize=18, ha='center')

    # Substrate annotation (top-right)
    if substrate:
        ax.text(0.98, 0.98, f"Substrate: {substrate}", transform=ax.transAxes,
                ha='right', va='top', fontsize=12)

    ax.set_ylim(0, 100)
    ax.set_xlim(-100, 1100)
    ax.axis('off')
    plt.tight_layout(pad=0.01)

    pdf_dir = Path(config['output_dir']) / "synteny_pdf"
    pdf_dir.mkdir(parents=True, exist_ok=True)
    cgcid_safe = cgcid.replace("|", "_")
    pulid_safe = pulid.replace("|", "_")
    out_path = pdf_dir / f"{cgcid_safe}_{pulid_safe}-syntenic.pdf"
    plt.savefig(out_path)
    plt.close()
    logger.debug(f"Saved synteny plot: {out_path}")


def read_blast_result_cgc(filename: str) -> Dict[str, List[blastp_hit]]:
    """Read BLAST results and group by 'CGC|contig:PUlid'."""
    querydict: Dict[str, List[blastp_hit]] = {}
    if not Path(filename).exists() or Path(filename).stat().st_size == 0:
        return querydict
    with open(filename) as fh:
        for line in fh:
            parts = line.strip().split()
            if len(parts) < 12:
                continue
            qseqid = parts[0]
            sseqid = parts[1]
            qids = qseqid.split("|")
            if len(qids) < 2:
                continue
            queryid = f"{qids[0]}|{qids[1]}"
            pulid, _ = parse_pul_ids(sseqid)
            key = f"{queryid}:{pulid}" if pulid else queryid
            querydict.setdefault(key, []).append(blastp_hit(parts))
    return querydict


def parse_pul_ids(sseqid: str) -> Tuple[str, str]:
    """ Parse PUL ID and protein ID from sseqid."""
    parts = sseqid.split(":")
    pulid = parts[1] if len(parts) > 1 else ""
    prot = parts[2] if len(parts) > 2 else (parts[-1] if parts else "")
    return pulid, prot


def Get_parameters_for_plot(CGC_stanrdard_list):
    starts = []
    ends = []
    strands = []
    types = []
    for gene in CGC_stanrdard_list:
        starts.append(gene.Gene_Start)
        ends.append(gene.Gene_END)
        strands.append(gene.Strand)
        types.append(gene.Gene_Type)
    return starts, ends, strands, types


def read_cgcgff(filename, geneid2gene):
    path = Path(filename)
    if not path.exists():
        return None
    # Derive PUL ID from parent directory name (e.g. PUL0648.out)
    pul_dir_name = path.parent.name  # e.g. PUL0648.out
    pul_id_from_dir = pul_dir_name.split('.')[0] if pul_dir_name.startswith("PUL") else pul_dir_name
    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 9:
                continue
            desc = fields[-1]
            protein_id = None
            feature = "Other"
            ann = "Other"
            for item in desc.split(";"):
                if item.startswith(plots_constants.GFF_PROTEIN_ID_COL):
                    protein_id = item.split("=", 1)[1]
                elif item.startswith(plots_constants.GFF_CGC_ANNOTATION_COL):
                    ann = item.split("=", 1)[1]
                    feature = (ann.split("+")[0].split("|")[0] if "|" in ann else ann.split("+")[0]) or "Other"
            if not protein_id:
                continue
            PULid = pul_id_from_dir
            newline = [
                PULid,
                feature,
                fields[0],
                protein_id,
                fields[3],
                fields[4],
                fields[6],
                ann
            ]
            geneid2gene[f"{PULid}:{protein_id}"] = CGC_stanrdard(newline)


class CGC_stanrdard(object):
    def __init__(self, lines):
        self.CGCID = lines[0]
        self.Gene_Type = lines[1]
        self.Contig_ID = lines[2]
        self.Protein_ID = lines[3]
        self.Gene_Start = int(lines[4])
        self.Gene_END = int(lines[5])
        self.Strand = lines[6]
        self.Protein_Family = lines[7]

    def __repr__(self):
        return "\t".join([str(self.__dict__[attr]) for attr in self.__dict__])


def read_UHGG_CGC_stanrdard_out(filename):
    geneid2gene = {}
    cgcid2gene = {}
    cgcid2geneid = {}
    with open(filename) as fh:
        for line in fh:
            if line.startswith(plots_constants.CGC_ID_FIELD) or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            fields[0] = fields[2] + "|" + fields[0]
            gene = CGC_stanrdard(fields)
            geneid2gene[gene.Protein_ID] = gene
            cgcid2gene.setdefault(gene.CGCID, []).append(gene)
            cgcid2geneid.setdefault(gene.CGCID, []).append(gene.Protein_ID)
    return geneid2gene, cgcid2gene, cgcid2geneid
