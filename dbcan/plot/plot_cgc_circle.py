from pathlib import Path
from pycirclize import Circos, config as circos_config
from pycirclize.parser import Gff
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as Patch
import csv
import logging

from dbcan.configs.base_config import CGCPlotConfig
from dbcan.constants.plot_cgc_circle_constants import (
    CGC_GFF_FILE, CGC_RESULT_FILE, CGC_CIRCOS_DIR,
    CGC_CIRCOS_PLOT_FILE, CGC_CIRCOS_CONTIG_FILE_TEMPLATE,
     CGC_ANNOTATION_ATTR, PROTEIN_ID_ATTR,
    CGC_ID_COLUMN, CONTIG_ID_COLUMN, CGC_PROTEIN_ID_FIELD,
    GENE_START_COLUMN, GENE_STOP_COLUMN,
    CGC_OUTER_TRACK_RANGE, CGC_CAZYME_TRACK_RANGE,
    CGC_FEATURE_TRACK_RANGE, CGC_RANGE_TRACK_RANGE,
    CGC_TRACK_PADDING, CGC_MAJOR_INTERVAL, CGC_MINOR_INTERVAL_DIVISOR,
    CGC_TRACK_BG_COLOR, CGC_GRID_COLOR, CGC_RANGE_COLOR,
    CGC_RANGE_BORDER_COLOR, CGC_AXIS_COLOR, CGC_LABEL_SIZE,
    CGC_LEGEND_POSITION, CGC_LEGEND_FONT_SIZE, CGC_TITLE_FONT_SIZE,
    CGC_FEATURE_COLORS, CGC_MIN_FIGURE_SIZE, CGC_MAX_FIGURE_SIZE,
    CGC_FIGURE_SIZE_SCALING_FACTOR, CGC_PLOT_TITLE,
    CGC_CONTIG_TITLE_TEMPLATE, CGC_LEGEND_TITLE,
    DEG_LOG2FC_RANGE, DEG_TRACK_RANGE, CGC_FEATURE_LEGEND, DEG_FILE
)

# Logging configuration is handled by main command or setup_logging()
# Removed hardcoded basicConfig to avoid interfering with global logging configuration

class CGCCircosPlot:
    def __init__(self, config: CGCPlotConfig):
        self.config = config

        if not self.gff_file.exists():
            raise FileNotFoundError(f"GFF file not found: {self.gff_file}")
        if not self.tsv_file.exists():
            raise FileNotFoundError(f"TSV file not found: {self.tsv_file}")
        if not self.output_dir.exists():
            self.output_dir.mkdir(parents=True, exist_ok=True)

        # Load GFF and sector sizes
        self.gff = Gff(str(self.gff_file))
        self.seqid2size = self.gff.get_seqid2size()
        contig_count = len(self.seqid2size)
        max_space = max(0, int(360 // contig_count) - 1)
        self.space = 0 if contig_count == 1 else min(2, max_space)
        self.circos = Circos(sectors=self.seqid2size, space=self.space)

        # Load all features first (do not pre-filter by a single type)
        self.seqid2features = self.gff.get_seqid2features(feature_type=["CDS", "gene"])
        detected_types = {f.type for feats in self.seqid2features.values() for f in feats}
        preferred_order = ["gene", "CDS"]
        self.allowed_feature_types = [t for t in preferred_order if t in detected_types]
        if not self.allowed_feature_types:
            logging.warning(f"[CGC Circle] No gene/CDS found. Falling back to all types: {detected_types}")
            self.allowed_feature_types = list(detected_types)
        else:
            logging.info(f"[CGC Circle] Using feature types: {self.allowed_feature_types} (detected: {detected_types})")
        for sid, feats in self.seqid2features.items():
            cnt_accept = sum(1 for f in feats if f.type in self.allowed_feature_types)
            logging.debug(f"[GFF] {sid} total={len(feats)} accepted={cnt_accept}")

        self.circos.text(CGC_PLOT_TITLE, size=CGC_TITLE_FONT_SIZE)

        # Load TSV
        try:
            self.tsv_data = pd.read_csv(self.tsv_file, sep='\t')
            required_columns = [
                CGC_ID_COLUMN, CONTIG_ID_COLUMN, CGC_PROTEIN_ID_FIELD,
                GENE_START_COLUMN, GENE_STOP_COLUMN
            ]
            missing_cols = [c for c in required_columns if c not in self.tsv_data.columns]
            if missing_cols:
                logging.warning(f"Missing required TSV columns: {missing_cols}")
        except Exception as e:
            logging.error(f"Error reading TSV file: {str(e)}")
            self.tsv_data = pd.DataFrame()

        # Load DEG (optional)
        if self.deg_tsv_file.exists():
            try:
                deg_df = pd.read_csv(self.deg_tsv_file, sep="\t", header=None, names=["protein_id", "log2FC"])
                self.deg_data = deg_df[deg_df["log2FC"].notnull()]
                self.deg_data["log2FC"] = pd.to_numeric(self.deg_data["log2FC"], errors='coerce')
            except Exception as e:
                logging.warning(f"Error reading DEG file {self.deg_tsv_file}: {e}")
                self.deg_data = None
        else:
            logging.info(f"DEG file not found: {self.deg_tsv_file}")
            self.deg_data = None

        self.load_substrate_labels()

    @property
    def input_dir(self) -> Path:
        return Path(self.config.output_dir.strip() if hasattr(self.config, 'output_dir') else "")

    @property
    def gff_file(self) -> Path:
        return self.input_dir / CGC_GFF_FILE

    @property
    def tsv_file(self) -> Path:
        return self.input_dir / CGC_RESULT_FILE

    @property
    def output_dir(self) -> Path:
        return self.input_dir / CGC_CIRCOS_DIR

    @property
    def deg_tsv_file(self) -> Path:
        return self.input_dir / DEG_FILE

    @property
    def substrate_file(self) -> Path:
        return self.input_dir / "substrate_prediction.tsv"

    def load_substrate_labels(self):
        self.cgcid2substrate = {}
        if not self.substrate_file.exists():
            logging.warning(f"Substrate prediction file not found: {self.substrate_file}")
            return
        try:
            with open(self.substrate_file, encoding="utf-8") as f:
                reader = csv.reader(f, delimiter="\t")
                header = next(reader, None)
                for row in reader:
                    if not row or row[0].startswith("#"):
                        continue
                    cgcid = row[0].strip() if len(row) > 0 else ""
                    pul_sub = row[2].strip() if len(row) > 2 else ""
                    dbsub_sub = row[5].strip() if len(row) > 5 else ""
                    labels = []
                    if pul_sub:
                        for sub in pul_sub.split(";"):
                            sub = sub.strip()
                            if sub:
                                labels.append(sub)
                    if dbsub_sub:
                        for sub in dbsub_sub.split(";"):
                            sub = sub.strip()
                            if sub:
                                labels.append(sub)
                    if cgcid and labels:
                        self.cgcid2substrate[cgcid] = "\n".join(labels)
        except Exception as e:
            logging.warning(f"Error reading substrate prediction file: {e}")

    def plot_feature_outer(self, circos=None):
        if circos is None:
            circos = self.circos
        for sector in circos.sectors:
            outer_track = sector.add_track(CGC_OUTER_TRACK_RANGE)
            outer_track.axis(fc=CGC_AXIS_COLOR)
            major_interval = CGC_MAJOR_INTERVAL
            minor_interval = int(major_interval / CGC_MINOR_INTERVAL_DIVISOR)
            if sector.size > minor_interval:
                outer_track.xticks_by_interval(
                    major_interval,
                    label_formatter=lambda v: f"{v / 1000:.0f} Kb"
                )
                outer_track.xticks_by_interval(
                    minor_interval,
                    tick_length=1,
                    show_label=False
                )

    def plot_features_cazyme(self, circos=None, sector_name=None):
        if circos is None:
            circos = self.circos
        for sector in circos.sectors:
            if sector_name and sector.name != sector_name:
                continue
            track = sector.add_track(CGC_CAZYME_TRACK_RANGE, r_pad_ratio=CGC_TRACK_PADDING)
            track.axis(fc=CGC_TRACK_BG_COLOR, ec="none")
            track.grid(2, color=CGC_GRID_COLOR)
            features = self.seqid2features[sector.name]
            for feature in features:
                if feature.type in self.allowed_feature_types:
                    cgc_type = feature.qualifiers.get(CGC_ANNOTATION_ATTR, ["unknown"])[0].split("|")[0]
                    if cgc_type == "CAZyme":
                        color = self.get_feature_color(cgc_type)
                        track.genomic_features(feature, fc=color)

    def plot_features_cgc(self, circos=None, sector_name=None):
        if circos is None:
            circos = self.circos
        for sector in circos.sectors:
            if sector_name and sector.name != sector_name:
                continue
            track = sector.add_track(CGC_FEATURE_TRACK_RANGE, r_pad_ratio=CGC_TRACK_PADDING)
            track.axis(fc=CGC_TRACK_BG_COLOR, ec="none")
            track.grid(2, color=CGC_GRID_COLOR)
            features = self.seqid2features[sector.name]
            if not self.tsv_data.empty and CGC_PROTEIN_ID_FIELD in self.tsv_data.columns:
                cgc_ids_list = self.tsv_data[CGC_PROTEIN_ID_FIELD].unique().astype(str)
                for feature in features:
                    if feature.type in self.allowed_feature_types:
                        cgc_type = feature.qualifiers.get(CGC_ANNOTATION_ATTR, ["unknown"])[0].split("|")[0]
                        cgc_id = str(feature.qualifiers.get(PROTEIN_ID_ATTR, ["unknown"])[0])
                        if cgc_id in cgc_ids_list:
                            color = self.get_feature_color(cgc_type)
                            track.genomic_features(feature, fc=color)

    def plot_cgc_range(self, circos=None, sector_name=None):
        if circos is None:
            circos = self.circos
        for sector in circos.sectors:
            if sector_name and sector.name != sector_name:
                continue
            cgc_track = sector.add_track(CGC_RANGE_TRACK_RANGE, r_pad_ratio=CGC_TRACK_PADDING)
            cgc_track.axis(fc=CGC_TRACK_BG_COLOR, ec="none")
            cgc_track.grid(2, color=CGC_GRID_COLOR)
            sector_size = self.seqid2size[sector.name]
            if self.tsv_data.empty or CONTIG_ID_COLUMN not in self.tsv_data.columns:
                continue
            sector_data = self.tsv_data[self.tsv_data[CONTIG_ID_COLUMN].astype(str) == sector.name]
            if CGC_ID_COLUMN in sector_data.columns:
                for cgc_id in sector_data[CGC_ID_COLUMN].unique():
                    rows = sector_data[sector_data[CGC_ID_COLUMN] == cgc_id]
                    if GENE_START_COLUMN in rows.columns and GENE_STOP_COLUMN in rows.columns:
                        try:
                            start = rows[GENE_START_COLUMN].min()
                            end = rows[GENE_STOP_COLUMN].max()
                            if start >= sector_size or end > sector_size:
                                logging.warning(
                                    f"Skipping CGC {cgc_id} ({start}-{end}) exceeds sector {sector.name} size {sector_size}"
                                )
                                continue
                            start = max(0, min(start, sector_size - 1))
                            end = max(1, min(end, sector_size))
                            cgc_track.rect(start, end, fc=CGC_RANGE_COLOR, ec=CGC_RANGE_BORDER_COLOR)
                            substrate_key = f"{sector.name}|{cgc_id}".strip()
                            substrate_label = self.cgcid2substrate.get(substrate_key, None)
                            reg_label = self.get_cgc_regulation_label(rows)
                            label = f"{cgc_id}"
                            if substrate_label:
                                label += f"\n{substrate_label}"
                            if reg_label:
                                label += f"\n{reg_label}"
                            cgc_track.annotate(
                                (start + end) / 2,
                                label,
                                label_size=CGC_LABEL_SIZE,
                                text_kws={"color": "red"} if substrate_label else {},
                                line_kws={"color": "red"} if reg_label == "up"
                                         else {"color": "blue"} if reg_label == "down"
                                         else {}
                            )
                        except Exception as e:
                            logging.warning(f"Error plotting CGC {cgc_id} on {sector.name}: {e}")

    def plot_log2fc_line(self, circos=None, sector_name=None):
        if circos is None:
            circos = self.circos
        if self.deg_data is None:
            logging.warning("No DEG data available")
            return
        for sector in circos.sectors:
            if sector_name and sector.name != sector_name:
                continue
            track = sector.add_track(DEG_LOG2FC_RANGE, r_pad_ratio=CGC_TRACK_PADDING)
            track.axis()
            track.grid(2, color=CGC_GRID_COLOR)
            features = self.seqid2features[sector.name]
            x, y = [], []
            for feature in features:
                if feature.type in self.allowed_feature_types:
                    protein_id = feature.qualifiers.get(PROTEIN_ID_ATTR, [""])[0]
                    pos = int((feature.location.start + feature.location.end) / 2)
                    if protein_id in self.deg_data['protein_id'].values:
                        log2fc = self.deg_data.loc[self.deg_data['protein_id'] == protein_id, 'log2FC'].iloc[0]
                        y_val = log2fc + 20
                    else:
                        y_val = 20
                    x.append(pos)
                    y.append(y_val)
            if len(x) > 1:
                x, y = zip(*sorted(zip(x, y)))
                vmin = min(y) - 1
                vmax = max(y) + 1
                track.line([min(x), max(x)], [20, 20], lw=1.5, ls="dotted", color="gray", vmin=vmin, vmax=vmax)
                track.line(x, y, color="pink", lw=1.5, vmin=vmin, vmax=vmax)

    def plot_deg_marker_circle(self, circos=None, sector_name=None):
        if circos is None:
            circos = self.circos
        if self.deg_data is None:
            logging.warning("No DEG data available")
            return
        for sector in circos.sectors:
            if sector_name and sector.name != sector_name:
                continue
            track = sector.add_track(DEG_TRACK_RANGE, r_pad_ratio=CGC_TRACK_PADDING)
            track.axis(fc=CGC_TRACK_BG_COLOR, ec="none")
            track.grid(2, color=CGC_GRID_COLOR)
            features = self.seqid2features[sector.name]
            for feature in features:
                if feature.type in self.allowed_feature_types:
                    protein_id = feature.qualifiers.get(PROTEIN_ID_ATTR, [""])[0]
                    if protein_id in self.deg_data['protein_id'].values:
                        log2fc = self.deg_data.loc[self.deg_data['protein_id'] == protein_id, 'log2FC'].iloc[0]
                        color = "#FF0000" if log2fc > 0 else "#4169E1"
                        track.genomic_features(feature, fc=color, ec="black", lw=0.2)

    def get_feature_color(self, cgc_type):
        return CGC_FEATURE_COLORS.get(cgc_type, "gray")

    def get_cgc_regulation_label(self, cgc_rows):
        if self.deg_data is None or cgc_rows.empty:
            return ""
        gene_ids = cgc_rows[CGC_PROTEIN_ID_FIELD].astype(str).tolist()
        deg_sub = self.deg_data[self.deg_data['protein_id'].astype(str).isin(gene_ids)]
        if deg_sub.empty:
            return ""
        up_count = (deg_sub['log2FC'] > 0).sum()
        down_count = (deg_sub['log2FC'] < 0).sum()
        total = len(gene_ids)
        if total == 0:
            return ""
        if up_count / total >= 0.5:
            return "up"
        if down_count / total >= 0.5:
            return "down"
        return ""

    def add_legend(self, circos=None):
        if circos is None:
            circos = self.circos
        legend_colors = [self.get_feature_color(label) for label in CGC_FEATURE_LEGEND]
        rect_handles = [Patch.Patch(color=color, label=CGC_FEATURE_LEGEND[idx])
                        for idx, color in enumerate(legend_colors)]
        rect_handles.append(Patch.Patch(color="#FF0000", label="DEG up regulated"))
        rect_handles.append(Patch.Patch(color="#4169E1", label="DEG down regulated"))
        circos.ax.legend(
            handles=rect_handles,
            bbox_to_anchor=CGC_LEGEND_POSITION,
            loc="center",
            fontsize=CGC_LEGEND_FONT_SIZE,
            title=CGC_LEGEND_TITLE,
            title_fontsize=CGC_LEGEND_FONT_SIZE,
            ncol=2,
        )

    def plot_single_contig(self, contig_name):
        try:
            contig_size = {contig_name: self.seqid2size[contig_name]}
            contig_circos = Circos(sectors=contig_size, space=0)
            contig_circos.text(
                CGC_CONTIG_TITLE_TEMPLATE.format(contig_name=contig_name),
                size=CGC_TITLE_FONT_SIZE
            )
            self.plot_feature_outer(contig_circos)
            self.plot_features_cazyme(contig_circos, contig_name)
            self.plot_features_cgc(contig_circos, contig_name)
            self.plot_cgc_range(contig_circos, contig_name)
            self.plot_log2fc_line(contig_circos, contig_name)
            self.plot_deg_marker_circle(contig_circos, contig_name)
            circos_config.ann_adjust.enable = True
            size = min(
                CGC_MAX_FIGURE_SIZE,
                max(CGC_MIN_FIGURE_SIZE,
                    CGC_MIN_FIGURE_SIZE + len(self.seqid2size) / CGC_FIGURE_SIZE_SCALING_FACTOR)
            )
            fig = contig_circos.plotfig(figsize=(size, size))
            self.add_legend(contig_circos)
            output_path = self.output_dir / CGC_CIRCOS_CONTIG_FILE_TEMPLATE.format(contig_name=contig_name)
            fig.savefig(output_path, format='svg', dpi=300)
            plt.close(fig)
            logging.info(f"Saved contig plot: {output_path}")
        except Exception as e:
            logging.error(f"Error plotting contig {contig_name}: {e}")

    def plot(self):
        try:
            self.plot_feature_outer()
            self.plot_features_cazyme()
            self.plot_features_cgc()
            self.plot_cgc_range()
            self.plot_log2fc_line()
            self.plot_deg_marker_circle()
            circos_config.ann_adjust.enable = True
            size = min(
                CGC_MAX_FIGURE_SIZE,
                max(CGC_MIN_FIGURE_SIZE,
                    CGC_MIN_FIGURE_SIZE + len(self.seqid2size) / CGC_FIGURE_SIZE_SCALING_FACTOR)
            )
            fig = self.circos.plotfig(figsize=(size, size))
            self.add_legend()
            output_path = self.output_dir / CGC_CIRCOS_PLOT_FILE
            fig.savefig(output_path, format='svg', dpi=300)
            plt.close(fig)
            logging.info(f"Saved combined plot: {output_path}")

            total_contigs = len(self.seqid2size)
            logging.info(f"Generating {total_contigs} individual contig plots...")
            for idx, contig_name in enumerate(sorted(self.seqid2size.keys()), 1):
                logging.info(f"[{idx}/{total_contigs}] {contig_name}")
                self.plot_single_contig(contig_name)
                if idx % 10 == 0:
                    plt.close('all')
                plt.close('all')
        except Exception as e:
            logging.error(f"Error during plotting: {e}")
            import traceback
            logging.error(traceback.format_exc())
