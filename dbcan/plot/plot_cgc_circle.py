from pycirclize import Circos, config as circos_config  # Rename to avoid confusion
from pycirclize.parser import Gff
import pandas as pd
import os
from matplotlib.patches import Patch
import matplotlib.pyplot as plt

from dbcan.parameter import CGCPlotConfig
from dbcan.constants import (CGC_GFF_FILE, CGC_RESULT_FILE, CGC_CIRCOS_DIR,
                           CGC_CIRCOS_PLOT_FILE, CGC_CIRCOS_CONTIG_FILE_TEMPLATE,
                           CGC_FEATURE_TYPE, CGC_ANNOTATION_ATTR, PROTEIN_ID_ATTR,
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
                           DEG_LOG2FC_RANGE, DEG_TRACK_RANGE,CGC_FEATURE_LEGEND,DEG_FILE)
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class CGCCircosPlot:
    def __init__(self, config: CGCPlotConfig):
        self.config = config
        self.input_dir = config.output_dir.strip() if hasattr(config, 'output_dir') else ""
        self.gff_file = os.path.join(self.input_dir, CGC_GFF_FILE)
        self.tsv_file = os.path.join(self.input_dir, CGC_RESULT_FILE)
        self.output_dir = os.path.join(self.input_dir, CGC_CIRCOS_DIR)
        self.deg_tsv_file = os.path.join(self.input_dir, DEG_FILE)  # Adjust as needed

        # Validate file existence
        if not os.path.exists(self.gff_file):
            raise FileNotFoundError(f"GFF file not found: {self.gff_file}")
        if not os.path.exists(self.tsv_file):
            raise FileNotFoundError(f"TSV file not found: {self.tsv_file}")
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir, exist_ok=True)

        # Load GFF data
        self.gff = Gff(self.gff_file)
        self.seqid2size = self.gff.get_seqid2size()
        contig_count = len(self.seqid2size)
        # Set space based on contig count
        max_space = max(0, int(360 // contig_count) - 1)
        self.space = 0 if contig_count == 1 else min(2, max_space)
        self.circos = Circos(sectors=self.seqid2size, space=self.space)
        self.feature_type = CGC_FEATURE_TYPE
        self.seqid2features = self.gff.get_seqid2features(feature_type=self.feature_type)
        self.circos.text(CGC_PLOT_TITLE, size=CGC_TITLE_FONT_SIZE)

        # Load TSV data with enhanced error handling
        try:
            self.tsv_data = pd.read_csv(self.tsv_file, sep='\t')
            # Validate required columns
            required_columns = [CGC_ID_COLUMN, CONTIG_ID_COLUMN, CGC_PROTEIN_ID_FIELD,
                            GENE_START_COLUMN, GENE_STOP_COLUMN]
            missing_cols = [col for col in required_columns if col not in self.tsv_data.columns]
            if missing_cols:
                logging.warning(f"Missing required columns in TSV file: {missing_cols}")
        except Exception as e:
            logging.error(f"Error reading TSV file: {str(e)}")
            self.tsv_data = pd.DataFrame()  # Create empty DataFrame

        # Load DEG data
        deg_df = pd.read_csv(self.deg_tsv_file, sep="\t", header=None, names=["protein_id", "log2FC"])
        self.deg_data = deg_df[deg_df["log2FC"].notnull()]
        self.deg_data["log2FC"] = pd.to_numeric(self.deg_data["log2FC"], errors='coerce')



    def plot_feature_outer(self, circos=None):
        """Plot outer track with position markers"""
        if circos is None:
            circos = self.circos

        for sector in circos.sectors:
            outer_track = sector.add_track(CGC_OUTER_TRACK_RANGE)
            outer_track.axis(fc=CGC_AXIS_COLOR)
            major_interval = CGC_MAJOR_INTERVAL
            minor_interval = int(major_interval / CGC_MINOR_INTERVAL_DIVISOR)
            if sector.size > minor_interval:
                outer_track.xticks_by_interval(major_interval, label_formatter=lambda v: f"{v / 1000:.0f} Kb")
                outer_track.xticks_by_interval(minor_interval, tick_length=1, show_label=False)

    def plot_features_cazyme(self, circos=None, sector_name=None):
        """Plot CAZyme features"""
        if circos is None:
            circos = self.circos

        for sector in circos.sectors:
            if sector_name and sector.name != sector_name:
                continue

            cds_track = sector.add_track(CGC_CAZYME_TRACK_RANGE, r_pad_ratio=CGC_TRACK_PADDING)
            cds_track.axis(fc=CGC_TRACK_BG_COLOR, ec="none")
            cds_track.grid(2, color=CGC_GRID_COLOR)
            features = self.seqid2features[sector.name]
            for feature in features:
                if feature.type == self.feature_type:
                    cgc_type = feature.qualifiers.get(CGC_ANNOTATION_ATTR, ["unknown"])[0].split("|")[0]
                    if cgc_type == "CAZyme":  # only plot CAZyme features
                        color = self.get_feature_color(cgc_type)
                        cds_track.genomic_features(feature, fc=color)

    def plot_features_cgc(self, circos=None, sector_name=None):
        """Plot CGC features"""
        if circos is None:
            circos = self.circos

        for sector in circos.sectors:
            if sector_name and sector.name != sector_name:
                continue

            cds_track = sector.add_track(CGC_FEATURE_TRACK_RANGE, r_pad_ratio=CGC_TRACK_PADDING)
            cds_track.axis(fc=CGC_TRACK_BG_COLOR, ec="none")
            cds_track.grid(2, color=CGC_GRID_COLOR)
            features = self.seqid2features[sector.name]

            # Protect against empty DataFrame
            if not self.tsv_data.empty and CGC_PROTEIN_ID_FIELD in self.tsv_data.columns:
                cgc_ids_list = self.tsv_data[CGC_PROTEIN_ID_FIELD].unique().astype(str)
                for feature in features:
                    if feature.type == self.feature_type:
                        cgc_type = feature.qualifiers.get(CGC_ANNOTATION_ATTR, ["unknown"])[0].split("|")[0]
                        cgc_id  = str(feature.qualifiers.get(PROTEIN_ID_ATTR, ["unknown"])[0])
                        if cgc_id in cgc_ids_list:
                            color = self.get_feature_color(cgc_type)
                            cds_track.genomic_features(feature, fc=color)

    def plot_cgc_range(self, circos=None, sector_name=None):
        """Plot CGC range as rectangles"""
        if circos is None:
            circos = self.circos

        for sector in circos.sectors:
            if sector_name and sector.name != sector_name:
                continue

            cgc_track = sector.add_track(CGC_RANGE_TRACK_RANGE, r_pad_ratio=CGC_TRACK_PADDING)
            cgc_track.axis(fc=CGC_TRACK_BG_COLOR, ec="none")
            cgc_track.grid(2, color=CGC_GRID_COLOR)

            # Get sector size for validation
            sector_size = self.seqid2size[sector.name]

            # Filter data for current sector
            if self.tsv_data.empty or CONTIG_ID_COLUMN not in self.tsv_data.columns:
                continue

            # use sector name as string for comparison
            sector_data = self.tsv_data[self.tsv_data[CONTIG_ID_COLUMN].astype(str) == sector.name]

            # Process CGC ranges
            if CGC_ID_COLUMN in sector_data.columns:
                for cgc_id in sector_data[CGC_ID_COLUMN].unique():
                    cgc_rows = sector_data[sector_data[CGC_ID_COLUMN] == cgc_id]
                    if GENE_START_COLUMN in cgc_rows.columns and GENE_STOP_COLUMN in cgc_rows.columns:
                        try:
                            start = cgc_rows[GENE_START_COLUMN].min()
                            end = cgc_rows[GENE_STOP_COLUMN].max()

                            # verify coordinates are within sector size
                            if start >= sector_size or end > sector_size:
                                logging.warning(
                                    f"Skipping CGC {cgc_id} with coordinates ({start}-{end}) "
                                    f"that exceed sector '{sector.name}' size ({sector_size})"
                                )
                                continue

                            # make sure start < end
                            start = max(0, min(start, sector_size-1))
                            end = max(1, min(end, sector_size))

                            cgc_track.rect(start, end, fc=CGC_RANGE_COLOR, ec=CGC_RANGE_BORDER_COLOR)

                            # if end - start > sector_size * 0.01:
                            cgc_track.annotate((start + end) / 2, cgc_id, label_size=CGC_LABEL_SIZE)

                        except Exception as e:
                            logging.warning(f"Error plotting CGC {cgc_id} on {sector.name}: {str(e)}")

    def plot_log2fc_line(self, circos=None, sector_name=None):
        """add log2FC line for DEGs"""
        if circos is None:
            circos = self.circos
        if self.deg_data is None:
            logging.warning("No DEG data available for plotting")
            return

        for sector in circos.sectors:
            if sector_name and sector.name != sector_name:
                continue

            line_track = sector.add_track(DEG_LOG2FC_RANGE, r_pad_ratio=CGC_TRACK_PADDING)
            line_track.axis()
            line_track.grid(2, color=CGC_GRID_COLOR)

            features = self.seqid2features[sector.name]
            x = []
            y = []
            for feature in features:
                if feature.type == self.feature_type:
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
                # basic line
                line_track.line([min(x), max(x)], [20, 20], lw=1.5, ls="dotted", color="gray", vmin=vmin, vmax=vmax)
                # log2fc line
                line_track.line(x, y, color="pink", lw=1.5, vmin=vmin, vmax=vmax)

    def plot_deg_marker_circle(self, circos=None, sector_name=None):
        "add DEGs marker circle"
        if circos is None:
            circos = self.circos
        if self.deg_data is None:
            logging.warning("No DEG data available for plotting")
            return



        for sector in circos.sectors:
            if sector_name and sector.name != sector_name:
                continue

            deg_track = sector.add_track(DEG_TRACK_RANGE, r_pad_ratio=CGC_TRACK_PADDING)
            deg_track.axis(fc=CGC_TRACK_BG_COLOR, ec="none")
            deg_track.grid(2, color=CGC_GRID_COLOR)
            features = self.seqid2features[sector.name]

            for feature in features:
                if feature.type == self.feature_type:
                    protein_id = feature.qualifiers.get(PROTEIN_ID_ATTR, [""])[0]
                    if protein_id in self.deg_data['protein_id'].values:
                        log2fc = self.deg_data.loc[self.deg_data['protein_id'] == protein_id, 'log2FC'].iloc[0]
                        color = "#FF0000" if log2fc > 0 else "#4169E1"  # red for up, blue for down
                        deg_track.genomic_features(feature, fc=color, ec="black", lw=0.2)

    def get_feature_color(self, cgc_type):
        """Get color for different feature types"""
        return CGC_FEATURE_COLORS.get(cgc_type, "gray")

    def add_legend(self, circos=None):
        """Add legend to the plot"""
        if circos is None:
            circos = self.circos

        legend_colors = [self.get_feature_color(label) for label in CGC_FEATURE_LEGEND]
        rect_handles = []
        for idx, color in enumerate(legend_colors):
            rect_handles.append(Patch(color=color, label=CGC_FEATURE_LEGEND[idx]))
        # add DEG legend
        rect_handles.append(Patch(color="#FF0000", label="DEG up regulated"))
        rect_handles.append(Patch(color="#4169E1", label="DEG down regulated"))
        _ = circos.ax.legend(
            handles=rect_handles,
            bbox_to_anchor=CGC_LEGEND_POSITION,
            loc="center",
            fontsize=CGC_LEGEND_FONT_SIZE,
            title=CGC_LEGEND_TITLE,
            title_fontsize=CGC_LEGEND_FONT_SIZE,
            ncol=2,
        )

    def plot_single_contig(self, contig_name):
        """Plot a single contig and save to individual file"""
        try:
            # Create independent Circos object for this contig
            contig_size = {contig_name: self.seqid2size[contig_name]}
            contig_circos = Circos(sectors=contig_size, space=0)
            contig_circos.text(CGC_CONTIG_TITLE_TEMPLATE.format(contig_name=contig_name), size=CGC_TITLE_FONT_SIZE)

            # Add various features
            self.plot_feature_outer(contig_circos)
            self.plot_features_cazyme(contig_circos, contig_name)
            self.plot_features_cgc(contig_circos, contig_name)
            self.plot_cgc_range(contig_circos, contig_name)
            self.plot_log2fc_line(contig_circos, contig_name)
            self.plot_deg_marker_circle(contig_circos, contig_name)

            # Enable annotation adjustment to avoid overlap
            circos_config.ann_adjust.enable = True

            # Dynamically adjust figure size based on contig size
            size = min(CGC_MAX_FIGURE_SIZE, max(CGC_MIN_FIGURE_SIZE, CGC_MIN_FIGURE_SIZE + len(self.seqid2size) / CGC_FIGURE_SIZE_SCALING_FACTOR)) # Scale based on contig length
            fig = contig_circos.plotfig(figsize=(size, size))
            self.add_legend(contig_circos)

            # Save to file
            output_path = os.path.join(self.output_dir, CGC_CIRCOS_CONTIG_FILE_TEMPLATE.format(contig_name=contig_name))
            fig.savefig(output_path, format='svg', dpi=300)
            plt.close(fig)
            logging.info(f"Individual contig plot saved to: {output_path}")

        except Exception as e:
            logging.error(f"Error plotting contig {contig_name}: {str(e)}")

    def plot(self):
        """Plot everything - combined and individual contigs"""
        try:
            # 1. First plot containing all contigs
            self.plot_feature_outer()
            self.plot_features_cazyme()
            self.plot_features_cgc()
            self.plot_cgc_range()
            self.plot_log2fc_line()
            self.plot_deg_marker_circle()
            circos_config.ann_adjust.enable = True  # Avoid annotation overlap

            # Adjust figure size based on number of contigs
            size = min(CGC_MAX_FIGURE_SIZE, max(CGC_MIN_FIGURE_SIZE, CGC_MIN_FIGURE_SIZE + len(self.seqid2size) / CGC_FIGURE_SIZE_SCALING_FACTOR))
            fig = self.circos.plotfig(figsize=(size, size))
            self.add_legend()

            output_path = os.path.join(self.output_dir, CGC_CIRCOS_PLOT_FILE)
            fig.savefig(output_path, format='svg', dpi=300)
            plt.close(fig)  # Close the figure to free memory
            logging.info(f"Combined circos plot saved to: {output_path}")

            # 2. Then plot each contig individually
            total_contigs = len(self.seqid2size)
            logging.info(f"Creating individual plots for {total_contigs} contigs...")

            for idx, contig_name in enumerate(sorted(self.seqid2size.keys()), 1):
                logging.info(f"Processing contig {idx}/{total_contigs}: {contig_name}")
                self.plot_single_contig(contig_name)
                if idx % 10 == 0:
                    plt.close('all')  # Close all figures to free memory every 10 plots

                plt.close('all')  # Close all figures to free memory after each plot

        except Exception as e:
            logging.error(f"An error occurred: {str(e)}")
            import traceback
            logging.error(traceback.format_exc())  # Print full traceback for debugging

# if __name__ == "__main__":
#     config_dict = {
#         "gff_file": "/mnt/array2/xinpeng/dbcan_nf/dbCAN-xinpeng/bcb-unl-github/run_dbcan_new/dbcan_test/test_data/old_output/cgc.gff",
#         "tsv_file": "/mnt/array2/xinpeng/dbcan_nf/dbCAN-xinpeng/bcb-unl-github/run_dbcan_new/dbcan_test/test_data/old_output/cgc_standard_out.tsv",
#         "output_dir": "/mnt/array2/xinpeng/dbcan_nf/dbCAN-xinpeng/bcb-unl-github/run_dbcan_new/dbcan_test/test_data/old_output",
#     }
#     config = CGCPlotConfig.from_dict(CGCPlotConfig, config_dict)
#     plotter = CGCCircosPlot(config)
#     plotter.plot()


# if __name__ == "__main__":
#     config_dict = {
#         "gff_file": "/mnt/array2/xinpeng/dbcan_nf/dbCAN-xinpeng/bcb-unl-github/run_dbcan_new/dbcan_test/test_data/cgc.gff",
#         "tsv_file": "/mnt/array2/xinpeng/dbcan_nf/dbCAN-xinpeng/bcb-unl-github/run_dbcan_new/dbcan_test/test_data/cgc_standard_out.tsv",
#         "output_dir": "/mnt/array2/xinpeng/dbcan_nf/dbCAN-xinpeng/bcb-unl-github/run_dbcan_new/dbcan_test/test_data/",
#     }
#     config = CGCPlotConfig.from_dict(CGCPlotConfig, config_dict)
#     plotter = CGCCircosPlot(config)
#     plotter.plot()
