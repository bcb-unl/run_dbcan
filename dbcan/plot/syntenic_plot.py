import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch
import matplotlib.patches as patches
from matplotlib.patches import FancyArrowPatch
from matplotlib.patches import Polygon
from matplotlib.lines import Line2D
from matplotlib import pyplot
from matplotlib.patches import Patch
import matplotlib
import warnings
import  os

# Example blocks: ["0-0-30","2-10-50"]
warnings.filterwarnings("ignore", category=UserWarning)

class SyntenicPlot:
    """Class for generating syntenic plots between CGCs and PULs"""
    
    def __init__(self, config):
        """
        Initialize syntenic plot generator with configuration
        
        Args:
            config: Configuration object with required paths and settings
        """
        # Store configuration reference
        self.config = config
        
        # Extract parameters from configuration
        self.output_dir = config.output_dir
        self.input_sub_out = getattr(config, 'input_sub_out', 
                                    os.path.join(self.output_dir, "substrate_prediction.tsv"))
        self.blastp = getattr(config, 'blastp', 
                            os.path.join(self.output_dir, "PUL_blast.out"))
        self.cgc = getattr(config, 'cgc', 
                         os.path.join(self.output_dir, "cgc_standard_out.tsv"))
        self.db_dir = getattr(config, 'db_dir', "db")
        
        # Create output directory for PDF files
        self.pdf_dir = os.path.join(self.output_dir, "synteny_pdf")
        os.makedirs(self.pdf_dir, exist_ok=True)
        
    def syntenic_plot_allpairs(self, config):
        """
        Generate syntenic plots for all CGC-PUL pairs
        
        Args:
            config: Configuration object (usually same as in __init__)
        """
        # Read BLAST results
        cgcpul_blastp = read_blast_result_cgc(self.blastp)
        
        # Read CGC and PUL gene data
        cgc_proteinid2gene, cgcid2gene, cgcid2geneid = read_UHGG_CGC_stanrdard_out(self.cgc)
        PULid_proteinid2gene, PULid2gene, PULid2geneid = self.read_PUL_cgcgff()
        
        # Create output directory for PDF files
        pdf_directory = os.path.join(self.output_dir, "synteny_pdf")
        os.makedirs(pdf_directory, exist_ok=True)
        
        # Process each CGC-PUL pair in input file
        plot_count = 0
        with open(self.input_sub_out) as file:
            for line in file.readlines()[1:]:  # Skip header line
                lines = line.rstrip().split("\t")
                if len(lines) < 2:
                    continue
                    
                cgc = lines[0]
                pul = lines[1]
                
                # Skip if no PUL match found
                if not pul:
                    continue
                # Skip if CGC or PUL data not found
                if cgc not in cgcid2gene or pul not in PULid2gene:
                    continue
                    
                cgcpul = cgc + ":" + pul
                
                # Skip if no BLAST results for this pair
                if cgcpul not in cgcpul_blastp:
                    continue
                    
                bed_cgc = cgcid2gene[cgc]
                bed_pul = PULid2gene[pul]

                # Extract plot parameters
                starts1, ends1, strands1, types1 = Get_parameters_for_plot(bed_cgc)
                starts2, ends2, strands2, types2 = Get_parameters_for_plot(bed_pul)
                genes1 = cgcid2geneid[cgc]
                genes2 = PULid2geneid[pul]
                
                # Build blocks for homologous regions
                blocks = []
                for record in cgcpul_blastp[cgcpul]:
                    query = record.qseqid
                    hit = record.sseqid
                    cgc_proteinid = query.split("|")[2]
                    
                    # Extract PUL protein ID
                    split_parts = hit.split(":")
                    pul_proteinid = split_parts[2] if len(split_parts) > 3 else split_parts[3]
                    #print(f"Processing {cgc} - {pul}: {cgc_proteinid} - {pul_proteinid}")
                    
                    try:
                        index1 = genes1.index(cgc_proteinid)
                        index2 = genes2.index(pul_proteinid)
                        blocks.append(f"{index1}-{index2}-{record.pident}")
                    except ValueError:
                        continue
                
                # Generate plot
                self._create_syntenic_plot(
                    starts1, starts2, ends1, ends2, 
                    strands1, strands2, types1, types2, 
                    blocks, cgc, pul
                )
                plot_count += 1
                
        print(f"Generated {plot_count} syntenic plots in {pdf_directory}")
    
    def read_PUL_cgcgff(self):
        """
        Read PUL genome annotation information
        
        Returns:
            Tuple (PULidgeneid2gene, cgcid2gene, cgcid2geneid)
        """
        PULidgeneid2gene = {}
        
        # Ensure database path is correctly formatted
        dbCANPUL_data_folder = self.db_dir if self.db_dir.endswith("/") else self.db_dir + "/"
        if not dbCANPUL_data_folder.startswith("/"):
            dbCANPUL_data_folder = os.path.abspath(dbCANPUL_data_folder) + "/dbCAN-PUL/"
        else:
            dbCANPUL_data_folder += "/dbCAN-PUL/"
        
        # Check if directory exists
        if not os.path.exists(dbCANPUL_data_folder):
            print(f"Warning: Could not find dbCAN-PUL data directory: {dbCANPUL_data_folder}")
            return {}, {}, {}
        
        # Process each PUL file
        folder = os.listdir(dbCANPUL_data_folder)
        for filename in folder:
            if filename.startswith("PUL") and filename.endswith(".out"):
                read_cgcgff(dbCANPUL_data_folder + filename + "/cgc.gff", PULidgeneid2gene)
        
        # Build lookup dictionaries
        cgcid2gene = {}
        cgcid2geneid = {}
        for PULidgeneid in PULidgeneid2gene:
            gene = PULidgeneid2gene[PULidgeneid]
            cgcid2gene.setdefault(gene.CGCID, []).append(gene)
            cgcid2geneid.setdefault(gene.CGCID, []).append(gene.Protein_ID)
        
        return PULidgeneid2gene, cgcid2gene, cgcid2geneid
        
    def _create_syntenic_plot(self, starts, starts1, ends, ends1, strands, strands1, types, types1, blocks, cgcid, pulid):
        """
        Generate a single syntenic plot
        
        Args:
            starts, starts1: Start positions of genes in CGC and PUL
            ends, ends1: End positions of genes in CGC and PUL
            strands, strands1: Strand directions of genes in CGC and PUL
            types, types1: Types of genes in CGC and PUL
            blocks: Homologous block information
            cgcid: CGC ID
            pulid: PUL ID
        """
        # Create plot config with output_dir
        plot_config = {'output_dir': self.output_dir}
        
        # Call external syntenic_plot function
        syntenic_plot(starts, starts1, ends, ends1, strands, strands1, 
                     types, types1, blocks, cgcid, pulid, plot_config)

class blastp_hit(object):
    def __init__(self,lines):
        self.qseqid = lines[0]
        self.sseqid = lines[1]
        self.pident = lines[2]
        self.length = int(lines[3])
        self.mismatch = int(lines[4])
        self.gapopen  = int(lines[5])
        self.qstart   = int(lines[6])
        self.qend     = int(lines[7])
        self.sstart   = int(lines[8])
        self.send     = int(lines[9])
        self.evalue   = float(lines[10])
        self.bitscore = float(lines[11])
        if len(lines) >= 13:
            self.qlen = int(lines[12])
        if len(lines) >= 14:
            self.slen = int(lines[13])
    
    def __repr__(self):
        return "\t".join([str(self.__dict__[attr]) for attr in self.__dict__])
    
    def __eq__(self, other):
        if self.evalue == self.evalue:
            if self.bitscore == self.bitscore:
                if self.pident == other.pident:
                    return 1
    
    def __le__(self,other):
        if self.evalue > other.evalue:
            return 1
        elif self.evalue == other.evalue:
            if self.bitscore < other.bitscore:
                return 1
            elif self.bitscore == other.bitscore:
                if self.pident < other.pident:
                    return 1


def identity_map(seqsim):
    if 80 <= seqsim <= 100:
        color = "red"
    elif 60 <= seqsim <= 80: 
        color = "blue"
    elif 40 <= seqsim <= 60: 
        color = "green"
    elif 20 <= seqsim <= 40: 
        color = "cyan"
    else:
        color = "gray"
    return color


def Get_Position(starts,ends,strands,maxbp,yshift=0,up=1):
    Width = 1000 ; Height = 160; px = 1/plt.rcParams['figure.dpi']
    poly_heigth = 5
    Triangle_length = 4
    plot_start_x, plot_start_y = [0,Height/2 - poly_heigth-yshift]
    polygens = []
    shfit_pos = starts[0]
    for i in range(len(starts)):
        starts[i] = starts[i] - shfit_pos
        ends[i]   = ends[i] - shfit_pos
    #maxbp = max(ends) - min(starts)  
    pixeachbp =  Width / maxbp  
    ###  5     4
    ###            3         
    ###  1     2

    ###      5     4
    ###  1
    ###      2     3

    blocks = [] ; lines = []
    for i in range(len(starts)):
        if strands[i] == "+":
            positions_str = str( starts[i] * pixeachbp) + " " + str(plot_start_y) + " " ## first point x,y
            positions_str += str( ends[i] * pixeachbp - Triangle_length) + " " + str(plot_start_y) + " "## second point
            if up == 1: ###cluster 1
                blocks.append(positions_str)
            positions_str += str( ends[i] * pixeachbp) + " " + str(plot_start_y + poly_heigth) + " " ## 3
            positions_str += str( ends[i] * pixeachbp - Triangle_length) + " " + str( plot_start_y + 2*poly_heigth) + " " ### 4
            positions_str += str( starts[i] * pixeachbp )+ " " + str(plot_start_y + 2*poly_heigth)
            
            positions_str1 = str( starts[i] * pixeachbp )+ " " + str(plot_start_y + 2*poly_heigth) + " "
            positions_str1 += str( ends[i] * pixeachbp - Triangle_length) + " " + str( plot_start_y + 2*poly_heigth) + " " ### 5
            if up == 2: ### cluster 2
                blocks.append(positions_str1)
        if strands[i] == "-":
            positions_str = str( starts[i] * pixeachbp ) + " " + str(plot_start_y + poly_heigth) + " "
            positions_str += str( starts[i] * pixeachbp + Triangle_length) + " " + str(plot_start_y) + " "
            positions_str += str(ends[i] * pixeachbp) + " " + str(plot_start_y) + " "
            positions_str1 = str(ends[i] * pixeachbp) + " " + str(plot_start_y) + " "
            positions_str1 += str( starts[i] * pixeachbp + Triangle_length) + " " + str(plot_start_y) + " "
            if up == 1:
                blocks.append(positions_str1)
            positions_str += str( ends[i] *pixeachbp ) + " " + str(plot_start_y + 2* poly_heigth) + " "
            positions_str += str( starts[i]* pixeachbp +Triangle_length) + " " + str(plot_start_y + 2* poly_heigth)
            positions_str1 = str( ends[i] *pixeachbp ) + " " + str(plot_start_y + 2* poly_heigth) + " "
            positions_str1 += str( starts[i]* pixeachbp +Triangle_length) + " " + str(plot_start_y + 2* poly_heigth)
            if up == 2:        
                blocks.append(positions_str1)
        #print (positions_str)
        polygens.append(positions_str)
        ### for genome line
        if i < len(starts) -1:
            positions_str = str( ends[i] *pixeachbp) + " " + str(plot_start_y + poly_heigth)  + " "
            positions_str += str( starts[i+1]*pixeachbp) + " " + str(plot_start_y + poly_heigth)
            lines.append(positions_str)
    
    scale_number = 10
    each_scale_bp = maxbp / scale_number
    each_scale_pix = each_scale_bp * pixeachbp
    
    plot_start_y -= 30
    scale_positions = []; scale_positions_texts = [] ; scale_text = []
    scale_positions.append("0 " + str(plot_start_y + 3*poly_heigth) + " " + str(10*each_scale_pix) + " " + str(plot_start_y + 3*poly_heigth))
    
    plot_start_y -= 1
    for i in range(scale_number+1):
        positions_str = str(i*each_scale_pix) + " "
        positions_str += str(plot_start_y + 3* poly_heigth) + " "
        positions_str += str(i*each_scale_pix) + " "
        positions_str += str(plot_start_y + 3*poly_heigth + 0.6* poly_heigth)
        scale_positions.append(positions_str)
        positions_str = str(i*each_scale_pix) + " " + str(plot_start_y + 3*poly_heigth + 0.6* poly_heigth)
        scale_positions_texts.append(positions_str)
        scale_text.append(str(int(each_scale_bp*i)+ shfit_pos))
    #print(scale_positions)
    #print(scale_text)
    return polygens,blocks,lines,scale_positions,scale_text


def plot_Polygon_homologous(polygens1, polygens2, types1, types2, size, ax):
    # colors for different types
    colors_map = {
        "CAZyme": "#E67E22",      # orange
        "TC": "#2ECC71",          # green
        "TF": "#9B59B6",          # purple
        "STP": "#F1C40F",         # golden yellow
        "PEPTIDASE": "#16A085",   # greenish
        "SULFATLAS": "#34495E",   # dark blue
        "Other": "#95A5A6"        # light gray
    }
    default_color = "#95A5A6"  # light gray

    for j in range(len(polygens1)):
        polygen = polygens1[j].split()
        points = []
        color = colors_map.get(types1[j].upper(), default_color) if types1[j] else default_color
        color = colors_map.get(types1[j], default_color)  # keep the original color mapping
        for i in range(int(len(polygen)/2)):
            points.append([float(polygen[2*i]), float(polygen[2*i+1])])
        ax.add_patch(
            Polygon(points, color=color, alpha=0.5, lw=0)
        )

    for j in range(len(polygens2)):
        polygen = polygens2[j].split()
        points = []
        color = colors_map.get(types2[j].upper(), default_color) if types2[j] else default_color
        color = colors_map.get(types2[j], default_color)
        for i in range(int(len(polygen)/2)):
            points.append([float(polygen[2*i]), float(polygen[2*i+1])])
        ax.add_patch(
            Polygon(points, color=color, alpha=0.5, lw=0)
        )


def decode_block(block):
    order1,order2,sim = block.split("-")
    return int(order1),int(order2),float(sim)


def points2(coord):
    x1,y1,x2,y2 = coord.split()
    return float(x1),float(y1),float(x2),float(y2)


def Shade_curve(x11,x12,y11,y12,x21,x22,y21,y22,xmid,ymid,color):
    M, C4, L, CP = Path.MOVETO, Path.CURVE4, Path.LINETO, Path.CLOSEPOLY
    pathdata = [ (M, (x11,y11)),
            (C4, (x11, ymid)),
            (C4, (x21, ymid)),
            (C4, (x21,y21)),
            (L, (x22,y22)),
            (C4, (x22, ymid)),
            (C4, (x12, ymid)),
            (C4, (x12, y12)),
            (CP, (x11,y11))]
    codes, verts = zip(*pathdata)
    path = Path(verts, codes)
    pp = PathPatch(path,color=color,alpha=0.2,lw=0)
    
    return pp

def Shade_Line(x11,x12,y11,y12,x21,x22,y21,y22,xmid,ymid):
    M, C4, L, CP = Path.MOVETO, Path.CURVE4, Path.LINETO, Path.CLOSEPOLY
    pathdata = [(M, (x11,y11)), (L, (x21,y21)), (L, (x22,y22)), (L, (x12,y12)), (CP, (x11,y11))]
    codes, verts = zip(*pathdata)
    path = Path(verts, codes)
    pp = PathPatch(path)
    return pp


def plot_syntenic_block(blocks,blocks1_coor,blocks2_coor,ax):
    for block in blocks:
        order1,order2,sim = decode_block(block)
        coord1 = blocks1_coor[order1]
        coord2 = blocks2_coor[order2]
        x11,y11,x12,y12 = points2(coord1)
        x21,y21,x22,y22 = points2(coord2)
        
        color = identity_map(sim)

        xmid = (x11 + x22 + x21 + x22 )/4
        ymid = (y11 + y22)/2
        
        pp = Shade_curve(x11,x12,y11,y12,x21,x22,y21,y22,xmid,ymid,color)
        #pp = Shade_Line(x11,x12,y11,y12,x21,x22,y21,y22,xmid,ymid)
        ax.add_patch(pp)

def plot_genome_line(lines_coor1,lines_coor2,ax):

    for line in lines_coor1:
        x1,y1,x2,y2 = points2(line)
        ax.add_patch(Polygon([(x1,y1),(x2,y2)], color="gray",lw=2))
    for line in lines_coor2:
        x1,y1,x2,y2 = points2(line)
        ax.add_patch(Polygon([(x1,y1),(x2,y2)], color="gray", lw=2))

### input: gene cluster1: all starts coordinate, end coodinate, strands,Types,
### input: gene cluster2: all starts coordinate, end coodinate, strands,Types,
### input: blocks information,
### the following is an example of inputs

starts = [63045,65290,65974,69423,71382,72920,73962,76286,78315,80633,82094,83100,85536,86051]
ends   = [65270,65826,69375,71255,72806,73918,76238,78289,80606,81952,83077,85391,86048,86809]
strands = "+,+,+,+,+,+,+,+,+,+,+,+,-,-".split(",")
Types = "CAZyme,TF,TC,TC,null,STP,CAZyme,CAZyme,CAZyme,null,CAZyme,CAZyme,null,TC".split(",")


starts1 = [0,2278,3299,4729,7100,8466,9697,10746,12338,15901,17511,20825,22497,
23601,24981,27822,29391,30525,31735,32771,35152,36308]
ends1   = [2271,3247,4631,6934,8450,9654,10720,12330,15101,17143,20805,22466,23592,
24792,27804,29367,30519,31707,32641,35093,36145,36917]
strands1 = "-,-,-,-,-,-,-,-,-,-,+,+,+,+,+,+,+,+,+,+,-,+".split(",")
Types1 = ['CAZyme', 'CAZyme', 'other', 'CAZyme', 'other', 'other', 'other', 'TC', 
'TC', 'other', 'TC', 'TC', 'other', 'other', 'TC', 'TC', 'other', 'other', 'other', 
'CAZyme',"STP","TF"]

### let;s try to give the syntenic block information

### format: cluster order-cluster1 order-sequence indenity

blocks = ["0-0-30",f"1-{len(starts1)-1}-50",f"5-{len(starts1)-2}-84","13-3-65","2-11-76"]

def syntenic_plot(starts, starts1, ends, ends1, strands, strands1, types, types1, blocks, cgcid, pulid, config):
    """
    Generate a syntenic plot for a CGC-PUL pair
    
    Args:
        starts, starts1: Start positions of genes in CGC and PUL
        ends, ends1: End positions of genes in CGC and PUL
        strands, strands1: Strand directions of genes in CGC and PUL
        types, types1: Types of genes in CGC and PUL
        blocks: Homologous block information
        cgcid: CGC ID
        pulid: PUL ID
        config: Configuration dictionary with output_dir key
    """
    # Set up legend elements
    custom_lines = [Line2D([0], [0], color="red", lw=4, alpha=0.5),
                   Line2D([0], [0], color="blue", lw=4, alpha=0.5),
                   Line2D([0], [0], color="green", lw=4, alpha=0.5),
                   Line2D([0], [0], color="cyan", lw=4, alpha=0.5),
                   Line2D([0], [0], color="gray", lw=4, alpha=0.5)]

    labelcolor = ["red", "blue", "green", "cyan", "gray"]
    labels = ["80-100", "60-80", "40-60", "20-40", "0-20"]

    genelabelcolor = [
        "#E67E22",   # CAZyme
        "#2ECC71",   # TC
        "#9B59B6",   # TF
        "#F1C40F",   # STP
        "#16A085",   # PEPTIDASE
        "#34495E",   # SULFATLAS
        "#95A5A6"    # Other
    ]
    geneslabels = [
        "CAZyme",
        "TC",
        "TF",
        "STP",
        "PEPTIDASE",
        "SULFATLAS",
        "Other"
    ]
    genecustom_lines = [Patch(color=c, alpha=0.5) for c in genelabelcolor]

    # Set figure dimensions
    px = 1/plt.rcParams['figure.dpi']
    width = 1600
    height = 320*2
    fig = plt.figure(figsize=(width*px, height*px*2/4))
    ax = fig.add_subplot(111)
    
    # Calculate maximum base pair length
    maxbp = max([max(ends) - min(starts), max(ends1) - min(starts1)])

    # Generate position data for plotting
    polygens, blocks_coor, lines_coor, _, _ = Get_Position(starts, ends, strands, maxbp, yshift=0, up=1)
    polygens1, blocks1_coor, lines_coor1, _, _ = Get_Position(starts1, ends1, strands1, maxbp, yshift=40, up=2)

    # Plot polygons for gene features
    plot_Polygon_homologous(polygens, polygens1, types, types1, 2, ax)
    
    # Plot syntenic relationships
    plot_syntenic_block(blocks, blocks_coor, blocks1_coor, ax)
    
    # Plot genome backbone lines
    plot_genome_line(lines_coor, lines_coor1, ax)

    # Add legends
    legend1 = pyplot.legend(custom_lines, labels, frameon=False, labelcolor=labelcolor,
                          loc='upper left', title="Identity", title_fontsize="x-large")
    ax.add_artist(legend1)

    legend2 = pyplot.legend(genecustom_lines, geneslabels, frameon=False,
                      labelcolor=genelabelcolor, loc='lower left', title="Gene", title_fontsize="x-large")
    ax.add_artist(legend2)

    # Add cluster labels
    plt.text(500, 90, cgcid, fontsize=30, horizontalalignment='center')
    plt.text(500, 0, pulid, fontsize=30, horizontalalignment='center')
    
    # Set axis limits and remove axes
    plt.ylim(0, 100)
    plt.xlim(-100, 1100)
    plt.axis('off')
    ax.plot()
    
    # Layout and save
    plt.tight_layout(pad=0.01)
    cgcid_safe = cgcid.replace("|", "_")  # Replace pipe characters for filesystem safety
    
    # Ensure output directory exists
    pdf_dir = os.path.join(config['output_dir'], "synteny_pdf")
    os.makedirs(pdf_dir, exist_ok=True)
    
    # Save figure to file
    plt.savefig(f"{pdf_dir}/{cgcid_safe}-syntenic.pdf")
    plt.close()
    
def read_blast_result_cgc(filename):
    querydict = {}
    for line in open(filename):
        lines = line.split()
        queryids = lines[0].split("|")
        queryid = queryids[0] + "|" + queryids[1]
        hitpulid = lines[1].split(":")[1]
        #lines[0] = lines[2]+"|"+lines[0] ### cgc id
        querydict.setdefault(queryid+":"+hitpulid,[]).append(blastp_hit(lines))
        #print (lines[0])  
    return querydict

def syntenic_plot_allpairs(args):
    error = 0
    cgcpul_blastp = read_blast_result_cgc(args.blastp)
    
    cgc_proteinid2gene,cgcid2gene,cgcid2geneid = read_UHGG_CGC_stanrdard_out(args.cgc)
    PULid_proteinid2gene,PULid2gene,PULid2geneid = read_PUL_cgcgff(args)
    
    pdf_directory = os.path.join(args.output_dir, "synteny_pdf")    
    os.makedirs(pdf_directory, exist_ok=True)
    for line in open(args.input_sub_out).readlines()[1:]: ### for each pairs
        lines = line.rstrip().split("\t")
        cgc = lines[0]
        pul = lines[1]
        if not pul: ### no homologous cgc found from dbCAN-PUL
            continue
        
        cgcpul = cgc+":"+pul
        bed_cgc = cgcid2gene[cgc]
        bed_pul = PULid2gene[pul]

        starts1,ends1,strands1,types1 = Get_parameters_for_plot(bed_cgc)
        starts2,ends2,strands2,types2 = Get_parameters_for_plot(bed_pul)
        genes1 = cgcid2geneid[cgc]
        genes2 = PULid2geneid[pul]
        # print (cgc,pul)
        # print (genes1)
        # print (genes2)
        blocks = []
        for record in cgcpul_blastp[cgcpul]: ### generate block information
            query = record.qseqid
            hit   = record.sseqid
            cgc_proteinid = query.split("|")[2]
            pul_proteinid = hit.split(":")[2]
            # if not pul_proteinid:
            #     pul_proteinid = hit.split(":")[3]
            try:
                index1 = genes1.index(cgc_proteinid)
                index2 = genes2.index(pul_proteinid)
                blocks.append(f"{index1}-{index2}-{record.pident}")
            except:
                # print (cgcpul,query,hit,cgc_proteinid,pul_proteinid,genes1,genes2)
                continue
            #print (cgc_proteinid2gene[cgc_proteinid],pul_proteinid2gene[pul_proteinid])
        syntenic_plot(starts1,starts2,ends1,ends2,strands1,strands2,types1,types2,blocks,cgc,pul,args)

def Get_parameters_for_plot(CGC_stanrdard_list):
    starts = [] ; ends = [] ; strands = []; types = []
    for gene in CGC_stanrdard_list:
        starts.append(gene.Gene_Start)
        ends.append(gene.Gene_END)
        strands.append(gene.Strand)
        types.append(gene.Gene_Type)
    return starts,ends,strands,types

def attribution(desc,mark="DB"):
    descs = desc.split(";")
    for des in descs:
        if des.startswith(mark):
            return des.split("=")[-1]
    return "Other"

def trim_proteinid(proteinid):
    clear_id = proteinid.split(".")
    return clear_id[0]
    #return ".".join(clear_id[0:-1])


def read_cgcgff(filename, geneid2gene):
    if not os.path.exists(filename):
        return None
    for line in open(filename):
        lines = line.rstrip("\n").split("\t")
        #  protein_id
        desc = lines[-1]
        protein_id = None
        feature = "Other"
        for item in desc.split(";"):
            if item.startswith("protein_id="):
                protein_id = item.split("=", 1)[1]
            elif item.startswith("CGC_annotation="):
                ann = item.split("=", 1)[1]
                # set + as delimiter and split by |
                feature = ann.split("+")[0].split("|")[0] if "|" in ann else ann.split("+")[0]
        if protein_id is None:
            continue  # skip if no protein_id found
        #protein_id = trim_proteinid(protein_id)
        PULid = filename.split('.')[0].split("/")[-1]
        # gene information 
        newline = [
            PULid,                # CGCID
            feature,              # Gene_Type
            lines[0],             # Contig_ID
            protein_id,           # Protein_ID
            lines[3],             # Gene_Start
            lines[4],             # Gene_END
            lines[6],             # Strand
            ann               # Protein_Family
        ]
        geneid2gene[PULid + ":" + protein_id] = CGC_stanrdard(newline)

def read_PUL_cgcgff(args):
    PULidgeneid2gene = {} ### "PUL00035:BT_3559"
    ### list all files

    dbCANPUL_data_folder = args.db_dir if args.db_dir.endswith("/") else args.db_dir+"/"
    if not dbCANPUL_data_folder.startswith("/"): ### absolute path
        dbCANPUL_data_folder = os.path.abspath(dbCANPUL_data_folder) + "/dbCAN-PUL/"
    else:
        dbCANPUL_data_folder += "/dbCAN-PUL/"
    folder = os.listdir(dbCANPUL_data_folder)
    for filename in folder:
        if filename.startswith("PUL") and filename.endswith(".out"): ### 
            read_cgcgff(dbCANPUL_data_folder+filename+"/cgc.gff",PULidgeneid2gene)
    
    cgcid2gene = {}; cgcid2geneid ={}
    for PULidgeneid in PULidgeneid2gene:
        gene = PULidgeneid2gene[PULidgeneid]
        cgcid2gene.setdefault(gene.CGCID,[]).append(gene)
        cgcid2geneid.setdefault(gene.CGCID,[]).append(gene.Protein_ID)
    
    return PULidgeneid2gene,cgcid2gene,cgcid2geneid

class CGC_stanrdard(object):
    def __init__(self,lines):
        self.CGCID = lines[0]
        self.Gene_Type = lines[1]
        self.Contig_ID =  lines[2]
        self.Protein_ID  = lines[3]
        self.Gene_Start  = int(lines[4])
        self.Gene_END   = int(lines[5])
        self.Strand = lines[6]
        self.Protein_Family = lines[7]
    def __repr__(self):
        return "\t".join([str(self.__dict__[attr]) for attr in self.__dict__])

def read_UHGG_CGC_stanrdard_out(filename):
    geneid2gene = {}; cgcid2gene = {};cgcid2geneid = {}
    
    for line in open(filename):
        if line[0:4] == "CGC#":
            continue
        lines = line.rstrip("\n").split("\t")
        lines[0] = lines[2] +"|" +lines[0]
        gene = CGC_stanrdard(lines)
        geneid2gene[gene.Protein_ID] = gene
        cgcid2gene.setdefault(gene.CGCID,[]).append(gene)
        cgcid2geneid.setdefault(gene.CGCID,[]).append(gene.Protein_ID)
    return geneid2gene,cgcid2gene,cgcid2geneid


# def parse_argv():
#     parser = argparse.ArgumentParser(description='syntenic plot for two homologous CGC')
#     parser.add_argument('function', help='what function will be used to analyze.')
#     parser.add_argument('-i','--input',help='cgc_finder output')
#     parser.add_argument('-b','--blastp',help='blastp result for cgc')
#     parser.add_argument('--cgc')
#     parser.add_argument('--pul')
#     parser.add_argument('--db_dir', default="db", help='Database directory')
#     args = parser.parse_args()
#     return args


# def main():
#     args = parse_argv()
#     syntenic_plot_allpairs(args)

# if __name__ == "__main__":
#     args = parse_argv()
#     #syntenic_plot(starts,starts1,ends,ends1,strands,strands1,Types,Types1,blocks)
#     if args.function == "syntenic_plot":
#         ## python3 /home/jinfang/libsvm_practise/syntenic.plot.py syntenic_plot -b PUL.pairs.blastp -i cgc.pul.hits.strict --cgc cgc_standard.out
#         ## python3 syntenic.plot.py syntenic_plot -b PUL.pairs.blastp --cgc cgc_standard.out -i cgc.pul.hits.strict --pul PUL.out
#         ## python3 syntenic.plot.py syntenic_plot -b PUL.pairs.blastp --cgc cgc_standard.out -i cgc.pul.hits.strict
#         ## python3 /array1/www/dbCAN3/ty/syntenic.plot.py syntenic_plot -b blastp.out --cgc cgc_standard.out -i sub.prediction.out
#         syntenic_plot_allpairs(args)
