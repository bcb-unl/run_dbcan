from dataclasses import dataclass
from typing import List, Optional

@dataclass
class PlotsConfig:
    # Common
    input_dir: str
    db_dir: str = "db"
    output: str = "output"

    # CGC specific
    cgcid: Optional[str] = None
    reads_count: Optional[str] = None
    bedtools: Optional[str] = None

    # Bar/Heatmap specific
    samples: Optional[List[str]] = None
    top: int = 20
    plot_style: str = "ggplot"
    vertical_bar: bool = False
    show_fig: bool = False
    show_abund: bool = False
    palette: Optional[str] = None
    cluster_map: bool = False
    filter_col: Optional[str] = None
    filter_value: Optional[str] = None
    pdf: str = "bar_plot.pdf"
