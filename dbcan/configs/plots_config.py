from dataclasses import dataclass, field
from typing import List, Optional, Tuple

import dbcan.constants.cgc_substrate_prediction_constants as CSP
import dbcan.constants.plots_constants as P


@dataclass
class PlotsConfig:
    # Common
    input_dir: str
    db_dir: str = "db"
    output_dir: str = "."
    formats: Tuple[str, ...] = P.PLOT_DEFAULT_FORMATS
    dpi: int = P.PLOT_DEFAULT_DPI

    # CGC specific
    cgcid: Optional[str] = None
    pul: Optional[str] = None
    reads_count: List[str] = field(default_factory=list)
    show_protein_ids: bool = False
    synteny_revcomp: bool = True
    evalue_cutoff: float = CSP.DEFAULT_EVALUE_CUTOFF
    identity_cutoff: float = CSP.DEFAULT_IDENTITY_CUTOFF
    coverage_cutoff: float = CSP.DEFAULT_COVERAGE_CUTOFF

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
    pdf: Optional[str] = None
