from dbcan.configs.base_config import BaseConfig
from dataclasses import dataclass, field
from typing import Optional, List
import dbcan.constants.cgcfinder_constants as C

@dataclass
class CGCFinderConfig(BaseConfig):
    output_dir: str
    additional_genes: List[str] = field(default_factory=list)

    gff_file: Optional[str] = None

    # thresholds and parameters for CGC identification
    num_null_gene: int = C.CGC_DEFAULT_NULL_GENE
    base_pair_distance: int = C.CGC_DEFAULT_BP_DISTANCE
    use_null_genes: bool = True
    use_distance: bool = False

    # extension behavior
    extend_mode: str = C.CGC_DEFAULT_EXTEND_MODE  # 'none' | 'bp' | 'gene'
    extend_bp: int = C.CGC_DEFAULT_EXTEND_BP
    extend_gene_count: int = C.CGC_DEFAULT_EXTEND_GENE_COUNT


    # additional criteria for CGC identification
    additional_logic: str = "all"  # 'all' | 'any'
    additional_min_categories: int = 1  # minimum number of different additional gene categories that must be present when using 'any' logic

    # New: minimum core CAZyme count, minimum cluster gene count, GFF feature filtering
    min_core_cazyme: int = C.CGC_DEFAULT_MIN_CORE_CAZYME
    min_cluster_genes: int = C.CGC_DEFAULT_MIN_CLUSTER_GENES
    feature_types: List[str] = field(default_factory=lambda: ["gene", "CDS"])  # GFF feature types to consider
