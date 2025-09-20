from dataclasses import dataclass
from typing import Optional
from dbcan.configs.base_config import BaseConfig
import dbcan.constants.cgc_substrate_prediction_constants as C


@dataclass
class HomologyParameters(BaseConfig):
    upghn: int = C.DEFAULT_UNIQ_PUL_GENE_HIT_NUM
    uqcgn: int = C.DEFAULT_UNIQ_QUERY_CGC_GENE_NUM
    cpn: int = C.DEFAULT_CAZYME_PAIR_NUM
    tpn: int = C.DEFAULT_TOTAL_PAIR_NUM
    identity_cutoff: float = C.DEFAULT_IDENTITY_CUTOFF
    coverage_cutoff: float = C.DEFAULT_COVERAGE_CUTOFF
    bitscore_cutoff: float = C.DEFAULT_BITSCORE_CUTOFF
    evalue_cutoff: float = C.DEFAULT_EVALUE_CUTOFF
    extra_pair_type: Optional[str] = C.DEFAULT_EXTRA_PAIR_TYPE
    extra_pair_type_num: Optional[str] = C.DEFAULT_EXTRA_PAIR_TYPE_NUM

    def get_extra_pairs(self):
        ept = self.extra_pair_type.split(",") if self.extra_pair_type else None
        eptn = self.extra_pair_type_num.split(",") if ept and self.extra_pair_type_num else None
        if ept and eptn and len(ept) != len(eptn):
            raise ValueError(f"({len(ept)})({len(eptn)}) extra_pair_type and extra_pair_type_num must have the same length.")
        return ept, eptn


@dataclass
class DBCANSubParameters(BaseConfig):
    hmmevalue: float = C.DEFAULT_SUB_HMM_EVALUE
    hmmcov: float = C.DEFAULT_SUB_HMM_COVERAGE
    num_of_protein_substrate_cutoff: int = C.DEFAULT_SUB_PROTEIN_SUBSTRATE_CUTOFF
    num_of_domains_substrate_cutoff: int = C.DEFAULT_SUB_DOMAIN_SUBSTRATE_CUTOFF
    substrate_scors: float = C.DEFAULT_SUBSTRATE_SCORE_THRESHOLD


@dataclass
class CGCSubstrateConfig(BaseConfig):
    output_dir: str
    db_dir: str
    threads: int = C.DEFAULT_THREADS
    odbcan_sub: bool = C.DEFAULT_OUTPUT_DBCAN_SUB
    odbcanpul: bool = C.DEFAULT_OUTPUT_DBCAN_PUL

    # Homology thresholds
    uniq_pul_gene_hit_num: int = C.DEFAULT_UNIQ_PUL_GENE_HIT_NUM
    uniq_query_cgc_gene_num: int = C.DEFAULT_UNIQ_QUERY_CGC_GENE_NUM
    CAZyme_pair_num: int = C.DEFAULT_CAZYME_PAIR_NUM
    total_pair_num: int = C.DEFAULT_TOTAL_PAIR_NUM
    extra_pair_type: Optional[str] = C.DEFAULT_EXTRA_PAIR_TYPE
    extra_pair_type_num: Optional[str] = C.DEFAULT_EXTRA_PAIR_TYPE_NUM
    identity_cutoff: float = C.DEFAULT_IDENTITY_CUTOFF
    coverage_cutoff: float = C.DEFAULT_COVERAGE_CUTOFF
    bitscore_cutoff: float = C.DEFAULT_BITSCORE_CUTOFF
    evalue_cutoff: float = C.DEFAULT_EVALUE_CUTOFF

    # dbCAN-sub thresholds
    hmmevalue: float = C.DEFAULT_SUB_HMM_EVALUE
    hmmcov: float = C.DEFAULT_SUB_HMM_COVERAGE
    num_of_domains_substrate_cutoff: int = C.DEFAULT_SUB_DOMAIN_SUBSTRATE_CUTOFF
    num_of_protein_substrate_cutoff: int = C.DEFAULT_SUB_PROTEIN_SUBSTRATE_CUTOFF
    substrate_scors: float = C.DEFAULT_SUBSTRATE_SCORE_THRESHOLD



@dataclass
class SynPlotConfig(BaseConfig):
    db_dir: str
    output_dir: str
