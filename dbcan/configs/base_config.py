from dataclasses import dataclass, field, fields
from typing import Optional, List, Tuple

from dbcan.constants.base_constants import GFF_MIN_MATCH_DEFAULT


@dataclass
class BaseConfig:
    @staticmethod
    def from_dict(config_class, config_dict):
        field_names = {f.name for f in fields(config_class)}
        filtered_dict = {k: v for k, v in config_dict.items() if k in field_names}
        return config_class(**filtered_dict)


@dataclass
class GeneralConfig(BaseConfig):
    input_raw_data: str = None
    output_dir: str = None
    mode: str = None
    db_dir: str = None
    # Worker threads for gene prediction (pyrodigal); None or <= 0 uses all CPUs.
    threads: Optional[int] = None

    #thread



@dataclass
class OverviewGeneratorConfig(BaseConfig):
    output_dir: str
    # Methods run in this CAZyme_annotation call (--methods). Only their result files are merged;
    # None merges every result file found (diamond.out, dbCAN_hmm_results.tsv, dbCANsub_hmm_results.tsv).
    methods: Optional[Tuple[str, ...]] = None


@dataclass
class GFFConfig(BaseConfig):
    output_dir: str
    input_gff: str
    gff_type: str
    # Prodigal: avoid BCBio GFF.parse on huge files (see annotate_prodigal_gff_streaming).
    prodigal_gff_streaming: str = "auto"  # auto | on | off
    prodigal_streaming_threshold_mb: int = 50
    # Minimum fraction of protein-coding GFF genes that must match a protein of uniInput.faa.
    gff_min_match: float = GFF_MIN_MATCH_DEFAULT




@dataclass
class CGCPlotConfig(BaseConfig):
    output_dir: str
    layout: str = "auto"   # auto | circos | overview (see plot_cgc_circle_constants)


@dataclass
class SyntenicPlotConfig:
    output_dir: str
    db_dir: str
    input_sub_out: Optional[str] = None
    blastp: Optional[str] = None
    cgc: Optional[str] = None


def create_config(config_class, **kwargs):
    return config_class.from_dict(config_class, kwargs)






