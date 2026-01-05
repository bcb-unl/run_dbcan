from dataclasses import dataclass, field, fields
from typing import Optional, List


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

    #thread



@dataclass
class OverviewGeneratorConfig(BaseConfig):
    output_dir: str


@dataclass
class GFFConfig(BaseConfig):
    output_dir: str
    input_gff: str
    gff_type: str




@dataclass
class CGCPlotConfig(BaseConfig):
    output_dir: str


@dataclass
class SyntenicPlotConfig:
    output_dir: str
    db_dir: str
    input_sub_out: Optional[str] = None
    blastp: Optional[str] = None
    cgc: Optional[str] = None


def create_config(config_class, **kwargs):
    return config_class.from_dict(config_class, kwargs)






