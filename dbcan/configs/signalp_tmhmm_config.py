from dbcan.configs.base_config import BaseConfig
from dataclasses import dataclass
from typing import Optional
import dbcan.constants.signalp_tmhmm_constants as S

@dataclass
class SignalPTMHMMConfig(BaseConfig):
    output_dir: str


    run_signalp: bool = False
    run_deeptmhmm: bool = False

    #parameters for SignalP
    signalp_org: str = S.DEFAULT_SIGNALP_ORG
    signalp_mode: str = S.DEFAULT_MODE
    signalp_format: str = S.DEFAULT_FORMAT

    # parameters for local DeepTMHMM wrapper
    deeptmhmm_dir: Optional[str] = None
    deeptmhmm_python: str = S.DEFAULT_DEEPTMHMM_PYTHON

    #general parameters
    force: bool = False
    threads: Optional[int] = None