from dataclasses import dataclass
from dbcan.configs.base_config import GeneralConfig

@dataclass
class DBDownloaderConfig(GeneralConfig):
    db_dir: str
    cgc: bool = True
    # align with DBDownloader robust networking
    timeout: int = 30
    retries: int = 3
    no_overwrite: bool = False
    resume: bool = True
    verify_ssl: bool = True
    aws_s3: bool = False
