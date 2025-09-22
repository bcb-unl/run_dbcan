from dataclasses import dataclass
from typing import Optional

# ...existing code...

@dataclass
class AbundanceConfig:
    # 统一入口/输出目录
    input_dir: str                     # run_dbcan 的输出目录
    bedtools_depth: str                # bedtools/自定义 cal_coverage 结果文件
    output_dir: str = "output"         # 结果输出目录（部分功能直接写固定文件名时仍写 CWD）
    # 可选输入
    gff: Optional[str] = None          # cal_coverage 使用
    R1: Optional[str] = None
    R2: Optional[str] = None
    db_dir: str = "db"
    # 归一化方法
    abundance: str = "RPM"             # FPKM/RPM/TPM
    # 比对过滤参数
    overlap_base_ratio: float = 0.2
    mapping_quality: int = 30
    identity: float = 0.98
    threads: int = 1
    hifi: bool = False
