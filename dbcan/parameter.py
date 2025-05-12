import rich_click as click
from dataclasses import dataclass, fields
from typing import Optional
import psutil

@dataclass
class BaseConfig:
    @staticmethod
    def from_dict(config_class, config_dict):
        field_names = {f.name for f in fields(config_class)}
        filtered_dict = {k: v for k, v in config_dict.items() if k in field_names}
        return config_class(**filtered_dict)

@dataclass
class GeneralConfig(BaseConfig):
    input_raw_data: str
    output_dir: str
    mode: str

@dataclass
class DBDownloaderConfig(BaseConfig):
    db_dir: str

@dataclass
class DiamondConfig(BaseConfig):
    db_dir: str
    threads: int
    output_dir: str

    e_value_threshold: float = 1e-102
    verbose_option: bool = False


@dataclass
class DiamondTCConfig(BaseConfig):
    db_dir: str
    threads: int
    output_dir: str

    e_value_threshold_tc: float = 1e-4
    coverage_threshold_tc: float = 0.35
    verbose_option: bool = False

@dataclass
class DiamondSulfataseConfig(BaseConfig):
    db_dir: str
    threads: int
    output_dir: str

    e_value_threshold_tc: float = 1e-4
    coverage_threshold_tc: float = 0.35
    verbose_option: bool = False


@dataclass
class DiamondPeptidaseConfig(BaseConfig):
    db_dir: str
    threads: int
    output_dir: str

    e_value_threshold_tc: float = 1e-4
    coverage_threshold_tc: float = 0.35
    verbose_option: bool = False



@dataclass
class PyHMMERConfig(BaseConfig):
    db_dir: str
    threads: int
    output_dir: str

    e_value_threshold_dbcan: float = 1e-15
    coverage_threshold_dbcan: float = 0.35


@dataclass
class DBCANSUBProcessorConfig(PyHMMERConfig):
    db_dir: str
    e_value_threshold_dbsub: float = 1e-15
    coverage_threshold_dbsub: float = 0.35



@dataclass
class OverviewGeneratorConfig(BaseConfig):
    output_dir: str


@dataclass
class PyHMMERTFConfig(PyHMMERConfig):
    db_dir: str
    e_value_threshold_tf: float  = 1e-4
    coverage_threshold_tf: float = 0.35


@dataclass
class PyHMMERSTPConfig(PyHMMERConfig):
    db_dir: str
    e_value_threshold_stp: float  = 1e-4
    coverage_threshold_stp: float  = 0.35



@dataclass
class GFFConfig(BaseConfig):
    output_dir: str
    input_gff: str
    gff_type: str

@dataclass
class CGCFinderConfig(BaseConfig):
    output_dir: str
    additional_genes: list

    num_null_gene: int =2
    base_pair_distance: int =15000
    use_null_genes: bool = True
    use_distance: bool = False

@dataclass
class HomologyParameters(BaseConfig):
    upghn: int  # uniq_pul_gene_hit_num
    uqcgn: int  # uniq_query_cgc_gene_num
    cpn: int    # CAZyme_pair_num
    tpn: int    # total_pair_num
    identity_cutoff: float
    coverage_cutoff: float
    bitscore_cutoff: float
    evalue_cutoff: float
    extra_pair_type: str = None
    extra_pair_type_num: str = None

    def get_extra_pairs(self):
        ept = self.extra_pair_type.split(",") if self.extra_pair_type else None
        eptn = self.extra_pair_type_num.split(",") if ept else None

        if ept and eptn and len(ept) != len(eptn):
            raise ValueError(f"({len(ept)}){len(eptn)}) extra_pair_type and extra_pair_type_num must have the same length.")

        return ept, eptn

@dataclass
class DBCANSubParameters(BaseConfig):
    hmmevalue: float
    hmmcov: float
    num_of_protein_substrate_cutoff: int
    num_of_domains_substrate_cutoff: int
    substrate_scors: float

@dataclass
class CGCSubstrateConfig(BaseConfig):
    output_dir: str
    workdir: str
    odbcan_sub: str
    db_dir: str
    output_dir: str
    pul: str
    out: str
    rerun: bool
    odbcanpul: bool

    uniq_pul_gene_hit_num: int=2
    uniq_query_cgc_gene_num: int=2
    CAZyme_pair_num: int=1
    total_pair_num: int=2
    extra_pair_type: str=None
    extra_pair_type_num: str="0"
    identity_cutoff: float=0
    coverage_cutoff: float=0
    bitscore_cutoff: float=50
    evalue_cutoff: float=0.01
    hmmcov: float=0
    hmmevalue: float=0.01
    num_of_domains_substrate_cutoff: int=2
    num_of_protein_substrate_cutoff: int=2
    substrate_scors: int=2


    def get_homology_params(self) -> HomologyParameters:

        return HomologyParameters(
            upghn=self.uniq_pul_gene_hit_num,
            uqcgn=self.uniq_query_cgc_gene_num,
            cpn=self.CAZyme_pair_num,
            tpn=self.total_pair_num,
            identity_cutoff=self.identity_cutoff,
            coverage_cutoff=self.coverage_cutoff,
            bitscore_cutoff=self.bitscore_cutoff,
            evalue_cutoff=self.evalue_cutoff,
            extra_pair_type=self.extra_pair_type,
            extra_pair_type_num=self.extra_pair_type_num
        )

    def get_dbsub_params(self) -> DBCANSubParameters:

        return DBCANSubParameters(
            hmmevalue=self.hmmevalue,
            hmmcov=self.hmmcov,
            num_of_protein_substrate_cutoff=self.num_of_protein_substrate_cutoff,
            num_of_domains_substrate_cutoff=self.num_of_domains_substrate_cutoff,
            substrate_scors=self.substrate_scors
        )

@dataclass
class SynPlotConfig(BaseConfig):
    db_dir: str
    output_dir: str


@dataclass
class CGCPlotConfig(BaseConfig):
    output_dir: str






def create_config(config_class, **kwargs):
    return config_class.from_dict(config_class, kwargs)

# Define shared options
output_dir_option = click.option('--output_dir', required=True, help='Directory for the output files')
threads_option = click.option('--threads',  type=int, help='Number of threads', default=psutil.cpu_count())

methods_option = click.option('--methods',
    default=['diamond', 'hmm', 'dbCANsub'],
    help='Specify the annotation methods to use. Options are diamond, hmm, and dbCANsub.',
    multiple=True)


# Define group options
def general_options(func):
    func = click.option('--input_raw_data', required=True, help='Path to the input raw data')(func)
    func = output_dir_option(func)
    func = click.option('--mode', default='prok', required=True, help='Mode of input sequence')(func)

    return func

def database_options(func):
    func = click.option('--db_dir', required=True, help='Directory for the database')(func)
    return func

def diamond_options(func):
    func = click.option('--e_value_threshold', type=float, help='E-value threshold for diamond', default=1e-102 )(func)
    func = click.option('--verbose_option', is_flag=True, help='Enable verbose option for diamond', default=False)(func)
    func = output_dir_option(func)
    return func

def diamond_tc_options(func):
    func = click.option('--e_value_threshold_tc', type=float, help='E-value threshold for TC' ,default=1e-4)(func)
    func = click.option('--coverage_threshold_tc', type=float, help='Coverage threshold for TC', default=0.35)(func)
    return func

def pyhmmer_dbcan_options(func):
    func = click.option('--e_value_threshold_dbcan',  type=float, help='E-value threshold for HMMER',  default=1e-15)(func)
    func = click.option('--coverage_threshold_dbcan',  type=float, help='Coverage threshold for HMMER', default=0.35)(func)
    func = output_dir_option(func)
    return func

def dbcansub_options(func):
    func = click.option('--e_value_threshold_dbsub',  type=float, help='E-value threshold for dbCAN-sub HMMER', default=1e-15)(func)
    func = click.option('--coverage_threshold_dbsub',  type=float, help='Coverage threshold for dbCAN-sub HMMER', default=0.35)(func)
    func = output_dir_option(func)
    return func

def pyhmmer_tf(func):
    func = click.option('--e_value_threshold_tf',  type=float, help='E-value threshold for TF HMMER', default=1e-4)(func)
    func = click.option('--coverage_threshold_tf',  type=float, help='Coverage threshold for TF HMMER', default=0.35)(func)
    func = output_dir_option(func)
    return func

def pyhmmer_stp(func):
    func = click.option('--e_value_threshold_stp',  type=float, help='E-value threshold for STP HMMER',default=1e-4)(func)
    func = click.option('--coverage_threshold_stp',  type=float, help='Coverage threshold for STP HMMER',default=0.35)(func)
    func = output_dir_option(func)
    return func


def cgc_gff_option(func):
    func = click.option('--input_gff', required=True, help='input GFF file')(func)
    func = click.option('--gff_type', required=True, help='GFF file type')(func)
    func = output_dir_option(func)
    return func

def cgc_options(func):
    func = click.option('--additional_genes', multiple=True, default=["TC"], help='Specify additional gene types for CGC annotation, including TC, TF, and STP')(func)
    func = click.option('--num_null_gene', type=int, default=2, help='Maximum number of null genes allowed between signature genes.')(func)
    func = click.option('--base_pair_distance', type=int, default=15000, help='Base pair distance of signature genes.')(func)
    func = click.option('--use_null_genes/--no-use_null_genes', is_flag=True, default=True, help='Use null genes in CGC annotation.')(func)
    func = click.option('--use_distance', is_flag=True, default=False, help='Use base pair distance in CGC annotation.')(func)
    func = output_dir_option(func)
    return func

def cgc_substrate_base_options(func):
    """base opiton"""
    func = output_dir_option(func)
    func = click.option('--pul', help="dbCAN-PUL PUL.faa")(func)
    func = click.option('-o', '--out', default="substrate.out", help="substrate prediction result")(func)
    func = click.option('-w', '--workdir', default=".", type=str, help="work directory")(func)
    func = click.option('-rerun', '--rerun', default=False, type=bool, help="re run the prediction")(func)
    func = click.option('-env', '--env', default="local", type=str, help="run environment")(func)
    func = click.option('-odbcan_sub', '--odbcan_sub', help="export dbcan-sub sub result")(func)
    func = click.option('-odbcanpul', '--odbcanpul', default=True, type=bool, help="export dbcan pul sub result")(func)
    func = click.option('--db_dir', default='./dbCAN_databases', required=True, help='database folder')(func)
    return func

def cgc_substrate_homology_params_options(func):
    """dbCAN-PUL approach homology parameters"""
    func = click.option('-upghn', '--uniq_pul_gene_hit_num', default=2, type=int, help="num of uniq gene hit of pul")(func)
    func = click.option('-uqcgn', '--uniq_query_cgc_gene_num', default=2, type=int, help="num of uniq gene hit of cgc")(func)
    func = click.option('-cpn', '--CAZyme_pair_num', default=1, type=int, help="num of CAZyme")(func)
    func = click.option('-tpn', '--total_pair_num', default=2, type=int, help="total pair number")(func)
    func = click.option('-ept', '--extra_pair_type', default=None, type=str, help="extra pair type")(func)
    func = click.option('-eptn', '--extra_pair_type_num', default="0", type=str, help="extra pair number")(func)
    func = click.option('-iden', '--identity_cutoff', default=0.0, type=float, help="identity ")(func)
    func = click.option('-cov', '--coverage_cutoff', default=0.0, type=float, help="coverage ")(func)
    func = click.option('-bsc', '--bitscore_cutoff', default=50.0, type=float, help="bit score")(func)
    func = click.option('-evalue', '--evalue_cutoff', default=0.01, type=float, help="evalue")(func)
    return func

def cgc_substrate_dbcan_sub_param_options(func):
    """dbCAN-sub substrate prediction parameters"""
    func = click.option('-hmmcov', '--hmmcov', default=0.0, type=float, help="hmm coverage")(func)
    func = click.option('-hmmevalue', '--hmmevalue', default=0.01, type=float, help="HMM evalue")(func)
    func = click.option('-ndsc', '--num_of_domains_substrate_cutoff', default=2, type=int, help="num of domains substrate")(func)
    func = click.option('-npsc', '--num_of_protein_substrate_cutoff', default=2, type=int, help="num of protein substrate")(func)
    func = click.option('-subs', '--substrate_scors', default=2, type=int, help="substrate score")(func)
    return func

def cgc_sub_options(func):
    """total option for cgc substrate prediction"""
    func = cgc_substrate_base_options(func)
    func = cgc_substrate_homology_params_options(func)
    func = cgc_substrate_dbcan_sub_param_options(func)
    return func

def syn_plot_options(func):
    func = click.option('--db_dir', required=True, help='Path to the database directory')(func)
    func = output_dir_option(func)
    return func

def cgc_circle_plot_options(func):
    func = output_dir_option(func)
    return func






