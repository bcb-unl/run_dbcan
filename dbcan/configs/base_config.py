import rich_click as click
from dataclasses import dataclass, field, fields
from typing import Optional, List

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
    return func

def diamond_tc_options(func):
    func = click.option('--e_value_threshold_tc', type=float, help='E-value threshold for TC' ,default=1e-4)(func)
    func = click.option('--coverage_threshold_tc', type=float, help='Coverage threshold for TC', default=0.35)(func)
    return func
def diamond_tf_options(func):
    func = click.option('--e_value_threshold_tf', type=float, help='E-value threshold for TF' ,default=1e-4)(func)
    func = click.option('--coverage_threshold_tf', type=float, help='Coverage threshold for TF', default=0.35)(func)
    func = click.option('--prokaryotic/--no-prokaryotic', is_flag=True, help='Enable prokaryotic mode for TF', default=True)(func)
    return func

def pyhmmer_dbcan_options(func):
    func = click.option('--e_value_threshold_dbcan',  type=float, help='E-value threshold for HMMER',  default=1e-15)(func)
    func = click.option('--coverage_threshold_dbcan',  type=float, help='Coverage threshold for HMMER', default=0.35)(func)
    return func

def dbcansub_options(func):
    func = click.option('--e_value_threshold_dbsub',  type=float, help='E-value threshold for dbCAN-sub HMMER', default=1e-15)(func)
    func = click.option('--coverage_threshold_dbsub',  type=float, help='Coverage threshold for dbCAN-sub HMMER', default=0.35)(func)
    return func

def pyhmmer_tf(func):
    func = click.option('--e_value_threshold_tf',  type=float, help='E-value threshold for TF HMMER', default=1e-4)(func)
    func = click.option('--coverage_threshold_tf',  type=float, help='Coverage threshold for TF HMMER', default=0.35)(func)
    func = click.option('--fungi/--no-fungi', is_flag=True, help='Enable fungi mode for TF HMMER', default=False)(func)
    return func

def pyhmmer_stp(func):
    func = click.option('--e_value_threshold_stp',  type=float, help='E-value threshold for STP HMMER',default=1e-4)(func)
    func = click.option('--coverage_threshold_stp',  type=float, help='Coverage threshold for STP HMMER',default=0.35)(func)
    return func

def pyhmmer_pfam(func):
    func = click.option('--run_pfam', help='Run Pfam HMMER for CGC null gene annotation', is_flag=True, default=False)(func)
    func = click.option('--e_value_threshold_pfam',  type=float, help='E-value threshold for Pfam HMMER',default=1e-4)(func)
    func = click.option('--coverage_threshold_pfam',  type=float, help='Coverage threshold for Pfam HMMER',default=0.35)(func)
    func = click.option('--null_from_gff', is_flag=True, default=False,
                        help='Extract null genes from cgc.gff instead of cgc_standard_out.tsv')(func)  # 新增
    return func

def cgc_gff_option(func):
    func = click.option('--input_gff', required=True, help='input GFF file')(func)
    func = click.option('--gff_type', required=True, help='GFF file type')(func)
    return func

def cgc_options(func):
    func = click.option('--additional_genes', multiple=True, default=["TC"], help='Specify additional gene types for CGC annotation, including TC, TF, and STP')(func)
    func = click.option('--additional_logic',
                        type=click.Choice(['all','any']),
                        default='all',
                        help="Logic for multiple --additional_genes: 'all' requires all present; 'any' requires at least one.")(func)
    func = click.option('--additional_min_categories',
                        type=int, default=1,
                        help="When --additional_logic=any, require at least this number of distinct additional categories.")(func)
    func = click.option('--num_null_gene', type=int, default=2, help='Maximum number of null genes allowed between signature genes.')(func)
    func = click.option('--base_pair_distance', type=int, default=15000, help='Base pair distance of signature genes.')(func)
    func = click.option('--use_null_genes/--no-use_null_genes', is_flag=True, default=True, help='Use null genes in CGC annotation.')(func)
    func = click.option('--use_distance', is_flag=True, default=False, help='Use base pair distance in CGC annotation.')(func)

    # 扩展参数
    func = click.option('--extend_mode',
                        type=click.Choice(['none', 'bp', 'gene']),
                        default='none',
                        help="Extend CGC region on both sides after identification. 'bp' extends by base pairs; 'gene' extends by gene count; 'none' disables extension.")(func)
    func = click.option('--extend_bp',
                        type=int,
                        default=0,
                        help='When --extend_mode=bp, extend this many base pairs on each side.')(func)
    func = click.option('--extend_gene_count',
                        type=int,
                        default=0,
                        help='When --extend_mode=gene, extend this many genes on each side.')(func)

    # 新增阈值与 GFF 过滤
    func = click.option('--min_core_cazyme',
                        type=int, default=1,
                        help='Minimum number of core CAZymes required per CGC.')(func)
    func = click.option('--min_cluster_genes',
                        type=int, default=2,
                        help='Minimum number of genes required per CGC.')(func)
    func = click.option('--feature_type', 'feature_types',
                        multiple=True, default=["CDS"],
                        help='GFF feature types to include (multiple allowed).')(func)
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
    return func

def cgc_circle_plot_options(func):
    func = output_dir_option(func)
    return func

def topology_annotation_options(func):
    func = click.option('--run_signalp/--no-run_signalp',
                        default=False,
                        help='Run SignalP6.0 (biolib) to predict signal peptides for all proteins in overview')(func)
    func = click.option('--run_deeptmhmm/--no-run_deeptmhmm',
                        default=False,
                        help='Run DeepTMHMM (biolib) to predict transmembrane helices for all proteins in overview')(func)
    func = click.option('--signalp_org',
                        default='other',
                        type=click.Choice(['other', 'euk']),
                        show_default=True,
                        help='Organism type passed to SignalP6')(func)
    func = click.option('--force_topology/--no-force_topology',
                        default=False,
                        help='Overwrite existing SignalP / DeepTMHMM columns instead of only filling empty cells')(func)
    return func






