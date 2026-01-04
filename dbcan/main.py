import rich_click as click
from dbcan.parameter import (
    create_config, general_options, database_options, output_dir_option, methods_option, threads_option, diamond_options, diamond_tc_options,diamond_tf_options,
    pyhmmer_dbcan_options, dbcansub_options ,pyhmmer_tf, pyhmmer_stp, cgc_gff_option, cgc_options, cgc_sub_options, syn_plot_options,
    cgc_circle_plot_options, cgc_substrate_base_options, cgc_substrate_homology_params_options, cgc_substrate_dbcan_sub_param_options, pyhmmer_pfam,
    topology_annotation_options   # <--- added
    , diamond_sulfatase_options, diamond_peptidase_options
    , database_download_options, logging_options
)
from pathlib import Path
import logging
import warnings
from typing import Optional

# Suppress numpy UserWarnings about subnormal values
warnings.filterwarnings("ignore", message=".*smallest subnormal.*", category=UserWarning)


def setup_logging(log_level: str, log_file: Optional[str] = None, verbose: bool = False):
    """Configure logging based on command line options"""
    level = logging.DEBUG if verbose else getattr(logging, log_level.upper(), logging.WARNING)
    
    handlers = [logging.StreamHandler()]
    if log_file:
        handlers.append(logging.FileHandler(log_file, mode='w'))
    
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=handlers,
        force=True  # Override any existing configuration
    )

def _invoke_subset(ctx, cmd, all_kwargs):
    """
    Invoke a click command with a subset of parameters from all_kwargs.
    This allows reusing commands with only the relevant parameters.
    """
    names = {p.name for p in cmd.params}
    sub_kwargs = {k: v for k, v in all_kwargs.items() if k in names}
    return ctx.invoke(cmd, **sub_kwargs)

@click.group()
@logging_options
@click.pass_context
def cli(ctx, log_level, log_file, verbose):
    """use dbCAN tools to annotate and analyze CAZymes and CGCs."""
    setup_logging(log_level, log_file, verbose)
    ctx.ensure_object(dict)
    ctx.obj['log_level'] = log_level
    ctx.obj['log_file'] = log_file
    ctx.obj['verbose'] = verbose

@cli.command('version')
@click.pass_context
def version_cmd(ctx):
    """show version information."""
    from dbcan._version import __version__
    click.echo(f"dbCAN version: {__version__}")

@cli.command('database')
@logging_options
@database_download_options
@click.pass_context
def database_cmd(ctx, log_level, log_file, verbose, **kwargs):
    """download dbCAN databases."""
    setup_logging(log_level, log_file, verbose)
    from dbcan.configs.database_config import DBDownloaderConfig
    from dbcan.core import run_dbCAN_database
    config = create_config(DBDownloaderConfig, **kwargs)
    run_dbCAN_database(config)

@cli.command('CAZyme_annotation')
@logging_options
@general_options
@database_options
@methods_option
@threads_option
@diamond_options
@pyhmmer_dbcan_options
@dbcansub_options
@topology_annotation_options
@click.pass_context
def cazyme_annotation_cmd(ctx, log_level, log_file, verbose, **kwargs):
    """annotate CAZyme using run_dbcan with prokaryotic, metagenomics, and protein sequences."""
    setup_logging(log_level, log_file, verbose)
    from dbcan.configs.signalp_tmhmm_config import SignalPTMHMMConfig
    from dbcan.configs.diamond_config import DiamondCAZyConfig
    from dbcan.configs.pyhmmer_config import PyHMMERDBCANConfig
    from dbcan.configs.pyhmmer_config import DBCANSUBConfig
    from dbcan.configs.base_config import OverviewGeneratorConfig
    from dbcan.configs.base_config import GeneralConfig
    from dbcan.core import run_dbCAN_input_process, run_dbCAN_CAZyme_annotation, run_dbCAN_topology_annotation

    # Step 1: Input processing
    config = create_config(GeneralConfig, **kwargs)
    run_dbCAN_input_process(config)

    # Step 2: CAZyme annotation
    diamond_config = create_config(DiamondCAZyConfig, **kwargs)
    pyhmmer_config = create_config(PyHMMERDBCANConfig, namespace="dbcan", **kwargs)
    dbcansubconfig = create_config(DBCANSUBConfig, namespace="dbsub", **kwargs)
    overviewconfig = create_config(OverviewGeneratorConfig, **kwargs)
    methods_option = kwargs.get('methods')
    run_dbCAN_CAZyme_annotation(diamond_config, pyhmmer_config, dbcansubconfig, overviewconfig, methods_option)

    # Step 3: Topology annotation (as supplement to CAZyme annotation)
    if kwargs.get('run_signalp', False):
        signalp_config = create_config(SignalPTMHMMConfig, **kwargs)
        run_dbCAN_topology_annotation(signalp_config)

@cli.command('gff_process')
@logging_options
@database_options
@output_dir_option
@threads_option
@pyhmmer_stp
@pyhmmer_tf
@diamond_tf_options
@diamond_tc_options
@diamond_sulfatase_options
@diamond_peptidase_options
@cgc_gff_option
@click.pass_context
def gff_process_cmd(ctx, log_level, log_file, verbose, **kwargs):
    """
    Generate GFF for CGC identification.
    need --input_gff when --input_raw_data is protein sequence.
    if --input_gff is not provided, will set default <output_dir>/uniInput.gff.

    """
    setup_logging(log_level, log_file, verbose)
    from dbcan.configs.diamond_config import DiamondTFConfig, DiamondTCConfig, DiamondSulfataseConfig, DiamondPeptidaseConfig
    from dbcan.configs.pyhmmer_config import PyHMMERTFConfig, PyHMMERSTPConfig
    from dbcan.configs.base_config import GFFConfig
    from dbcan.core import run_dbCAN_CGCFinder_preprocess

    mode = kwargs.get('mode', 'prok')
    if not kwargs.get('input_gff'):
        if mode != 'protein':
            auto_gff = Path(kwargs['output_dir']).joinpath('uniInput.gff')
            kwargs['input_gff'] = str(auto_gff)
            if not kwargs.get('gff_type'):
                kwargs['gff_type'] = 'prodigal'
            logging.info(f"[gff_process] Auto-set input_gff={kwargs['input_gff']} gff_type={kwargs['gff_type']} (mode={mode})")
        else:
            raise click.UsageError("--input_gff is required when --mode=protein")

    if not kwargs.get('gff_type'):
        # Fallback if user gave a custom GFF in non-protein mode but forgot type
        kwargs['gff_type'] = 'prodigal'
        logging.info(f"[gff_process] Defaulting gff_type=prodigal")

    diamond_tc_config = create_config(DiamondTCConfig, namespace="tc", **kwargs)
    diamond_tf_config = create_config(DiamondTFConfig, namespace="tf_diamond", **kwargs)
    pyhmmer_tf_config = create_config(PyHMMERTFConfig, namespace="tf", **kwargs)
    pyhmmer_stp_config = create_config(PyHMMERSTPConfig, namespace="stp", **kwargs)
    diamond_sulfatlas_config = create_config(DiamondSulfataseConfig, namespace="sulfatase", **kwargs)
    diamond_peptidase_config = create_config(DiamondPeptidaseConfig, namespace="peptidase", **kwargs)

    gff_config = create_config(GFFConfig, **kwargs)
    run_dbCAN_CGCFinder_preprocess(
        diamond_tc_config,
        diamond_tf_config,
        pyhmmer_tf_config,
        pyhmmer_stp_config,
        diamond_sulfatlas_config,
        diamond_peptidase_config,
        gff_config
    )

@cli.command('cgc_finder')
@logging_options
@output_dir_option
@cgc_options
@click.pass_context
def cgc_finder_cmd(ctx, log_level, log_file, verbose, **kwargs):
    """identify CAZyme Gene Clusters(CGCs)"""
    setup_logging(log_level, log_file, verbose)
    from dbcan.configs.cgcfinder_config import CGCFinderConfig
    from dbcan.core import run_dbCAN_CGCFinder
    config = create_config(CGCFinderConfig, **kwargs)
    run_dbCAN_CGCFinder(config)


@cli.command('Pfam_null_cgc')
@logging_options
@threads_option
@database_options
@output_dir_option
@pyhmmer_pfam
@click.pass_context
def pfam_null_cgc_cmd(ctx, log_level, log_file, verbose, **kwargs):
    """identify CAZyme Gene Clusters(CGCs)"""
    setup_logging(log_level, log_file, verbose)
    from dbcan.configs.pyhmmer_config import PyHMMERPfamConfig
    from dbcan.core import run_dbCAN_Pfam_null_cgc
    config = create_config(PyHMMERPfamConfig, namespace="pfam", **kwargs)
    run_dbCAN_Pfam_null_cgc(config)


@cli.command('substrate_prediction')
@logging_options
@cgc_substrate_base_options
@cgc_substrate_homology_params_options
@cgc_substrate_dbcan_sub_param_options
@click.pass_context
def substrate_prediction_cmd(ctx, log_level, log_file, verbose, **kwargs):
    setup_logging(log_level, log_file, verbose)
    from dbcan.configs.cgc_substrate_config import CGCSubstrateConfig, SynPlotConfig
    from dbcan.core import run_dbCAN_CGCFinder_substrate, run_dbcan_syn_plot
    """predict substrate specificities of CAZyme Gene Clusters(CGCs)."""
    cgcsubconfig = create_config(CGCSubstrateConfig, **kwargs)
    run_dbCAN_CGCFinder_substrate(cgcsubconfig)
    synplotconfig = create_config(SynPlotConfig, **kwargs)
    run_dbcan_syn_plot(synplotconfig)



@cli.command('cgc_circle_plot')
@logging_options
@cgc_circle_plot_options
@click.pass_context
def cgc_circle_plot_cmd(ctx, log_level, log_file, verbose, **kwargs):
    """generate circular plots for CAZyme Gene Clusters(CGCs)."""
    setup_logging(log_level, log_file, verbose)
    from dbcan.configs.base_config import CGCPlotConfig
    from dbcan.core import run_dbCAN_cgc_circle
    config = create_config(CGCPlotConfig, **kwargs)
    run_dbCAN_cgc_circle(config)


@cli.command('easy_CGC')
@logging_options
@general_options
@database_options
@methods_option
@threads_option
@diamond_options
@pyhmmer_dbcan_options
@dbcansub_options
@pyhmmer_stp
@pyhmmer_tf
@diamond_tf_options
@diamond_tc_options
@cgc_gff_option
@cgc_options
@click.pass_context
def easy_cgc_cmd(ctx, log_level, log_file, verbose, **kwargs):
    """Perform complete CGC analysis: CAZyme annotation, GFF processing, and CGC identification in one step."""
    setup_logging(log_level, log_file, verbose)
    try:
        # step 1: CAZyme annotation
        click.echo("step 1/3  CAZyme annotation...")
        _invoke_subset(ctx, cazyme_annotation_cmd, kwargs)
        # step 2: GFF processing
        click.echo("step 2/3  GFF processing...")
        _invoke_subset(ctx, gff_process_cmd, kwargs)
        # step 3: CGC identification
        click.echo("step 3/3  CGC identification...")
        _invoke_subset(ctx, cgc_finder_cmd, kwargs)
        click.echo("CGC analysis completed.")
    except Exception as e:
        import traceback
        click.echo(f"error: {str(e)}")
        click.echo(traceback.format_exc())
        ctx.exit(1)


@cli.command('easy_substrate')
@logging_options
@cgc_sub_options
@methods_option
@threads_option
@diamond_options
@pyhmmer_dbcan_options
@dbcansub_options
@pyhmmer_stp
@pyhmmer_tf
@diamond_tf_options
@diamond_tc_options
@cgc_gff_option
@cgc_options
@click.pass_context
def easy_substrate_cmd(ctx, log_level, log_file, verbose, **kwargs):
    """Perform complete CGC analysis: CAZyme annotation, GFF processing, CGC identification, and substrate prediction in one step."""
    setup_logging(log_level, log_file, verbose)
    try:
        # step 1: CAZyme annotation
        click.echo("step 1/4  CAZyme annotation...")
        _invoke_subset(ctx, cazyme_annotation_cmd, kwargs)

        # step 2: GFF processing
        click.echo("step 2/4  GFF processing...")
        _invoke_subset(ctx, gff_process_cmd, kwargs)

        # step 3: CGC identification
        click.echo("step 3/4  CGC identification...")
        _invoke_subset(ctx, cgc_finder_cmd, kwargs)

        # step 4: Substrate prediction
        click.echo("step 4/4  Substrate prediction...")
        _invoke_subset(ctx, substrate_prediction_cmd, kwargs)


    except Exception as e:
        import traceback
        click.echo(f"error: {str(e)}")
        click.echo(traceback.format_exc())
        ctx.exit(1)

    click.echo("CGC substrate analysis completed.")

if __name__ == "__main__":
    cli()
