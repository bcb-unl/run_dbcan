import rich_click as click
from dbcan.parameter import (
    create_config, GeneralConfig, DBDownloaderConfig, DiamondConfig, DiamondTCConfig, DiamondPeptidaseConfig, DiamondSulfataseConfig,
    PyHMMERConfig, DBCANSUBProcessorConfig, PyHMMERTFConfig, PyHMMERSTPConfig,OverviewGeneratorConfig, GFFConfig, CGCFinderConfig, CGCSubstrateConfig,
    SynPlotConfig, CGCPlotConfig,

    general_options, database_options, output_dir_option, methods_option, threads_option, diamond_options, diamond_tc_options,
    pyhmmer_dbcan_options, dbcansub_options ,pyhmmer_tf, pyhmmer_stp, cgc_gff_option, cgc_options, cgc_sub_options, syn_plot_options,
    cgc_circle_plot_options, cgc_substrate_base_options, cgc_substrate_homology_params_options, cgc_substrate_dbcan_sub_param_options
)
from dbcan.core import (
    run_dbCAN_database, run_dbCAN_input_process, run_dbCAN_CAZyme_annotation,
    run_dbCAN_CGCFinder_preprocess, run_dbCAN_CGCFinder,
    run_dbCAN_CGCFinder_substrate, run_dbcan_syn_plot, run_dbCAN_cgc_circle
)

@click.group()
def cli():
    """use dbCAN tools to annotate and analyze CAZymes and CGCs."""
    pass

@cli.command('version')
@click.pass_context
def version_cmd(ctx):
    """show version information."""
    from dbcan._version import __version__
    click.echo(f"dbCAN version: {__version__}")

@cli.command('database')
@database_options
@click.pass_context
def database_cmd(ctx, **kwargs):
    """download dbCAN databases."""
    config = create_config(DBDownloaderConfig, **kwargs)
    run_dbCAN_database(config)

@cli.command('CAZyme_annotation')
@general_options
@database_options
@output_dir_option
@methods_option
@threads_option
@diamond_options
@pyhmmer_dbcan_options
@dbcansub_options
@click.pass_context
def cazyme_annotation_cmd(ctx, **kwargs):
    """annotate CAZyme using run_dbcan with prokaryotic, metagenomics, and protein sequences."""
    config = create_config(GeneralConfig, **kwargs)
    run_dbCAN_input_process(config)
    diamond_config = create_config(DiamondConfig,  **kwargs)
    pyhmmer_config = create_config(PyHMMERConfig, **kwargs)
    dbcansubconfig = create_config(DBCANSUBProcessorConfig, **kwargs)
    overviewconfig = create_config(OverviewGeneratorConfig, **kwargs)
    methods_option = kwargs.get('methods')
    run_dbCAN_CAZyme_annotation(diamond_config, pyhmmer_config, dbcansubconfig, overviewconfig, methods_option)




@cli.command('gff_process')
@database_options
@output_dir_option
@threads_option
@pyhmmer_stp
@pyhmmer_tf
@diamond_tc_options
@cgc_gff_option
@click.pass_context
def gff_process_cmd(ctx, **kwargs):
    """generate GFF file for CAZyme Gene Cluster identification."""
    diamond_tc_config = create_config(DiamondTCConfig, **kwargs)
    pyhmmer_tf_config = create_config(PyHMMERTFConfig, **kwargs)
    pyhmmer_stp_config = create_config(PyHMMERSTPConfig, **kwargs)
    diamond_sulfatlas_config = create_config(DiamondSulfataseConfig, **kwargs)
    diamond_peptidase_config = create_config(DiamondPeptidaseConfig, **kwargs)

    gff_config = create_config(GFFConfig, **kwargs)
    run_dbCAN_CGCFinder_preprocess(diamond_tc_config, pyhmmer_tf_config, pyhmmer_stp_config, diamond_sulfatlas_config, diamond_peptidase_config, gff_config)



@cli.command('cgc_finder')
@cgc_options
@click.pass_context
def cgc_finder_cmd(ctx, **kwargs):
    """identify CAZyme Gene Clusters(CGCs)"""
    config = create_config(CGCFinderConfig, **kwargs)
    run_dbCAN_CGCFinder(config)




@cli.command('substrate_prediction')
@syn_plot_options
@cgc_substrate_base_options
@cgc_substrate_homology_params_options
@cgc_substrate_dbcan_sub_param_options
@click.pass_context
def substrate_prediction_cmd(ctx, **kwargs):
    """predict substrate specificities of CAZyme Gene Clusters(CGCs)."""
    cgcsubconfig = create_config(CGCSubstrateConfig, **kwargs)
    run_dbCAN_CGCFinder_substrate(cgcsubconfig)
    synplotconfig = create_config(SynPlotConfig, **kwargs)
    run_dbcan_syn_plot(synplotconfig)



@cli.command('cgc_circle_plot')
@cgc_circle_plot_options
@click.pass_context
def cgc_circle_plot_cmd(ctx, **kwargs):
    """generate circular plots for CAZyme Gene Clusters(CGCs)."""
    config = create_config(CGCPlotConfig, **kwargs)
    run_dbCAN_cgc_circle(config)


@cli.command('easy_CGC')
@general_options
@database_options
@output_dir_option
@methods_option
@threads_option
@diamond_options
@pyhmmer_dbcan_options
@dbcansub_options
@pyhmmer_stp
@pyhmmer_tf
@diamond_tc_options
@cgc_gff_option
@cgc_options
@click.pass_context
def easy_cgc_cmd(ctx, **kwargs):
    """Perform complete CGC analysis: CAZyme annotation, GFF processing, and CGC identification in one step."""
    try:
        # step 1: CAZyme annotation
        click.echo("step 1/3  CAZyme annotation...")
        config = create_config(GeneralConfig, **kwargs)
        run_dbCAN_input_process(config)
        diamond_config = create_config(DiamondConfig, **kwargs)
        pyhmmer_config = create_config(PyHMMERConfig, **kwargs)
        dbcansubconfig = create_config(DBCANSUBProcessorConfig, **kwargs)
        overviewconfig = create_config(OverviewGeneratorConfig, **kwargs)
        methods_option = kwargs.get('methods')
        run_dbCAN_CAZyme_annotation(diamond_config, pyhmmer_config, dbcansubconfig, overviewconfig, methods_option)

        # step 2: GFF processing
        click.echo("step 2/3  GFF processing...")
        diamond_tc_config = create_config(DiamondTCConfig, **kwargs)
        pyhmmer_tf_config = create_config(PyHMMERTFConfig, **kwargs)
        pyhmmer_stp_config = create_config(PyHMMERSTPConfig, **kwargs)
        diamond_sulfatlas_config = create_config(DiamondSulfataseConfig, **kwargs)
        diamond_peptidase_config = create_config(DiamondPeptidaseConfig, **kwargs)
        gff_config = create_config(GFFConfig, **kwargs)
        run_dbCAN_CGCFinder_preprocess(diamond_tc_config, pyhmmer_tf_config, pyhmmer_stp_config, diamond_sulfatlas_config, diamond_peptidase_config, gff_config)
        # step 3: CGC identification
        click.echo("step 3/3  CGC identification...")
        cgc_config = create_config(CGCFinderConfig, **kwargs)
        run_dbCAN_CGCFinder(cgc_config)

        click.echo("CGC analysis completed.")
    except Exception as e:
        import traceback
        click.echo(f"error: {str(e)}")
        click.echo(traceback.format_exc())
        ctx.exit(1)




@cli.command('easy_substrate')
@general_options
@database_options
@output_dir_option
@methods_option
@threads_option
@diamond_options
@pyhmmer_dbcan_options
@dbcansub_options
@pyhmmer_stp
@pyhmmer_tf
@diamond_tc_options
@cgc_gff_option
@cgc_options
@syn_plot_options
@cgc_substrate_base_options
@cgc_substrate_homology_params_options
@cgc_substrate_dbcan_sub_param_options
@click.pass_context
def easy_substrate_cmd(ctx, **kwargs):
    """Perform complete CGC analysis: CAZyme annotation, GFF processing, CGC identification, and substrate prediction in one step."""
    try:
        # step 1: CAZyme annotation
        click.echo("step 1/4  CAZyme annotation...")
        config = create_config(GeneralConfig, **kwargs)
        run_dbCAN_input_process(config)
        diamond_config = create_config(DiamondConfig, **kwargs)
        pyhmmer_config = create_config(PyHMMERConfig, **kwargs)
        dbcansubconfig = create_config(DBCANSUBProcessorConfig, **kwargs)
        overviewconfig = create_config(OverviewGeneratorConfig, **kwargs)
        methods_option = kwargs.get('methods')
        run_dbCAN_CAZyme_annotation(diamond_config, pyhmmer_config, dbcansubconfig, overviewconfig, methods_option)

        # step 2: GFF processing
        click.echo("step 2/4  GFF processing...")
        diamond_tc_config = create_config(DiamondTCConfig, **kwargs)
        pyhmmer_tf_config = create_config(PyHMMERTFConfig, **kwargs)
        pyhmmer_stp_config = create_config(PyHMMERSTPConfig, **kwargs)
        diamond_sulfatlas_config = create_config(DiamondSulfataseConfig, **kwargs)
        diamond_peptidase_config = create_config(DiamondPeptidaseConfig, **kwargs)
        gff_config = create_config(GFFConfig, **kwargs)
        run_dbCAN_CGCFinder_preprocess(diamond_tc_config, pyhmmer_tf_config, pyhmmer_stp_config, diamond_sulfatlas_config, diamond_peptidase_config, gff_config)

        # step 3: CGC identification
        click.echo("step 3/4  CGC identification...")
        cgc_config = create_config(CGCFinderConfig, **kwargs)
        run_dbCAN_CGCFinder(cgc_config)

        # step 4: Substrate prediction
        click.echo("step 4/4  Substrate prediction...")
        cgcsubconfig = create_config(CGCSubstrateConfig, **kwargs)
        run_dbCAN_CGCFinder_substrate(cgcsubconfig)

        synplotconfig = create_config(SynPlotConfig, **kwargs)
        run_dbcan_syn_plot(synplotconfig)


    except Exception as e:
        import traceback
        click.echo(f"error: {str(e)}")
        click.echo(traceback.format_exc())
        ctx.exit(1)

    click.echo("CGC analysis completed.")

if __name__ == "__main__":
    cli()
