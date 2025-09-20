import os

def run_dbCAN_database(config):
    from dbcan.utils.database import DBDownloader
    downloader = DBDownloader(config)
    downloader.download_file()

def run_dbCAN_input_process(config):
    from dbcan.IO.fasta import get_processor
    processor = get_processor(config)
    processor.process_input()

def run_dbCAN_cazy_diamond(config):
    from dbcan.annotation.diamond import CAZYDiamondProcessor
    processor = CAZYDiamondProcessor(config)
    processor.run()
    #processor.format_results()

def run_dbCAN_hmmer(config):
    from dbcan.annotation.pyhmmer_search import PyHMMERDBCANProcessor
    processor = PyHMMERDBCANProcessor(config)
    processor.run()

def run_dbCAN_dbcansub_hmmer(config):
    from dbcan.annotation.pyhmmer_search import PyHMMERDBCANSUBProcessor
    processor = PyHMMERDBCANSUBProcessor(config)
    processor.run()


def run_dbCAN_CAZyme_overview(config):
    from dbcan.IO.OverviewGenerator import OverviewGenerator
    generator = OverviewGenerator(config)
    generator.run()


def run_dbCAN_CAZyme_annotation(diamondconfig, dbcanconfig, dbcansubconfig, overviewconfig, methods):
    import logging
    if 'diamond' in methods:
        logging.info("DIAMOND CAZy...")
        try:
            run_dbCAN_cazy_diamond(diamondconfig)
        except Exception as e:
            logging.error(f"DIAMOND CAZy failed: {e}")

    if 'hmm' in methods:
        logging.info("pyhmmer vs dbCAN-HMM...")
        try:
            run_dbCAN_hmmer(dbcanconfig)
            logging.info("HMMER dbCAN done")
        except Exception as e:
            logging.error(f"HMMER dbCAN failed: {e}")

    if 'dbCANsub' in methods:
        logging.info("pyhmmer vs dbCAN-sub-HMM...")
        try:
            run_dbCAN_dbcansub_hmmer(dbcansubconfig)
            logging.info("dbCAN-sub HMM done")
        except Exception as e:
            logging.error(f"dbCAN-sub HMM failed: {e}")

    logging.info("generate overview of CAZymes...")
    try:
        run_dbCAN_CAZyme_overview(overviewconfig)
        logging.info("CAZyme overview generated")
    except Exception as e:
        logging.error(f"CAZyme overview failed: {e}")
#    else:
#        logging.warning("No CAZyme results to generate overview.")


def run_dbCAN_tcdb_diamond(config):
    from dbcan.annotation.diamond import TCDBDiamondProcessor
    processor = TCDBDiamondProcessor(config)
    processor.run()
    #processor.format_results()

def run_dbCAN_sulfatlas_diamond(config):
    from dbcan.annotation.diamond import SulfatlasDiamondProcessor
    processor = SulfatlasDiamondProcessor(config)
    processor.run()
    #processor.format_results()

def run_dbCAN_peptidase_diamond(config):
    from dbcan.annotation.diamond import PeptidaseDiamondProcessor
    processor = PeptidaseDiamondProcessor(config)
    processor.run()
    #processor.format_results()

def run_dbCAN_diamond_tf(config):
    from dbcan.annotation.diamond import TFDiamondProcessor
    processor = TFDiamondProcessor(config)
    processor.run()
    #processor.format_results()

def run_dbCAN_hmmer_tf(config):
    from dbcan.annotation.pyhmmer_search import PyHMMERTFProcessor
    processor = PyHMMERTFProcessor(config)
    processor.run()

def run_dbCAN_hmmer_stp(config):
    from dbcan.annotation.pyhmmer_search import PyHMMERSTPProcessor
    processor = PyHMMERSTPProcessor(config)
    processor.run()

def run_dbCAN_CGCFinder_preprocess(tcdbconfig, tfdiamondconfig, tfconfig, stpconfig, sulfatlasconfig, peptidaseconfig, cgcgffconfig):
    run_dbCAN_tcdb_diamond(tcdbconfig)
    if getattr(tfdiamondconfig, 'prokaryotic', True):
        run_dbCAN_diamond_tf(tfdiamondconfig)
    if getattr(tfconfig, 'fungi', False):
        run_dbCAN_hmmer_tf(tfconfig)
    run_dbCAN_hmmer_stp(stpconfig)
    run_dbCAN_sulfatlas_diamond(sulfatlasconfig)
    run_dbCAN_peptidase_diamond(peptidaseconfig)


    from dbcan.process.process_utils import process_cgc_sig_results
    process_cgc_sig_results(
        tcdbconfig,
        tfdiamondconfig if getattr(tfdiamondconfig, 'prokaryotic', True) else None,
        tfconfig if getattr(tfconfig, 'fungi', False) else None,
        stpconfig,
        sulfatlasconfig,
        peptidaseconfig
    )
    from dbcan.IO.gff import get_gff_processor
    processor = get_gff_processor(cgcgffconfig)
    processor.process_gff()

def run_dbCAN_CGCFinder(config):
    from dbcan.annotation.CGCFinder import CGCFinder
    cgc_finder = CGCFinder(config)
    cgc_finder.run()

def run_dbCAN_Pfam_null_cgc(config):
    from dbcan.process.process_utils import (
        process_cgc_null_pfam_annotation,
        extract_null_fasta_from_cgc,
        annotate_cgc_null_with_pfam_and_gff,
        extract_null_fasta_from_gff
    )
    from dbcan.annotation.pyhmmer_search import PyHMMERPfamProcessor

    # choose the source of null genes
    if getattr(config, 'null_from_gff', False):
        extract_null_fasta_from_gff(
            os.path.join(config.output_dir, 'cgc.gff'),
            os.path.join(config.output_dir, 'uniInput.faa'),
            os.path.join(config.output_dir, 'null_proteins.faa')
        )
    else:
        extract_null_fasta_from_cgc(
            os.path.join(config.output_dir, 'cgc_standard_out.tsv'),
            os.path.join(config.output_dir, 'uniInput.faa'),
            os.path.join(config.output_dir, 'null_proteins.faa')
        )

    pfam_processor = PyHMMERPfamProcessor(config)
    pfam_processor.run()
    process_cgc_null_pfam_annotation(config)
    annotate_cgc_null_with_pfam_and_gff(
        os.path.join(config.output_dir, 'cgc_standard_out.tsv'),
        os.path.join(config.output_dir, 'Pfam_hmm_results.tsv'),
        os.path.join(config.output_dir, 'cgc.gff'),
        os.path.join(config.output_dir, 'cgc_standard_out.pfam_annotated.tsv'),
        os.path.join(config.output_dir, 'cgc.pfam_annotated.gff')
    )

def run_dbCAN_CGCFinder_substrate(config):
    from dbcan.annotation.cgc_substrate_prediction import cgc_substrate_prediction
    cgc_substrate_prediction(config)



def run_dbcan_syn_plot(config):
    from dbcan.plot.syntenic_plot import SyntenicPlot
    syntenic_plot = SyntenicPlot(config)

    syntenic_plot.syntenic_plot_allpairs()

def run_dbCAN_cgc_circle(config):
    from dbcan.plot.plot_cgc_circle import CGCCircosPlot
    cgc_plot = CGCCircosPlot(config)
    cgc_plot.plot()

def run_dbCAN_topology_annotation(config):
    """
    Run SignalP6 to annotate proteins in overview.tsv with signal peptide information.
    DeepTMHMM has been removed due to licensing issues.
    """
    import logging
    if not config.run_signalp:
        logging.info("No SignalP requested; skipping.")
        return
    try:
        from dbcan.annotation.signalp_tmhmm import SignalPTMHMMProcessor
        processor = SignalPTMHMMProcessor(config)
        results = processor.run()
        if results and 'signalp_out' in results:
            logging.info(f"SignalP results: {results['signalp_out']}")
        else:
            logging.warning("SignalP produced no results")
    except ImportError as e:
        logging.error(f"SignalP module import failed: {e}")
    except Exception as e:
        logging.error(f"SignalP annotation failed: {e}")
        import traceback
        traceback.print_exc()


