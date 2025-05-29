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
    from dbcan.annotation.diamond import CAZyDiamondProcessor
    processor = CAZyDiamondProcessor(config)
    processor.run()
    processor.format_results()

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
    generator.generate_non_cazyme_faa()


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
    processor.format_results()

def run_dbCAN_sulfatlas_diamond(config):
    from dbcan.annotation.diamond import SulfatlasDiamondProcessor
    processor = SulfatlasDiamondProcessor(config)
    processor.run()
    processor.format_results()

def run_dbCAN_peptidase_diamond(config):
    from dbcan.annotation.diamond import PeptidaseDiamondProcessor
    processor = PeptidaseDiamondProcessor(config)
    processor.run()
    processor.format_results()

def run_dbCAN_hmmer_tf(config):
    from dbcan.annotation.pyhmmer_search import PyHMMERTFProcessor
    processor = PyHMMERTFProcessor(config)
    processor.run()

def run_dbCAN_hmmer_stp(config):
    from dbcan.annotation.pyhmmer_search import PyHMMERSTPProcessor
    processor = PyHMMERSTPProcessor(config)
    processor.run()

def run_dbCAN_CGCFinder_preprocess(tcdbconfig, tfconfig, stpconfig, sulfatlasconfig,peptidaseconfig,cgcgffconfig):
    run_dbCAN_tcdb_diamond(tcdbconfig)
    run_dbCAN_hmmer_tf(tfconfig)
    run_dbCAN_hmmer_stp(stpconfig)
    run_dbCAN_sulfatlas_diamond(sulfatlasconfig)
    run_dbCAN_peptidase_diamond(peptidaseconfig)


    from dbcan.process.process_utils import process_cgc_sig_results
    process_cgc_sig_results(tcdbconfig, tfconfig, stpconfig, sulfatlasconfig, peptidaseconfig)
    from dbcan.IO.gff import get_gff_processor
    processor = get_gff_processor(cgcgffconfig)
    processor.process_gff()

def run_dbCAN_CGCFinder(config):
    from dbcan.annotation.CGCFinder import CGCFinder
    cgc_finder = CGCFinder(config)
    cgc_finder.run()



def run_dbCAN_CGCFinder_substrate(config):
    from dbcan.annotation.cgc_substrate_prediction import cgc_substrate_prediction
    cgc_substrate_prediction(config)



def run_dbcan_syn_plot(config):
    from dbcan.plot.syntenic_plot import SyntenicPlot
    syntenic_plot = SyntenicPlot(config)

    syntenic_plot.syntenic_plot_allpairs(config)

def run_dbCAN_cgc_circle(config):
    from dbcan.plot.plot_cgc_circle import CGCCircosPlot
    cgc_plot = CGCCircosPlot(config)
    cgc_plot.plot()


