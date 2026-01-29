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

def _ensure_empty_diamond_file(config):
    """Ensure empty diamond results file exists even if not run."""
    import logging
    from pathlib import Path
    import pandas as pd
    import dbcan.constants.diamond_constants as D
    
    output_path = Path(config.output_dir) / config.output_file
    if not output_path.exists():
        out_cols = getattr(config, "column_names", None) or D.CAZY_COLUMN_NAMES
        if out_cols:
            try:
                empty_df = pd.DataFrame(columns=out_cols)
                output_path.parent.mkdir(parents=True, exist_ok=True)
                empty_df.to_csv(output_path, sep='\t', index=False)
                logging.info(f"Created empty diamond results file (not in methods) -> {output_path.name}")
            except Exception as e:
                logging.error(f"Failed to create empty diamond file: {e}", exc_info=True)

def _ensure_empty_dbcan_hmm_file(config):
    """Ensure empty dbCAN HMM results file exists even if not run."""
    import logging
    from pathlib import Path
    import pandas as pd
    import dbcan.constants.process_utils_constants as P
    
    output_path = Path(config.output_dir) / config.output_file
    if not output_path.exists():
        out_cols = P.HMMER_COLUMN_NAMES
        if out_cols:
            try:
                empty_df = pd.DataFrame(columns=out_cols)
                output_path.parent.mkdir(parents=True, exist_ok=True)
                empty_df.to_csv(output_path, sep='\t', index=False)
                logging.info(f"Created empty dbCAN HMM results file (not in methods) -> {output_path.name}")
            except Exception as e:
                logging.error(f"Failed to create empty dbCAN HMM file: {e}", exc_info=True)

def _ensure_empty_dbcansub_file(config):
    """Ensure empty dbCAN-sub results file and raw file exist even if not run."""
    import logging
    from pathlib import Path
    import pandas as pd
    import dbcan.constants.process_dbcan_sub_constants as P
    import dbcan.constants.process_utils_constants as P_UTILS
    
    # Create empty final results file
    output_path = Path(config.output_dir) / P.DBCAN_SUB_HMM_RESULT_FILE
    if not output_path.exists():
        out_cols = getattr(P, "DBCAN_SUB_COLUMN_NAMES", None)
        if out_cols:
            try:
                empty_df = pd.DataFrame(columns=out_cols)
                output_path.parent.mkdir(parents=True, exist_ok=True)
                empty_df.to_csv(output_path, sep='\t', index=False)
                logging.info(f"Created empty dbCAN-sub results file (not in methods) -> {output_path.name}")
            except Exception as e:
                logging.error(f"Failed to create empty dbCAN-sub results file: {e}", exc_info=True)
    
    # Create empty raw file (dbCANsub_hmm_raw.tsv.raw.tsv)
    raw_output_path = Path(config.output_dir) / config.output_file
    raw_output_path = raw_output_path.with_suffix(raw_output_path.suffix + ".raw.tsv")
    if not raw_output_path.exists():
        raw_cols = P_UTILS.HMMER_COLUMN_NAMES
        if raw_cols:
            try:
                empty_df = pd.DataFrame(columns=raw_cols)
                raw_output_path.parent.mkdir(parents=True, exist_ok=True)
                empty_df.to_csv(raw_output_path, sep='\t', index=False)
                logging.info(f"Created empty dbCAN-sub raw file (not in methods) -> {raw_output_path.name}")
            except Exception as e:
                logging.error(f"Failed to create empty dbCAN-sub raw file: {e}", exc_info=True)

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
            logging.info("DIAMOND CAZy done")
        except Exception as e:
            logging.error(f"DIAMOND CAZy failed: {e}")
            # Even if diamond fails, try to create empty result file
            _ensure_empty_diamond_file(diamondconfig)
    else:
        # Even if diamond is not in methods, create empty result file to avoid warnings
        _ensure_empty_diamond_file(diamondconfig)

    if 'hmm' in methods:
        logging.info("pyhmmer vs dbCAN-HMM...")
        try:
            run_dbCAN_hmmer(dbcanconfig)
            logging.info("HMMER dbCAN done")
        except Exception as e:
            logging.error(f"HMMER dbCAN failed: {e}")
            # Even if hmm fails, try to create empty result file
            _ensure_empty_dbcan_hmm_file(dbcanconfig)
    else:
        # Even if hmm is not in methods, create empty result file to avoid warnings
        _ensure_empty_dbcan_hmm_file(dbcanconfig)

    if 'dbCANsub' in methods:
        logging.info("pyhmmer vs dbCAN-sub-HMM...")
        try:
            run_dbCAN_dbcansub_hmmer(dbcansubconfig)
            logging.info("dbCAN-sub HMM done")
        except Exception as e:
            logging.error(f"dbCAN-sub HMM failed: {e}", exc_info=True)
            # Even if HMM search fails, we should still try to create empty result file
            # This is handled inside PyHMMERDBCANSUBProcessor.run()
            # But we also need to ensure raw file exists
            _ensure_empty_dbcansub_file(dbcansubconfig)
    else:
        # Even if dbCANsub is not in methods, create empty result file to avoid warnings
        _ensure_empty_dbcansub_file(dbcansubconfig)

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


