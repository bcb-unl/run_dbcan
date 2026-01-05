import logging
import subprocess
from pathlib import Path
import psutil
import pandas as pd
from Bio import SeqIO
import time

from dbcan.configs.signalp_tmhmm_config import SignalPTMHMMConfig
import dbcan.constants.signalp_tmhmm_constants as C

logger = logging.getLogger(__name__)

"""
Processor for running SignalP6 and annotating overview.tsv with results.
Currently only SignalP6 is implemented; DeepTMHMM was removed due to licensing issues.
"""

class SignalPTMHMMProcessor:
    def __init__(self, config: SignalPTMHMMConfig):
        self.config = config
        self._validate()

    @property
    def output_dir(self) -> Path:
        return Path(self.config.output_dir)

    @property
    def run_signalp(self) -> bool:
        return self.config.run_signalp

    @property
    def signalp_org(self) -> str:
        return self.config.signalp_org

    @property
    def signalp_mode(self) -> str:
        return self.config.signalp_mode

    @property
    def signalp_format(self) -> str:
        return self.config.signalp_format

    @property
    def force(self) -> bool:
        return self.config.force

    @property
    def threads(self) -> int:
        return self.config.threads or psutil.cpu_count()

    def _validate(self):
        if not self.output_dir.exists():
            raise FileNotFoundError(f"Output directory not found: {self.output_dir}")
        if self.signalp_org not in C.SIGNALP_ORGANISMS:
            raise ValueError(f"signalp_org must be one of: {C.SIGNALP_ORGANISMS}")
        if self.signalp_mode not in C.SIGNALP_MODES:
            raise ValueError(f"signalp_mode must be one of: {C.SIGNALP_MODES}")
        if self.signalp_format not in C.OUTPUT_FORMATS:
            raise ValueError(f"signalp_format must be one of: {C.OUTPUT_FORMATS}")
        if not self.run_signalp:
            logger.warning("No SignalP requested; nothing to run.")

    def load_overview(self) -> pd.DataFrame:
        path = self.output_dir / C.OVERVIEW_FILE
        if not path.exists():
            raise FileNotFoundError(f"overview file not found: {path}")
        # Use chunked reading for large files (>100MB)
        file_size_mb = path.stat().st_size / (1024 * 1024)
        if file_size_mb > 100:
            logger.info(f"Large overview file detected ({file_size_mb:.1f}MB), using chunked reading")
            chunks = []
            for chunk in pd.read_csv(path, sep="\t", chunksize=100000):
                chunks.append(chunk)
            return pd.concat(chunks, ignore_index=True)
        else:
            return pd.read_csv(path, sep="\t")

    def collect_all_gene_ids(self, df: pd.DataFrame) -> set:
        if C.GENE_ID_COL not in df.columns:
            raise ValueError(f"overview.tsv missing '{C.GENE_ID_COL}' column")
        return set(df[C.GENE_ID_COL].astype(str))

    def extract_fasta(self, gene_ids: set) -> Path:
        overview_path = self.output_dir / C.OVERVIEW_FILE
        if not overview_path.exists():
            raise FileNotFoundError(f"overview file not found: {overview_path}")
        try:
            # Use chunked reading for large files (>100MB)
            file_size_mb = overview_path.stat().st_size / (1024 * 1024)
            if file_size_mb > 100:
                logger.info(f"Large overview file detected ({file_size_mb:.1f}MB), using chunked reading")
                chunks = []
                for chunk in pd.read_csv(overview_path, sep="\t", chunksize=100000):
                    chunks.append(chunk)
                df = pd.concat(chunks, ignore_index=True)
            else:
                df = pd.read_csv(overview_path, sep="\t")
        except Exception as e:
            raise RuntimeError(f"failed to read overview.tsv: {e}")
        if C.GENE_ID_COL not in df.columns:
            raise ValueError(f"overview.tsv missing '{C.GENE_ID_COL}' column")

        if C.RECOMMEND_RESULTS_COL in df.columns:
            filtered = df[df[C.RECOMMEND_RESULTS_COL].astype(str) != C.DEFAULT_EMPTY]
            selected_ids = set(filtered[C.GENE_ID_COL].astype(str))
        else:
            selected_ids = set(df[C.GENE_ID_COL].astype(str))

        src = self.output_dir / C.INPUT_PROTEIN_NAME
        if not src.exists():
            raise FileNotFoundError(f"input protein FASTA not found: {src}")

        out_path = self.output_dir / C.TOPOLOGY_INPUT_FILE
        count = 0
        with open(src) as handle, open(out_path, "w") as out:
            for rec in SeqIO.parse(handle, "fasta"):
                rid = rec.id.split()[0]
                if rid in selected_ids:
                    rec.id = rid
                    rec.description = ""
                    SeqIO.write(rec, out, "fasta")
                    count += 1
        logger.info(f"Extracted {count} sequences for SignalP -> {out_path}")
        return out_path


    def run_signalp_predict(self, faa: Path) -> Path | None:
        out_dir = self.output_dir / C.SIGNALP_OUT_DIR
        out_dir.mkdir(exist_ok=True)
        cmd = [
            "signalp6",
            "--fastafile", str(faa),
            "--output_dir", str(out_dir),
            "--mode", self.signalp_mode,
            "--organism", self.signalp_org,
            "--format", self.signalp_format,
            "--torch_num_threads", str(self.threads)
        ]
        logger.info("Running SignalP6: %s", " ".join(cmd))
        start = time.time()
        try:
            result = subprocess.run(cmd, check=False, capture_output=True, text=True)
            elapsed = time.time() - start
            if result.returncode != 0:
                logger.error("SignalP6 failed (exit %s)", result.returncode)
                if result.stderr:
                    logger.error("SignalP6 stderr: %s", result.stderr.strip())
                raise subprocess.CalledProcessError(result.returncode, cmd, output=result.stdout, stderr=result.stderr)

            logger.info("SignalP6 finished in %.2fs: %s", elapsed, out_dir)
            return out_dir
        except FileNotFoundError:
            logger.error("signalp6 executable not found in PATH")
        except subprocess.CalledProcessError as e:
            logger.error("SignalP6 failed (exit %s)", e.returncode)
        except Exception as e:
            logger.error("SignalP6 unexpected error: %s", e)
        return None

    def parse_signalp_results(self, signalp_out_dir: Path) -> dict:
        results = {}
        result_file = signalp_out_dir / C.SIGNALP_RESULT_TSV
        if not result_file.exists():
            logger.warning(f"No SignalP result file found at {result_file}")
            return results
        with open(result_file) as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    results[parts[0]] = parts[1]
        return results

    def update_overview_with_signalp(self, signalp_results: dict):
        overview_path = self.output_dir / C.OVERVIEW_FILE
        try:
            # Use chunked reading for large files (>100MB)
            file_size_mb = overview_path.stat().st_size / (1024 * 1024)
            if file_size_mb > 100:
                logger.info(f"Large overview file detected ({file_size_mb:.1f}MB), using chunked reading")
                chunks = []
                for chunk in pd.read_csv(overview_path, sep="\t", chunksize=100000):
                    chunks.append(chunk)
                df = pd.concat(chunks, ignore_index=True)
            else:
                df = pd.read_csv(overview_path, sep="\t")
            if signalp_results:
                df[C.SIGNALP_COL] = df[C.GENE_ID_COL].astype(str).map(lambda x: signalp_results.get(x, C.DEFAULT_EMPTY))
            else:
                if C.SIGNALP_COL not in df.columns:
                    df[C.SIGNALP_COL] = C.DEFAULT_EMPTY
            df.to_csv(overview_path, sep="\t", index=False)
            logger.info("Updated overview with SignalP predictions")
        except Exception as e:
            logger.error(f"Error updating overview: {e}")

    def run(self):
        if not self.run_signalp:
            logger.info("SignalP not requested; skipping.")
            return {}
        try:
            df = self.load_overview()
        except Exception as e:
            logger.error("Failed to load overview: %s", e)
            return {}
        gene_ids = self.collect_all_gene_ids(df)
        if not gene_ids:
            logger.info("No gene IDs to process for SignalP")
            return {}
        fasta = self.extract_fasta(gene_ids)
        out_dir = self.run_signalp_predict(fasta)
        signalp_results = {}
        if out_dir:
            signalp_results = self.parse_signalp_results(out_dir)
        self.update_overview_with_signalp(signalp_results)
        return {"signalp_out": out_dir} if out_dir else {}
