import logging
import os
import shutil
import subprocess
import sys
from pathlib import Path
import psutil
import pandas as pd
from Bio import SeqIO
import time

from dbcan.configs.signalp_tmhmm_config import SignalPTMHMMConfig
import dbcan.constants.signalp_tmhmm_constants as C

logger = logging.getLogger(__name__)

"""
Processor for running optional local topology tools and annotating overview.tsv.
SignalP6 and DeepTMHMM must be installed separately by users.
"""


def _subprocess_env_with_python_lib(python_exe: str | None = None) -> dict:
    """Prepend ``<prefix>/lib`` to LD_LIBRARY_PATH for a Python executable.

    Conda-forge wheels (e.g. matplotlib) may require a newer ``libstdc++`` than
    the system copy under ``/lib/x86_64-linux-gnu``; putting the env's ``lib``
    first avoids ``CXXABI_* not found`` when the dynamic loader picks the wrong
    ``libstdc++.so.6``.
    """
    env = os.environ.copy()
    exe_path: Path | None = None
    if python_exe and python_exe.strip():
        p = Path(python_exe.strip())
        if p.is_file():
            exe_path = p.resolve()
        else:
            resolved = shutil.which(python_exe.strip())
            if resolved:
                exe_path = Path(resolved).resolve()
    if exe_path is None:
        exe_path = Path(sys.executable).resolve()

    lib_dir = exe_path.parent.parent / "lib"
    if not lib_dir.is_dir():
        return env

    prefix = str(lib_dir)
    prev = env.get("LD_LIBRARY_PATH", "")
    if prev:
        parts = prev.split(os.pathsep)
        if prefix not in parts:
            env["LD_LIBRARY_PATH"] = prefix + os.pathsep + prev
    else:
        env["LD_LIBRARY_PATH"] = prefix
    return env

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
    def run_deeptmhmm(self) -> bool:
        return self.config.run_deeptmhmm

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

    @property
    def deeptmhmm_dir(self) -> Path | None:
        if not self.config.deeptmhmm_dir:
            return None
        return Path(self.config.deeptmhmm_dir).expanduser().resolve()

    @property
    def deeptmhmm_python(self) -> str:
        return self.config.deeptmhmm_python or C.DEFAULT_DEEPTMHMM_PYTHON

    def _validate(self):
        if not self.output_dir.exists():
            raise FileNotFoundError(f"Output directory not found: {self.output_dir}")
        if self.signalp_org not in C.SIGNALP_ORGANISMS:
            raise ValueError(f"signalp_org must be one of: {C.SIGNALP_ORGANISMS}")
        if self.signalp_mode not in C.SIGNALP_MODES:
            raise ValueError(f"signalp_mode must be one of: {C.SIGNALP_MODES}")
        if self.signalp_format not in C.OUTPUT_FORMATS:
            raise ValueError(f"signalp_format must be one of: {C.OUTPUT_FORMATS}")
        if self.run_deeptmhmm:
            if self.deeptmhmm_dir is None:
                raise ValueError("--deeptmhmm_dir is required when --run_deeptmhmm is enabled")
            if not self.deeptmhmm_dir.exists():
                raise FileNotFoundError(f"DeepTMHMM directory not found: {self.deeptmhmm_dir}")
            predict_script = self.deeptmhmm_dir / C.DEEPTMHMM_PREDICT_SCRIPT
            if not predict_script.exists():
                raise FileNotFoundError(f"DeepTMHMM predict.py not found: {predict_script}")
        if not self.run_signalp and not self.run_deeptmhmm:
            logger.warning("No topology annotation tool requested; nothing to run.")

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
        logger.info(f"Extracted {count} sequences for topology annotation -> {out_path}")
        return out_path


    def run_signalp_predict(self, faa: Path) -> Path | None:
        out_dir = self.output_dir / C.SIGNALP_OUT_DIR
        out_dir.mkdir(exist_ok=True)
        faa_abs = faa.expanduser().resolve()
        out_dir_abs = out_dir.expanduser().resolve()
        cmd = [
            "signalp6",
            "--fastafile", str(faa_abs),
            "--output_dir", str(out_dir_abs),
            "--mode", self.signalp_mode,
            "--organism", self.signalp_org,
            "--format", self.signalp_format,
            "--torch_num_threads", str(self.threads)
        ]
        logger.info("Running SignalP6: %s", " ".join(cmd))
        start = time.time()
        try:
            result = subprocess.run(
                cmd,
                check=False,
                capture_output=True,
                text=True,
                env=_subprocess_env_with_python_lib(),
            )
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

    def _next_deeptmhmm_output_dir(self) -> Path:
        base = self.output_dir / C.DEEPTMHMM_OUT_DIR
        if not base.exists():
            return base
        idx = 1
        while True:
            candidate = self.output_dir / f"{C.DEEPTMHMM_OUT_DIR}_{idx}"
            if not candidate.exists():
                return candidate
            idx += 1

    def run_deeptmhmm_predict(self, faa: Path) -> Path | None:
        if self.deeptmhmm_dir is None:
            logger.error("DeepTMHMM directory is not configured")
            return None
        out_dir = self._next_deeptmhmm_output_dir()
        predict_script = self.deeptmhmm_dir / C.DEEPTMHMM_PREDICT_SCRIPT
        # Use absolute paths: subprocess cwd is deeptmhmm_dir, so relative paths would
        # resolve inside that directory and break (e.g. test_sp/... not found).
        faa_abs = faa.expanduser().resolve()
        out_dir_abs = out_dir.expanduser().resolve()
        cmd = [
            self.deeptmhmm_python,
            str(predict_script),
            "--fasta",
            str(faa_abs),
            "--output-dir",
            str(out_dir_abs),
        ]
        logger.info("Running DeepTMHMM: %s", " ".join(cmd))
        start = time.time()
        try:
            result = subprocess.run(
                cmd,
                check=False,
                capture_output=True,
                text=True,
                cwd=str(self.deeptmhmm_dir),
                env=_subprocess_env_with_python_lib(self.deeptmhmm_python),
            )
            elapsed = time.time() - start
            if result.returncode != 0:
                logger.error("DeepTMHMM failed (exit %s)", result.returncode)
                if result.stderr:
                    logger.error("DeepTMHMM stderr: %s", result.stderr.strip())
                raise subprocess.CalledProcessError(result.returncode, cmd, output=result.stdout, stderr=result.stderr)
            logger.info("DeepTMHMM finished in %.2fs: %s", elapsed, out_dir)
            return out_dir
        except FileNotFoundError:
            logger.error("DeepTMHMM Python executable not found: %s", self.deeptmhmm_python)
        except subprocess.CalledProcessError as e:
            logger.error("DeepTMHMM failed (exit %s)", e.returncode)
        except Exception as e:
            logger.error("DeepTMHMM unexpected error: %s", e)
        return None

    def parse_deeptmhmm_results(self, deeptmhmm_out_dir: Path) -> dict:
        results = {}
        result_file = deeptmhmm_out_dir / C.DEEPTMHMM_RESULT_3LINE
        if not result_file.exists():
            logger.warning(f"No DeepTMHMM result file found at {result_file}")
            return results

        with open(result_file) as f:
            lines = [line.strip() for line in f if line.strip()]

        for i in range(0, len(lines), 3):
            header = lines[i]
            if not header.startswith(">"):
                continue
            header = header[1:]
            if "|" in header:
                gene_id, prediction = header.split("|", 1)
                results[gene_id.strip()] = prediction.strip()
            else:
                results[header.split()[0]] = C.DEFAULT_EMPTY
        return results

    def update_overview_with_topology(self, signalp_results: dict | None = None, deeptmhmm_results: dict | None = None):
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
            elif self.run_signalp:
                if C.SIGNALP_COL not in df.columns:
                    df[C.SIGNALP_COL] = C.DEFAULT_EMPTY
            if deeptmhmm_results:
                df[C.DEEPTMHMM_COL] = df[C.GENE_ID_COL].astype(str).map(
                    lambda x: deeptmhmm_results.get(x, C.DEFAULT_EMPTY)
                )
            elif self.run_deeptmhmm:
                if C.DEEPTMHMM_COL not in df.columns:
                    df[C.DEEPTMHMM_COL] = C.DEFAULT_EMPTY
            df.to_csv(overview_path, sep="\t", index=False)
            logger.info("Updated overview with topology predictions")
        except Exception as e:
            logger.error(f"Error updating overview: {e}")

    def run(self):
        if not self.run_signalp and not self.run_deeptmhmm:
            logger.info("No topology annotation tool requested; skipping.")
            return {}
        try:
            df = self.load_overview()
        except Exception as e:
            logger.error("Failed to load overview: %s", e)
            return {}
        gene_ids = self.collect_all_gene_ids(df)
        if not gene_ids:
            logger.info("No gene IDs to process for topology annotation")
            return {}

        fasta = self.extract_fasta(gene_ids)
        signalp_out_dir = None
        deeptmhmm_out_dir = None
        signalp_results = {}
        deeptmhmm_results = {}

        if self.run_signalp:
            signalp_out_dir = self.run_signalp_predict(fasta)
        if signalp_out_dir:
            signalp_results = self.parse_signalp_results(signalp_out_dir)

        if self.run_deeptmhmm:
            deeptmhmm_out_dir = self.run_deeptmhmm_predict(fasta)
        if deeptmhmm_out_dir:
            deeptmhmm_results = self.parse_deeptmhmm_results(deeptmhmm_out_dir)

        self.update_overview_with_topology(signalp_results, deeptmhmm_results)

        results = {}
        if signalp_out_dir:
            results["signalp_out"] = signalp_out_dir
        if deeptmhmm_out_dir:
            results["deeptmhmm_out"] = deeptmhmm_out_dir
        return results