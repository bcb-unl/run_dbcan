"""
Structure-based CAZyme annotation using Foldseek against CAZyme3D database.

Searches user protein sequences (FASTA) against a pre-built Foldseek database
(e.g. CAZyme3D), then maps target IDs to CAZy families via a mapping table.
Output format is compatible with Overview (Gene ID, CAZy ID).
"""
import logging
import shutil
import subprocess
import time
from pathlib import Path

import pandas as pd

from dbcan.configs.structure_search_config import StructureSearchConfig
import dbcan.constants.structure_search_constants as S

logger = logging.getLogger(__name__)

# Foldseek m8 format columns (default easy-search output)
FOLDSEEK_M8_COLS = [
    "query", "target", "ident", "alnlen", "mismatch", "gapopen",
    "qstart", "qend", "tstart", "tend", "evalue", "bits",
]


def _find_foldseek() -> str:
    """Return path to foldseek binary or raise if not found."""
    exe = shutil.which(S.FOLDSEEK_CMD)
    if not exe:
        raise FileNotFoundError(
            f"'{S.FOLDSEEK_CMD}' not found in PATH. "
            "Install Foldseek (https://github.com/steineggerlab/foldseek) to use structure search."
        )
    return exe


class StructureSearchProcessor:
    """Run Foldseek search against CAZyme3D (or compatible) structure DB and output CAZy IDs."""

    def __init__(self, config: StructureSearchConfig):
        self.config = config
        self._validate()

    @property
    def db_path(self) -> Path:
        """Path to Foldseek database (prefix; Foldseek uses <prefix>_* files)."""
        return Path(self.config.db_dir) / self.config.db_file

    @property
    def mapping_path(self) -> Path:
        """Path to target_id -> CAZy family mapping TSV."""
        return Path(self.config.db_dir) / self.config.mapping_file

    @property
    def input_faa(self) -> Path:
        return Path(self.config.output_dir) / self.config.input_faa

    @property
    def output_file(self) -> Path:
        return Path(self.config.output_dir) / self.config.output_file

    def _validate(self) -> None:
        if not self.input_faa.exists():
            raise FileNotFoundError(f"Input protein file not found: {self.input_faa}")
        # Foldseek DB can be a directory or a prefix; allow _a3m.ffindex or the prefix itself
        db_dir = Path(self.config.db_dir)
        if not db_dir.exists():
            raise FileNotFoundError(f"Database directory not found: {self.config.db_dir}")
        # Foldseek DB is a prefix; createdb produces <prefix>, <prefix>_h, <prefix>_a3m.ffindex, etc.
        prefix_path = db_dir / self.config.db_file
        if not prefix_path.exists() and not (db_dir / f"{self.config.db_file}_a3m.ffindex").exists():
            raise FileNotFoundError(
                f"Foldseek database not found under {self.config.db_dir} with prefix '{self.config.db_file}'. "
                "Provide a pre-built CAZyme3D Foldseek DB or build one with 'foldseek createdb'."
            )
        if not self.mapping_path.exists():
            raise FileNotFoundError(
                f"CAZyme3D-to-CAZy mapping file not found: {self.mapping_path}. "
                "Required columns: target_id, CAZy ID (or cazy_family)."
            )
        self.output_file.parent.mkdir(parents=True, exist_ok=True)

    def _load_mapping(self) -> pd.DataFrame:
        """Load target_id -> CAZy ID mapping. Expects TSV with at least target id and CAZy column."""
        df = pd.read_csv(self.mapping_path, sep="\t")
        # Accept common column names
        id_col = None
        cazy_col = None
        for c in df.columns:
            c_l = c.strip().lower()
            if c_l in ("target", "target_id", "targetid", "id", "accession", "cazyid", "uniprot_mapped"):
                id_col = c
            if c_l in ("cazy", "cazy_id", "cazy_family", "family", "cazy id"):
                cazy_col = c
        if id_col is None:
            id_col = df.columns[0]
        if cazy_col is None:
            cazy_col = df.columns[1] if len(df.columns) > 1 else df.columns[0]
        return df[[id_col, cazy_col]].rename(columns={id_col: "target", cazy_col: "cazy_id"})

    def _run_foldseek(self, out_m8: Path, tmp_dir: Path) -> None:
        """Run foldseek easy-search. Query = input faa, target = config DB."""
        foldseek_bin = _find_foldseek()
        # easy-search: query target_db output tmp
        # Query: FASTA is supported in recent Foldseek for sequence-to-structure search
        cmd = [
            foldseek_bin, "easy-search",
            str(self.input_faa),
            str(self.db_path),
            str(out_m8),
            str(tmp_dir),
            "-e", str(self.config.e_value_threshold),
            "--threads", str(self.config.threads),
            "--format-mode", "3",  # m8 tabular
        ]
        logger.info("Running Foldseek: %s", " ".join(cmd))
        start = time.time()
        result = subprocess.run(cmd, capture_output=True, text=True)
        elapsed = time.time() - start
        if result.returncode != 0:
            logger.error("Foldseek stderr: %s", result.stderr and result.stderr.strip())
            raise RuntimeError(
                f"Foldseek easy-search failed (exit {result.returncode}). "
                "Ensure Foldseek supports sequence query (FASTA) or provide structure-predicted inputs."
            )
        logger.info("Foldseek completed in %.2fs", elapsed)

    def _hits_to_overview(self, m8_path: Path, mapping: pd.DataFrame) -> pd.DataFrame:
        """Parse m8 and map targets to CAZy IDs; aggregate per query to Gene ID, CAZy ID."""
        if not m8_path.exists() or m8_path.stat().st_size == 0:
            return pd.DataFrame(columns=S.STRUCTURE_SEARCH_COLUMN_NAMES_OVERVIEW)
        df = pd.read_csv(m8_path, sep="\t", header=None, names=FOLDSEEK_M8_COLS)
        df["evalue"] = pd.to_numeric(df["evalue"], errors="coerce").fillna(1.0)
        df = df[df["evalue"] <= float(self.config.e_value_threshold)]
        if df.empty:
            return pd.DataFrame(columns=S.STRUCTURE_SEARCH_COLUMN_NAMES_OVERVIEW)
        # Normalize target to first token (Foldseek may append _1 etc.)
        df["target_base"] = df["target"].astype(str).str.split("_").str[0]
        mapping = mapping.copy()
        mapping["target_base"] = mapping["target"].astype(str).str.split("_").str[0]
        merged = df.merge(mapping[["target_base", "cazy_id"]], on="target_base", how="left")
        merged = merged[merged["cazy_id"].notna() & (merged["cazy_id"].astype(str).str.strip() != "")]
        if merged.empty:
            return pd.DataFrame(columns=S.STRUCTURE_SEARCH_COLUMN_NAMES_OVERVIEW)
        # Best hit per query (by evalue), then aggregate multi-domain as CAZy+CAZy
        merged = merged.sort_values("evalue")
        best = merged.groupby("query", as_index=False).agg({
            "cazy_id": lambda x: "+".join(x.drop_duplicates().astype(str).tolist()),
        })
        best = best.rename(columns={"query": "Gene ID", "cazy_id": "CAZy ID"})
        return best[S.STRUCTURE_SEARCH_COLUMN_NAMES_OVERVIEW]

    def run(self) -> None:
        """Run structure search and write overview-compatible TSV."""
        import tempfile
        out_m8 = self.output_file.with_suffix(self.output_file.suffix + ".m8")
        with tempfile.TemporaryDirectory(prefix="foldseek_") as tmp_dir:
            self._run_foldseek(out_m8, Path(tmp_dir))
        mapping = self._load_mapping()
        overview_df = self._hits_to_overview(out_m8, mapping)
        overview_df.to_csv(self.output_file, sep="\t", index=False)
        logger.info("Structure search results written to %s (%d rows)", self.output_file, len(overview_df))
        # Clean up raw m8 if desired (keep for debugging; optional delete)
        try:
            out_m8.unlink()
        except OSError:
            pass
