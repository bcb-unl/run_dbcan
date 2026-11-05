"""Helpers for reading the TSV files exchanged between run_dbcan stages.

Protein / gene / contig identifiers must always be treated as opaque strings.
With pandas defaults, IDs such as ``0``, ``007`` or ``1e5`` are parsed as numbers
(``007`` -> ``7``) and IDs such as ``NA``, ``nan`` or ``None`` become missing
values, which breaks the joins between result tables (GitHub issue #82).
"""
from __future__ import annotations

from typing import Iterable, Optional

import pandas as pd

# Columns that hold identifiers in run_dbcan intermediate/result tables.
ID_COLUMNS = (
    "Gene ID",
    "Target Name",
    "HMM Name",
    "Subfam Name",
    "CAZy ID",
    "TCDB ID",
    "TF ID",
    "Sul ID",
    "Peptidase ID",
    "Target ID",
    "Protein ID",
    "protein_id",
    "Contig ID",
    "Contig",
)


def read_tsv(path, id_columns: Optional[Iterable[str]] = None, **kwargs):
    """``pandas.read_csv`` for run_dbcan TSV files with identifier-safe defaults.

    * ``id_columns`` (default :data:`ID_COLUMNS`) are read as ``str`` when present.
    * Only empty fields are treated as missing values, so strings such as
      ``NA``/``nan``/``None`` are kept verbatim in every column.

    Other keyword arguments are passed to :func:`pandas.read_csv` (for example
    ``header``, ``names``, ``usecols`` or ``chunksize``).
    """
    cols = tuple(ID_COLUMNS if id_columns is None else id_columns)
    dtype = kwargs.pop("dtype", None)
    if dtype is None or isinstance(dtype, dict):
        merged = {c: str for c in cols}
        if dtype:
            merged.update(dtype)
        dtype = merged
    kwargs.setdefault("sep", "\t")
    kwargs.setdefault("keep_default_na", False)
    kwargs.setdefault("na_values", [""])
    return pd.read_csv(path, dtype=dtype, **kwargs)
