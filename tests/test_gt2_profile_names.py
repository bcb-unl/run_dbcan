"""GT2_* profiles are merged to GT2 only for dbCAN.hmm, never for dbCAN-sub.

dbCAN.hmm up to V14 split GT2 into several GT2_* profiles (GT2_Glycos_transf_2, ...)
that are reported as GT2. dbCAN-sub eHMMs such as ``GT2_e2022.hmm|GT2:518|2.4.1.12:6``
carry the subfamily, composition and EC numbers in their name; renaming them to
``GT2.hmm`` erased that information for every GT2 eHMM.
"""
import random

import pytest

pyhmmer = pytest.importorskip("pyhmmer")
from pyhmmer.easel import Alphabet, TextMSA, TextSequence  # noqa: E402
from pyhmmer.plan7 import Background, Builder  # noqa: E402

from dbcan.annotation.pyhmmer_search import (  # noqa: E402
    PyHMMERDBCANProcessor,
    PyHMMERDBCANSUBProcessor,
)
from dbcan.configs.pyhmmer_config import DBCANSUBConfig, PyHMMERDBCANConfig  # noqa: E402
from dbcan.IO.table import read_tsv  # noqa: E402

ALPHABET = Alphabet.amino()
_rng = random.Random(1)
SEQ = "".join(_rng.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(250))


def _write_hmm(path, name):
    msa = TextMSA(name=name.encode(), sequences=[TextSequence(name=b"s%d" % i, sequence=SEQ) for i in range(3)])
    hmm, _, _ = Builder(ALPHABET).build_msa(msa.digitize(ALPHABET), Background(ALPHABET))
    with open(path, "wb") as fh:
        hmm.write(fh)


@pytest.fixture
def workdir(tmp_path):
    db = tmp_path / "db"
    db.mkdir()
    _write_hmm(db / "dbCAN.hmm", "GT2_Glycos_transf_2.hmm")
    _write_hmm(db / "dbCAN-sub.hmm", "GT2_e2022.hmm|GT2:518|2.4.1.12:6")
    (db / "fam-substrate-mapping.tsv").write_text(
        "Substrate_high_level\tSubstrate_curated\tFamily\tName\tEC_Number\n"
    )
    (tmp_path / "uniInput.faa").write_text(f">p1\n{SEQ}\n")
    return tmp_path, db


def test_dbcan_hmm_merges_gt2_profiles(workdir):
    out, db = workdir
    PyHMMERDBCANProcessor(PyHMMERDBCANConfig(db_dir=str(db), output_dir=str(out), threads=1)).run()
    df = read_tsv(out / "dbCAN_hmm_results.tsv")
    assert df["HMM Name"].tolist() == ["GT2.hmm"]


def test_dbcan_sub_keeps_gt2_subfamily_and_ec(workdir):
    out, db = workdir
    PyHMMERDBCANSUBProcessor(DBCANSUBConfig(db_dir=str(db), output_dir=str(out), threads=1)).run()
    df = read_tsv(out / "dbCANsub_hmm_results.tsv")
    assert df["Subfam Name"].tolist() == ["GT2_e2022"]
    assert df["Subfam Composition"].tolist() == ["GT2:518"]
    assert df["Subfam EC"].tolist() == ["2.4.1.12:6"]


def test_output_hmm_name_only_merges_for_dbcan_hmm():
    assert PyHMMERDBCANProcessor.MERGE_GT2_PROFILES is True
    assert PyHMMERDBCANSUBProcessor.MERGE_GT2_PROFILES is False
