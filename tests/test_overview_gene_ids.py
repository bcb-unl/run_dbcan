"""Regression tests for GitHub issue #82: protein IDs must be handled as opaque strings.

Purely numeric IDs (``0``, ``1``), IDs with leading zeros (``007``), scientific-looking IDs
(``1e5``) and IDs such as ``NA``/``nan``/``None`` used to be parsed by pandas as numbers or
missing values, which produced an overview.tsv without annotations (or no overview at all).
"""
from pathlib import Path

import pandas as pd
import pytest

import dbcan.constants.diamond_constants as D
from dbcan.configs.base_config import OverviewGeneratorConfig
from dbcan.configs.pyhmmer_config import DBCANSUBConfig
from dbcan.IO.OverviewGenerator import OverviewGenerator
from dbcan.IO.table import read_tsv
from dbcan.process.process_dbcan_sub import DBCANSUBProcessor
from dbcan.process.process_utils import process_results

HMM_ROW = "GH5_2.hmm\t300\t{id}\t400\t1e-50\t1\t290\t10\t300\t0.96\tdbCAN\n"
SUB_ROW = "GH5_e12.hmm|GH5:10|3.2.1.4:3\t300\t{id}\t400\t1e-50\t1\t290\t10\t300\t0.96\tdbCAN-sub\n"


def _run_pipeline(out: Path, hmm_ids, sub_ids, diamond_ids) -> pd.DataFrame:
    """Drive the real post-processing code paths from raw hit tables to overview.tsv."""
    out.mkdir(parents=True, exist_ok=True)
    raw = out / "hmm.raw.tsv"
    raw.write_text("".join(HMM_ROW.format(id=i) for i in hmm_ids))
    process_results(None, str(out / "dbCAN_hmm_results.tsv"), temp_hits_file=raw)

    raw_sub = out / "sub.raw.tsv"
    raw_sub.write_text("".join(SUB_ROW.format(id=i) for i in sub_ids))
    process_results(None, str(out / "dbCANsub_hmm_raw.tsv"), temp_hits_file=raw_sub)
    DBCANSUBProcessor(DBCANSUBConfig(db_dir=str(out), output_dir=str(out))).process_dbcan_sub()

    rows = [[i, f"ABC{n}.1|GH5_2", 50, 300, 1, 0, 1, 300, 1, 300, 1e-120, 500] for n, i in enumerate(diamond_ids)]
    pd.DataFrame(rows, columns=D.CAZY_COLUMN_NAMES).to_csv(out / "diamond.out", sep="\t", index=False)

    OverviewGenerator(OverviewGeneratorConfig(output_dir=str(out))).run()
    overview = out / "overview.tsv"
    assert overview.exists(), "overview.tsv was not written"
    return read_tsv(overview).set_index("Gene ID")


@pytest.mark.parametrize(
    "ids",
    [
        ["0", "1", "2"],
        ["007", "1e5", "10"],
        ["NA", "nan", "None", "NULL"],
        ["5", "abc", "x|y"],
    ],
    ids=["numeric", "numeric-looking", "na-like", "mixed"],
)
def test_overview_keeps_ids_verbatim_and_annotated(tmp_path, ids):
    ov = _run_pipeline(tmp_path, ids, ids, ids)
    assert sorted(ov.index) == sorted(ids)
    for gid in ids:
        assert ov.loc[gid, "#ofTools"] == 3, gid
        assert ov.loc[gid, "dbCAN_hmm"].startswith("GH5_2("), gid
        assert ov.loc[gid, "dbCAN_sub"].startswith("GH5_e12("), gid
        assert ov.loc[gid, "DIAMOND"] == "GH5_2", gid


def test_overview_mixed_tools_with_numeric_and_text_ids(tmp_path):
    """Numeric IDs in the HMM tables and text IDs in DIAMOND used to crash sorted()."""
    ov = _run_pipeline(tmp_path, ["5"], ["5"], ["abc"])
    assert set(ov.index) == {"5", "abc"}
    assert ov.loc["5", "#ofTools"] == 2
    assert ov.loc["abc", "#ofTools"] == 1


def test_read_tsv_keeps_identifier_strings(tmp_path):
    p = tmp_path / "t.tsv"
    p.write_text("Target Name\tTarget From\ti-Evalue\n007\t1\t1e-5\nNA\t2\t1e-9\n1e5\t3\t0.1\n")
    df = read_tsv(p)
    assert df["Target Name"].tolist() == ["007", "NA", "1e5"]
    assert df["Target From"].dtype.kind == "i"
    assert df["i-Evalue"].dtype.kind == "f"


def test_gff_processing_keeps_numeric_looking_protein_ids(tmp_path):
    """cgc.gff annotation must match overview/total_cgc_info rows for IDs like 00101 or 1.10."""
    from dbcan.configs.base_config import GFFConfig
    from dbcan.IO.gff import get_gff_processor

    (tmp_path / "in.gff").write_text(
        "##gff-version 3\n"
        "c1\tx\tCDS\t100\t1000\t.\t+\t0\tID=00101\n"
        "c1\tx\tCDS\t1100\t2000\t.\t-\t0\tID=00102\n"
        "c1\tx\tCDS\t2100\t3000\t.\t-\t0\tID=1.10\n"
    )
    (tmp_path / "overview.tsv").write_text(
        "Gene ID\tEC#\tdbCAN_hmm\tdbCAN_sub\tDIAMOND\t#ofTools\tRecommend Results\n"
        "00101\t-\tGH5\t-\tGH5\t2\tGH5\n"
        "1.10\t-\tGH10\t-\tGH10\t2\tGH10\n"
    )
    cols = ["Annotate Name", "Annotate Length", "Target Name", "Target Length", "i-Evalue", "Annotate From",
            "Annotate To", "Target From", "Target To", "Coverage", "Annotate File Name", "Type"]
    (tmp_path / "total_cgc_info.tsv").write_text(
        "\t".join(cols) + "\n" + "2.A.1.1.1\t100\t00102\t300\t1e-10\t1\t100\t1\t100\t0.9\tx\tTC\n"
    )
    cfg = GFFConfig(output_dir=str(tmp_path), input_gff=str(tmp_path / "in.gff"), gff_type="prodigal")
    assert get_gff_processor(cfg).process_gff()
    text = (tmp_path / "cgc.gff").read_text()
    assert "protein_id=00101;CGC_annotation=CAZyme|GH5" in text
    assert "protein_id=1.10;CGC_annotation=CAZyme|GH10" in text
    assert "protein_id=00102;CGC_annotation=TC|2.A.1.1.1" in text
