"""--extend_mode gene/bp: extensions never overlap a neighbouring CGC, so each gene is in at most one CGC."""
import pandas as pd
import pytest

from dbcan.annotation.CGCFinder import CGCFinder
from dbcan.configs.cgcfinder_config import CGCFinderConfig

OVERVIEW_HEADER = "Gene ID\tEC#\tdbCAN_hmm\tdbCAN_sub\tDIAMOND\t#ofTools\tRecommend Results\tSubstrate\n"
SIG_HEADER = ("Annotate Name\tAnnotate Length\tTarget Name\tTarget Length\ti-Evalue\tAnnotate From\tAnnotate To\t"
              "Target From\tTarget To\tCoverage\tAnnotate File Name\tType\tSource\n")


def _setup(tmp_path, layout):
    """layout: string, one char per gene 1 kb apart on contig c1: C = CAZyme, T = TC, . = null."""
    ann = {"C": "CAZyme|GH5_e1", "T": "TC|2.A.1.1.1", ".": "null"}
    gff = "".join(
        f"c1\t.\tCDS\t{1000 * i + 1}\t{1000 * i + 900}\t.\t+\t.\tprotein_id=g{i};CGC_annotation={ann[c]}\n"
        for i, c in enumerate(layout)
    )
    (tmp_path / "cgc.gff").write_text(gff)
    (tmp_path / "overview.tsv").write_text(OVERVIEW_HEADER)
    (tmp_path / "total_cgc_info.tsv").write_text(SIG_HEADER)


def _cgcs(tmp_path, **cfg):
    assert CGCFinder(CGCFinderConfig(output_dir=str(tmp_path), additional_genes=["TC"], **cfg)).run()
    out = pd.read_csv(tmp_path / "cgc_standard_out.tsv", sep="\t", dtype=str)
    return [[int(p[1:]) for p in genes] for genes in out.groupby("CGC#", sort=False)["Protein ID"].apply(list)]


#           0123456789
LAYOUT = "CT...CT..."  # two CGCs (genes 0-1 and 5-6) separated by 3 null genes


def test_no_extension(tmp_path):
    _setup(tmp_path, LAYOUT)
    assert _cgcs(tmp_path) == [[0, 1], [5, 6]]


@pytest.mark.parametrize("cfg, expected", [
    # 2 genes each side: CGC1 would take 2-3, CGC2 would take 3-4 -> gene 3 is shared; equal distance -> CGC1
    (dict(extend_mode="gene", extend_gene_count=2), [[0, 1, 2, 3], [4, 5, 6, 7, 8]]),
    # 10 genes each side: CGC1 must stop before CGC2's core (gene 5) and CGC2 before CGC1's core
    (dict(extend_mode="gene", extend_gene_count=10), [[0, 1, 2, 3], [4, 5, 6, 7, 8, 9]]),
    # bp: genes that start/end within 2.5 kb of the core (3 genes each side here); shared genes 2-4 split
    (dict(extend_mode="bp", extend_bp=2500), [[0, 1, 2, 3], [4, 5, 6, 7, 8, 9]]),
])
def test_extensions_do_not_overlap(tmp_path, cfg, expected):
    _setup(tmp_path, LAYOUT)
    cgcs = _cgcs(tmp_path, **cfg)
    assert cgcs == expected
    genes = [g for cgc in cgcs for g in cgc]
    assert len(genes) == len(set(genes))


def test_extension_not_shared_keeps_full_length(tmp_path):
    #           0         1
    #           0123456789012345
    _setup(tmp_path, "CT.........CT...")
    assert _cgcs(tmp_path, extend_mode="gene", extend_gene_count=2) == [[0, 1, 2, 3], [9, 10, 11, 12, 13, 14]]
