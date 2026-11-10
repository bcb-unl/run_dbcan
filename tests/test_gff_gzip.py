"""gff_process reads gzip GFF (``.gz``) for every legacy --gff_type, like the gene/protein match pre-check."""
import gzip
from pathlib import Path

import pytest

import dbcan.constants.gff_constants as G
from dbcan.configs.base_config import GFFConfig
from dbcan.IO.gff import _open_gff_text, get_gff_processor
from dbcan.IO.gff_resolver import check_gff_resolution, read_fasta_ids

FX = Path(__file__).parent / "_data" / "gff_dialects"

NCBI_PROK_GFF = """##gff-version 3
NC_1.1\tRefSeq\tgene\t500\t900\t.\t-\t.\tID=gene-b2;locus_tag=b2
NC_1.1\tRefSeq\tCDS\t500\t900\t.\t-\t0\tID=cds-WP_2.1;Parent=gene-b2;protein_id=WP_2.1
NC_1.1\tRefSeq\tgene\t10\t300\t.\t+\t.\tID=gene-b1;locus_tag=b1
NC_1.1\tRefSeq\tCDS\t10\t300\t.\t+\t0\tID=cds-WP_1.1;Parent=gene-b1;protein_id=WP_1.1
NC_1.1\tRefSeq\tgene\t1000\t1076\t.\t+\t.\tID=gene-b3;locus_tag=b3
NC_1.1\tRefSeq\ttRNA\t1000\t1076\t.\t+\t.\tID=rna-b3;Parent=gene-b3
"""

PRODIGAL_GFF = """##gff-version 3
##sequence-region c1 1 5000
c1\tProdigal_v2.6.3\tCDS\t2\t100\t.\t+\t0\tID=p1;partial=00;
c1\tProdigal_v2.6.3\tCDS\t200\t350\t.\t-\t0\tID=p2;partial=00;
"""

JGI_GFF = """##gff-version 3
scaffold_1\tJGI\tgene\t10\t300\t.\t+\t.\tID=g1;proteinId=242238
scaffold_1\tJGI\tgene\t500\t900\t.\t-\t.\tID=g2;proteinId=242239
"""

CGC = {"WP_1.1": "CAZyme|GH5", "p1": "CAZyme|GH5", "jgi-Xylhe1-242238-fg_1": "TC|1.A.1"}


def _write(tmp_path: Path, name: str, text: str, gz: bool) -> Path:
    p = tmp_path / name
    if gz:
        with gzip.open(p, "wt") as fh:
            fh.write(text)
    else:
        p.write_text(text)
    return p


def _run(tmp_path: Path, gff: Path, gff_type: str, faa_ids, **cfg) -> str:
    out = tmp_path / ("out_" + gff.name.replace(".", "_"))
    out.mkdir()
    (out / G.GFF_INPUT_PROTEIN_FILE).write_text("".join(f">{i}\nM\n" for i in faa_ids))
    proc = get_gff_processor(GFFConfig(output_dir=str(out), input_gff=str(gff), gff_type=gff_type, **cfg))
    proc.load_cgc_type = lambda: {k: {G.GFF_CGC_ANNOTATION_COL: v} for k, v in CGC.items()}
    assert proc.process_gff()
    return (out / G.GFF_OUTPUT_FILE).read_text()


@pytest.mark.parametrize("name,text,gff_type,ids,cfg", [
    ("prok.gff", NCBI_PROK_GFF, G.GFF_FORMAT_NCBI_PROK, ["WP_1.1", "WP_2.1"], {}),
    ("prodigal.gff", PRODIGAL_GFF, G.GFF_FORMAT_PRODIGAL, ["p1", "p2"], {"prodigal_gff_streaming": "off"}),
    ("prodigal_s.gff", PRODIGAL_GFF, G.GFF_FORMAT_PRODIGAL, ["p1", "p2"], {"prodigal_gff_streaming": "on"}),
    ("jgi.gff", JGI_GFF, G.GFF_FORMAT_JGI, ["jgi-Xylhe1-242238-fg_1", "jgi-Xylhe1-242239-fg_2"], {}),
])
def test_gz_gff_gives_same_cgc_gff_as_plain(tmp_path, name, text, gff_type, ids, cfg):
    plain = _run(tmp_path, _write(tmp_path, name, text, gz=False), gff_type, ids, **cfg)
    packed = _run(tmp_path, _write(tmp_path, name + ".gz", text, gz=True), gff_type, ids, **cfg)
    assert plain and packed == plain
    assert "CGC_annotation=CAZyme|GH5" in packed or "CGC_annotation=TC|1.A.1" in packed


def test_ncbi_prok_gz_annotates_cds_and_sorts(tmp_path):
    out = _run(tmp_path, _write(tmp_path, "prok.gff.gz", NCBI_PROK_GFF, gz=True),
               G.GFF_FORMAT_NCBI_PROK, ["WP_1.1", "WP_2.1"])
    assert out.splitlines() == [
        "NC_1.1\t.\tgene\t10\t300\t.\t+\t.\tprotein_id=WP_1.1;CGC_annotation=CAZyme|GH5",
        "NC_1.1\t.\tgene\t500\t900\t.\t-\t.\tprotein_id=WP_2.1;CGC_annotation=null",
        "NC_1.1\t.\tgene\t1000\t1076\t.\t+\t.\tprotein_id=NA;CGC_annotation=Other|tRNA",
    ]


def test_ncbi_euk_gz_fixture(tmp_path):
    ids = read_fasta_ids(FX / "ncbi_euk.faa")
    text = (FX / "ncbi_euk.gff3").read_text()
    plain = _run(tmp_path, _write(tmp_path, "euk.gff3", text, gz=False), G.GFF_FORMAT_NCBI_EUK, ids)
    packed = _run(tmp_path, _write(tmp_path, "euk.gff3.gz", text, gz=True), G.GFF_FORMAT_NCBI_EUK, ids)
    assert packed == plain and packed.count("\n") == 6


def test_gz_passes_precheck_and_gff_process(tmp_path):
    """The pre-check already read .gz; gff_process now reads the same file."""
    gff = _write(tmp_path, "prok.gff.gz", NCBI_PROK_GFF, gz=True)
    res = check_gff_resolution(gff, G.GFF_FORMAT_NCBI_PROK, ["WP_1.1", "WP_2.1"])
    assert res.n_matched == 2
    assert "protein_id=WP_2.1" in _run(tmp_path, gff, G.GFF_FORMAT_NCBI_PROK, ["WP_1.1", "WP_2.1"])


def test_open_gff_text_detects_gzip_by_magic_bytes(tmp_path):
    gff = _write(tmp_path, "no_suffix.gff", PRODIGAL_GFF, gz=True)
    with _open_gff_text(gff) as fh:
        assert fh.read() == PRODIGAL_GFF
    plain = _write(tmp_path, "plain.gff", PRODIGAL_GFF, gz=False)
    with _open_gff_text(plain) as fh:
        assert fh.read() == PRODIGAL_GFF
    with pytest.raises(FileNotFoundError):
        _open_gff_text(tmp_path / "missing.gff")
