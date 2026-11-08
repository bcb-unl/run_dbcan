"""Eukaryotic GFF3/GTF gene -> protein resolution (--gff_type funannotate/BRAKER/GFF3/GTF).

Fixtures in tests/_data/gff_dialects: 60 kb of Xylona heveae KV407454.1 (NCBI annotation,
5 protein-coding genes + 1 tRNA) rewritten in the ID schemes of funannotate2 (gfftk), BRAKER
(gtf2gff.pl --gff3, with a second, shorter isoform of the first gene), AUGUSTUS GTF, Ensembl,
Helixer and a FASTA named by gene ID; plus 4 genes of a real AUGUSTUS 3.5.0 GFF3 (``transcript``
rows, CDS ``g1.t1.cds``) and the original NCBI GFF3 of the same region.
"""
import csv
from pathlib import Path

import pytest

import dbcan.constants.gff_constants as G
from dbcan.configs.base_config import GFFConfig
from dbcan.core import DbcanStepError, check_gff_protein_match
from dbcan.IO.gff import ResolverGFFProcessor, get_gff_processor
from dbcan.IO.gff_resolver import (
    GffResolutionError,
    ProteinIndex,
    check_gff_resolution,
    normalize_gff_type,
    rank_gff_types,
    read_fasta_ids,
    resolve_gff,
    write_cgc_gff,
)

DATA = Path(__file__).parent / "_data"
FX = DATA / "gff_dialects"

DIALECTS = [
    ("funannotate2", "gff3", G.GFF_FORMAT_FUNANNOTATE),
    ("braker", "gff3", G.GFF_FORMAT_BRAKER),
    ("augustus_gtf", "gtf", G.GFF_FORMAT_GTF),
    ("ensembl", "gff3", G.GFF_FORMAT_GFF3),
    ("helixer", "gff3", G.GFF_FORMAT_GFF3),
    ("gene_named", "gff3", G.GFF_FORMAT_FUNANNOTATE),
]


def _truth(name):
    with open(FX / f"{name}.truth.tsv") as fh:
        return {r["gene_id"]: r["protein_id"] for r in csv.DictReader(fh, delimiter="\t")}


@pytest.mark.parametrize("name,ext,gff_type", DIALECTS)
def test_dialect_resolves_every_coding_gene(name, ext, gff_type):
    res = resolve_gff(FX / f"{name}.{ext}", gff_type, read_fasta_ids(FX / f"{name}.faa"))
    got = {r.gene_id: r.protein_id for r in res.rows}
    truth = _truth(name)
    assert len(truth) == 5
    assert all(got[g] == p for g, p in truth.items())
    assert res.n_coding == res.n_matched == 5
    # the tRNA gene is kept as a null row, as NCBI_euk does
    assert res.n_noncoding == 1
    assert [r.feature_type for r in res.rows if r.status == G.GFF_STATUS_NONCODING] == (
        ["tRNA"] if ext == "gff3" else ["transcript"]
    )


@pytest.mark.parametrize("name,ext,_t", DIALECTS)
def test_generic_gff3_type_also_resolves_every_dialect(name, ext, _t):
    res = resolve_gff(FX / f"{name}.{ext}", "GFF3", read_fasta_ids(FX / f"{name}.faa"))
    assert res.match_rate == 1.0


def test_real_augustus_gff3_with_transcript_rows():
    res = resolve_gff(FX / "augustus_native.gff3", "BRAKER", read_fasta_ids(FX / "augustus_native.faa"))
    assert [(r.gene_id, r.protein_id) for r in res.rows] == [(f"g{i}", f"g{i}.t1") for i in range(1, 5)]
    assert {r.matched_by for r in res.rows} == {"transcript.ID"}


def test_isoforms_one_row_per_gene_longest_cds():
    res = resolve_gff(FX / "braker.gff3", "BRAKER", read_fasta_ids(FX / "braker.faa"))
    multi = [r for r in res.rows if r.n_transcripts > 1]
    assert len(multi) == 1 and len(res.rows) == 6
    assert multi[0].protein_id == "g75.t1"
    assert multi[0].isoform_proteins == ("g75.t1", "g75.t2")


def test_isoform_annotated_as_cazyme_is_preferred():
    ids = read_fasta_ids(FX / "braker.faa")
    res = resolve_gff(FX / "braker.gff3", "BRAKER", ids, annotations={"g75.t2": "CAZyme|GH5"})
    assert [r.protein_id for r in res.rows if r.gene_id == "g75"] == ["g75.t2"]
    # a CAZyme beats another signature on the other isoform
    res = resolve_gff(FX / "braker.gff3", "BRAKER", ids,
                      annotations={"g75.t1": "TC|2.A.1.1.1", "g75.t2": "CAZyme|GH5"})
    assert [r.protein_id for r in res.rows if r.gene_id == "g75"] == ["g75.t2"]
    # both isoforms are CAZymes: the longest CDS wins, even if the shorter one has an extra hit
    res = resolve_gff(FX / "braker.gff3", "BRAKER", ids,
                      annotations={"g75.t1": "CAZyme|GH5", "g75.t2": "CAZyme|GH5+TC|9.B.464.1.2"})
    assert [r.protein_id for r in res.rows if r.gene_id == "g75"] == ["g75.t1"]


def test_ncbi_euk_parity_cgc_gff_byte_identical(tmp_path):
    """The generic GFF3 resolver writes the same cgc.gff as the NCBI_euk processor."""
    ids = read_fasta_ids(FX / "ncbi_euk.faa")
    (tmp_path / "uniInput.faa").write_text("".join(f">{i}\nM\n" for i in ids))
    proc = get_gff_processor(GFFConfig(output_dir=str(tmp_path), input_gff=str(FX / "ncbi_euk.gff3"),
                                       gff_type="NCBI_euk"))
    assert proc.process_gff()
    write_cgc_gff(resolve_gff(FX / "ncbi_euk.gff3", "GFF3", ids), tmp_path / "resolver.gff")
    expected = (tmp_path / "cgc.gff").read_text()
    assert expected.count("\n") == 6 and "Other|tRNA" in expected
    assert (tmp_path / "resolver.gff").read_text() == expected


@pytest.mark.parametrize("gff,faa,gff_type", [
    ("ncbi_euk.gff3", "ncbi_euk.faa", "NCBI_euk"),
    ("funannotate2.gff3", "funannotate2.faa", "prodigal"),
])
def test_report_for_original_types_matches_processor_output(tmp_path, gff, faa, gff_type):
    """For prodigal/NCBI/JGI the report mirrors the unchanged processor (same protein IDs)."""
    ids = read_fasta_ids(FX / faa)
    (tmp_path / "uniInput.faa").write_text("".join(f">{i}\nM\n" for i in ids))
    proc = get_gff_processor(GFFConfig(output_dir=str(tmp_path), input_gff=str(FX / gff), gff_type=gff_type))
    assert proc.process_gff()
    written = sorted(l.split("\t")[8].split(";")[0].split("=", 1)[1]
                     for l in (tmp_path / "cgc.gff").read_text().splitlines())
    tsv = tmp_path / "report.tsv"
    resolve_gff(FX / gff, gff_type, ids, out_tsv=tsv)
    reported = sorted(r["protein_id"] for r in csv.DictReader(open(tsv), delimiter="\t"))
    assert reported == written


def test_jgi_report_uses_protein_id_third_field():
    idx = ProteinIndex(["jgi-Xylhe1-242238-fgenesh1_kg.1_#_17"])
    assert idx.lookup("242238") == "jgi-Xylhe1-242238-fgenesh1_kg.1_#_17"
    # keys go through the --mode protein ID cleaning ('|' -> '-' for JGI IDs)
    assert idx.lookup("jgi|Xylhe1|242238|fgenesh1_kg.1_#_17") == "jgi-Xylhe1-242238-fgenesh1_kg.1_#_17"


def test_version_alias_is_dropped_when_ambiguous():
    idx = ProteinIndex(["XP_1.1", "XP_1.2"])
    assert idx.lookup("XP_1") is None and idx.lookup("XP_1.2") == "XP_1.2"
    assert ProteinIndex(["XP_1.1"]).lookup("XP_1") == "XP_1.1"


def test_url_decoding_multi_parent_and_cds_prefix(tmp_path):
    gff = tmp_path / "x.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "c1\tt\tgene\t1\t900\t.\t+\t.\tID=g%3B1\n"
        "c1\tt\tmRNA\t1\t900\t.\t+\t.\tID=m1;Parent=g%3B1\n"
        "c1\tt\tmRNA\t1\t600\t.\t+\t.\tID=m2;Parent=g%3B1\n"
        "c1\tt\tCDS\t10\t300\t.\t+\t0\tID=cds-P1;Parent=m1,m2\n"
        "c1\tt\tCDS\t400\t600\t.\t+\t0\tID=cds-P1;Parent=m1,m2\n"
        "c1\tt\tCDS\t700\t890\t.\t+\t0\tID=cds-P1;Parent=m1\n")
    faa = tmp_path / "x.faa"
    faa.write_text(">P1\nMK\n>m2\nMK\n")
    res = resolve_gff(gff, "GFF3", read_fasta_ids(faa))
    assert [(r.gene_id, r.protein_id, r.n_transcripts, r.matched_by) for r in res.rows] == [
        ("g;1", "P1", 2, "CDS.ID~")]


def test_gtf_without_gene_rows_and_noncoding_transcript(tmp_path):
    """gffread/NCBI-style GTF: no gene rows, protein_id on CDS, a tRNA transcript without CDS."""
    gtf = tmp_path / "x.gtf"
    gtf.write_text(
        'c1\tsrc\ttranscript\t1\t500\t.\t+\t.\ttranscript_id "t1"; gene_id "G1";\n'
        'c1\tsrc\tCDS\t1\t300\t.\t+\t0\ttranscript_id "t1"; gene_id "G1"; protein_id "XP_9.1";\n'
        'c1\tsrc\ttranscript\t600\t680\t.\t-\t.\ttranscript_id "t2"; gene_id "G2"; transcript_biotype "tRNA";\n'
        'c1\tsrc\texon\t600\t680\t.\t-\t.\ttranscript_id "t2"; gene_id "G2";\n')
    faa = tmp_path / "x.faa"
    faa.write_text(">XP_9.1 some protein\nMK\n")
    res = resolve_gff(gtf, "gtf", read_fasta_ids(faa))
    assert [(r.gene_id, r.protein_id, r.status, r.feature_type) for r in res.rows] == [
        ("G1", "XP_9.1", "matched", "gene"), ("G2", "NA", "noncoding", "tRNA")]


def test_gff_type_is_case_insensitive_and_validated():
    assert normalize_gff_type("braker") == "BRAKER"
    assert normalize_gff_type("ncbi_EUK") == "NCBI_euk"
    with pytest.raises(ValueError, match="Choose one of"):
        normalize_gff_type("maker")


def test_wrong_type_fails_with_suggestion(tmp_path):
    """funannotate2 GFF read as NCBI_euk: 0 genes matched -> error naming the right type."""
    ids = read_fasta_ids(FX / "funannotate2.faa")
    with pytest.raises(GffResolutionError) as exc:
        check_gff_resolution(FX / "funannotate2.gff3", "NCBI_euk", ids, 0.9, out_tsv=tmp_path / "r.tsv")
    msg = str(exc.value)
    assert "0/5 protein-coding genes" in msg
    assert "FUN_000075" in msg                        # an unmatched example
    assert "Suggested: --gff_type" in msg
    assert isinstance(exc.value, DbcanStepError)
    assert (tmp_path / "r.tsv").read_text().count("unmatched") == 5
    best = rank_gff_types(FX / "funannotate2.gff3", ids)[0]
    assert best[1] == 1.0


def test_fasta_from_another_annotation_fails(tmp_path):
    with pytest.raises(GffResolutionError, match="No --gff_type matches"):
        check_gff_resolution(FX / "braker.gff3", "BRAKER", read_fasta_ids(FX / "helixer.faa"), 0.9)


def test_min_match_zero_only_warns(tmp_path):
    res = check_gff_resolution(FX / "braker.gff3", "BRAKER", read_fasta_ids(FX / "helixer.faa"), 0)
    assert res.n_matched == 0


def test_resolution_tsv_columns(tmp_path):
    res = check_gff_resolution(FX / "braker.gff3", "BRAKER", read_fasta_ids(FX / "braker.faa"),
                               out_tsv=tmp_path / G.GFF_RESOLUTION_FILE)
    rows = list(csv.DictReader(open(tmp_path / G.GFF_RESOLUTION_FILE), delimiter="\t"))
    assert list(rows[0].keys()) == G.GFF_RESOLUTION_COLUMNS
    assert rows[0]["gene_id"] == "g75" and rows[0]["isoform_proteins"] == "g75.t1,g75.t2"
    assert rows[0]["matched_by"] == "transcript.ID"
    assert [r["status"] for r in rows].count("noncoding") == 1
    assert res.match_rate == 1.0


def test_resolver_processor_writes_cgc_gff_and_report(tmp_path):
    ids = read_fasta_ids(FX / "funannotate2.faa")
    (tmp_path / "uniInput.faa").write_text("".join(f">{i}\nM\n" for i in ids))
    (tmp_path / "overview.tsv").write_text(
        "Gene ID\tEC#\tdbCAN_hmm\tdbCAN_sub\tDIAMOND\t#ofTools\tRecommend Results\n"
        "FUN_000077-T1\t-\tGH5(1-200)\t-\tGH5\t2\tGH5\n")
    (tmp_path / "total_cgc_info.tsv").write_text(
        "Annotate Name\tAnnotate Length\tTarget Name\tTarget Length\ti-Evalue\tAnnotate From\tAnnotate To\t"
        "Target From\tTarget To\tCoverage\tAnnotate File Name\tType\n"
        "2.A.1.1.1\t500\tFUN_000078-T1\t400\t1e-50\t1\t400\t1\t400\t0.9\tdiamond.out.tc\tTC\n")
    proc = get_gff_processor(GFFConfig(output_dir=str(tmp_path), input_gff=str(FX / "funannotate2.gff3"),
                                       gff_type="FUNANNOTATE"))
    assert isinstance(proc, ResolverGFFProcessor) and proc.config.gff_type == "funannotate"
    assert proc.process_gff()
    lines = (tmp_path / "cgc.gff").read_text().splitlines()
    assert len(lines) == 6
    assert any("protein_id=FUN_000077-T1;CGC_annotation=CAZyme|GH5" in l for l in lines)
    assert any("protein_id=FUN_000078-T1;CGC_annotation=TC|2.A.1.1.1" in l for l in lines)
    assert any("protein_id=NA;CGC_annotation=Other|tRNA" in l for l in lines)
    assert (tmp_path / G.GFF_RESOLUTION_FILE).exists()


def test_check_skipped_for_pyrodigal_gff(tmp_path):
    (tmp_path / "uniInput.gff").write_text("c1\tpyrodigal\tCDS\t1\t300\t.\t+\t0\tID=c1_1;\n")
    (tmp_path / "uniInput.faa").write_text(">other\nMK\n")
    cfg = GFFConfig(output_dir=str(tmp_path), input_gff=str(tmp_path / "uniInput.gff"), gff_type="prodigal")
    assert check_gff_protein_match(cfg) is None


def test_check_runs_for_user_gff(tmp_path):
    (tmp_path / "uniInput.faa").write_text(">other\nMK\n")
    cfg = GFFConfig(output_dir=str(tmp_path), input_gff=str(FX / "braker.gff3"), gff_type="BRAKER")
    with pytest.raises(GffResolutionError):
        check_gff_protein_match(cfg)
    cfg.gff_min_match = 0.0
    assert check_gff_protein_match(cfg).n_matched == 0


@pytest.mark.parametrize("gff,faa,gff_type,matched,noncoding", [
    ("EscheriaColiK12MG1655.gff", "EscheriaColiK12MG1655.faa", "NCBI_prok", 4140, 173),
    ("Xylhe1_GeneCatalog_proteins_20130827.gff", "Xylhe1_GeneCatalog_proteins_20130827.aa.fasta", "JGI", 8205, 0),
])
def test_report_for_original_types_on_full_genomes(tmp_path, gff, faa, gff_type, matched, noncoding):
    ids = read_fasta_ids(DATA / faa)
    (tmp_path / "uniInput.faa").write_text("".join(f">{i}\nM\n" for i in ids))
    proc = get_gff_processor(GFFConfig(output_dir=str(tmp_path), input_gff=str(DATA / gff), gff_type=gff_type))
    assert proc.process_gff()
    written = sorted(l.split("\t")[8].split(";")[0].split("=", 1)[1]
                     for l in (tmp_path / "cgc.gff").read_text().splitlines())
    res = resolve_gff(DATA / gff, gff_type, ids)
    assert sorted(r.protein_id for r in res.rows) == written
    assert (res.n_matched, res.n_noncoding) == (matched, noncoding)
    # the generic GFF3 resolver finds the same proteins on the JGI file (proteinId on mRNA)
    if gff_type == "JGI":
        assert resolve_gff(DATA / gff, "GFF3", ids).n_matched == matched


def test_prodigal_report_is_streamed(tmp_path):
    gff = tmp_path / "p.gff"
    gff.write_text("# Sequence Data\n"
                   "c1\tpyrodigal\tCDS\t1\t300\t.\t+\t0\tID=c1_1;partial=00\n"
                   "c1\tpyrodigal\tCDS\t400\t900\t.\t-\t0\tID=c1_2;partial=00\n"
                   "c2\tpyrodigal\tCDS\t1\t90\t.\t+\t0\tID=c2_1;partial=10\n")
    tsv = tmp_path / G.GFF_RESOLUTION_FILE
    res = resolve_gff(gff, "PRODIGAL", ["c1_1", "c1_2", "x"], out_tsv=tsv)
    assert (res.n_coding, res.n_matched, res.n_noncoding, res.fasta_placed) == (3, 2, 0, 2)
    assert res.rows == [] and [r.gene_id for r in res.unmatched_examples()] == ["c2_1"]
    lines = tsv.read_text().splitlines()
    assert len(lines) == 4 and lines[3].split("\t")[7] == "unmatched"
    with pytest.raises(GffResolutionError, match="1 FASTA proteins not placed|Only 2/3"):
        check_gff_resolution(gff, "prodigal", ["c1_1", "c1_2", "x"], 0.9)


def test_gffread_transcripts_without_gene_rows_are_grouped_by_gene_id(tmp_path):
    """gffread -E/--keep-genes off: mRNA rows have no Parent; geneID=/gene_id=/gene= name the gene."""
    gff = tmp_path / "x.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "c1\tgffread\tmRNA\t100\t900\t.\t+\t.\tID=t1;geneID=G1\n"
        "c1\tgffread\tCDS\t100\t900\t.\t+\t0\tParent=t1\n"
        "c1\tgffread\tmRNA\t50\t600\t.\t+\t.\tID=t2;geneID=G1\n"
        "c1\tgffread\tCDS\t150\t600\t.\t+\t0\tParent=t2\n"
        "c1\tgffread\ttranscript\t40\t950\t.\t+\t.\tID=t3;geneID=G1\n"          # non-coding isoform
        "c1\tgffread\tmRNA\t2000\t2600\t.\t-\t.\tID=t4;gene_id=G2\n"
        "c1\tgffread\tCDS\t2000\t2600\t.\t-\t0\tParent=t4\n"
        "c1\tgffread\tmRNA\t3000\t3600\t.\t-\t.\tID=t5;gene=abcD\n"
        "c1\tgffread\tCDS\t3000\t3600\t.\t-\t0\tParent=t5\n"
        "c1\tgffread\ttranscript\t4000\t4100\t.\t+\t.\tID=t6;geneID=G4\n"
        "c1\tgffread\ttranscript\t4000\t4120\t.\t+\t.\tID=t7;geneID=G4\n")
    faa = tmp_path / "x.faa"
    faa.write_text(">t1\nMK\n>t2\nMK\n>t4\nMK\n>t5\nMK\n")
    res = resolve_gff(gff, "GFF3", read_fasta_ids(faa))
    assert [(r.gene_id, r.protein_id, r.status, r.n_transcripts, r.start, r.end) for r in res.rows] == [
        ("G1", "t1", "matched", 2, 50, 900),
        ("G2", "t4", "matched", 1, 2000, 2600),
        ("abcD", "t5", "matched", 1, 3000, 3600),
        ("G4", "NA", "noncoding", 0, 4000, 4100),
    ]
    assert res.rows[0].isoform_proteins == ("t1", "t2")


def test_pseudogenes_are_noncoding_rows(tmp_path):
    gff = tmp_path / "x.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "c1\tRefSeq\tgene\t1\t300\t.\t+\t.\tID=gene1;Name=aaa\n"
        "c1\tRefSeq\tCDS\t1\t300\t.\t+\t0\tID=cds1;Parent=gene1;protein_id=WP_1.1\n"
        "c1\tRefSeq\tpseudogene\t400\t900\t.\t+\t.\tID=gene2;Name=bbb;pseudo=true\n"
        "c1\tRefSeq\tCDS\t400\t900\t.\t+\t0\tID=cds2;Parent=gene2;pseudo=true\n"
        "c1\tRefSeq\tgene\t1000\t1300\t.\t-\t.\tID=gene3;Name=ccc\n"
        "c1\tRefSeq\tCDS\t1000\t1300\t.\t-\t0\tID=cds3;Parent=gene3;pseudo=true\n"
        "c1\tRefSeq\tpseudogene\t1400\t1500\t.\t-\t.\tID=gene4;Name=ddd\n")
    faa = tmp_path / "x.faa"
    faa.write_text(">WP_1.1\nMK\n>bbb\nMK\n")     # 'bbb' must not be taken through the Name of a pseudogene
    res = resolve_gff(gff, "GFF3", read_fasta_ids(faa))
    assert [(r.gene_id, r.protein_id, r.status, r.feature_type) for r in res.rows] == [
        ("gene1", "WP_1.1", "matched", "gene"),
        ("gene2", "NA", "noncoding", "pseudogene"),
        ("gene3", "NA", "noncoding", "pseudogene"),
        ("gene4", "NA", "noncoding", "pseudogene"),
    ]
    assert res.match_rate == 1.0
    out = tmp_path / "cgc.gff"
    write_cgc_gff(res, out)
    assert out.read_text().count("CGC_annotation=Other|pseudogene") == 3


def test_ecoli_ncbi_gff_with_generic_gff3_type_matches_every_coding_gene():
    ids = read_fasta_ids(DATA / "EscheriaColiK12MG1655.faa")
    res = resolve_gff(DATA / "EscheriaColiK12MG1655.gff", "GFF3", ids)
    assert res.n_matched == res.n_coding == 4140
    assert sum(r.feature_type == "pseudogene" for r in res.rows) == 184


def test_exact_protein_id_wins_over_an_earlier_weak_key(tmp_path):
    """gene A's Name equals the protein_id of gene B (later in the file): B keeps its protein."""
    gff = tmp_path / "x.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "c1\tt\tgene\t1\t300\t.\t+\t.\tID=gA;Name=PB\n"
        "c1\tt\tmRNA\t1\t300\t.\t+\t.\tID=mA;Parent=gA;Name=PB\n"
        "c1\tt\tCDS\t1\t300\t.\t+\t0\tID=cA;Parent=mA;protein_id=PA\n"
        "c1\tt\tgene\t400\t900\t.\t+\t.\tID=gB\n"
        "c1\tt\tmRNA\t400\t900\t.\t+\t.\tID=mB;Parent=gB\n"
        "c1\tt\tCDS\t400\t900\t.\t+\t0\tID=cB;Parent=mB;protein_id=PB\n"
        "c1\tt\tgene\t1000\t1300\t.\t+\t.\tID=gC;Name=QC\n"
        "c1\tt\tmRNA\t1000\t1300\t.\t+\t.\tID=mC;Parent=gC\n"
        "c1\tt\tCDS\t1000\t1300\t.\t+\t0\tID=cC;Parent=mC\n")
    faa = tmp_path / "x.faa"
    faa.write_text(">PB\nMK\n>QC\nMK\n")     # PA is missing from the FASTA
    for gff_type in ("GFF3", "funannotate", "BRAKER"):
        res = resolve_gff(gff, gff_type, read_fasta_ids(faa))
        got = {r.gene_id: (r.protein_id, r.matched_by) for r in res.rows}
        assert got["gB"] == ("PB", "CDS.protein_id")
        assert got["gA"][0] == "unknown"
        assert got["gC"][0] == "QC"           # weak keys still match in the second pass
