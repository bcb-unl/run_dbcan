"""Gene IDs of the GFF used by cal_coverage (dbcan_utils and run_dbcan expression).

The GFF_record / dedupe_gff_genes helpers were replaced by read_gene_regions(); these
tests keep the original expectations: protein_id wins over ID, ID is the fallback, and
rows sharing one protein ID are one gene.
"""
from dbcan.utils.utils import read_gene_regions


def _gff(tmp_path, rows):
    path = tmp_path / "genes.gff"
    path.write_text("".join("\t".join(map(str, r)) + "\n" for r in rows))
    return str(path)


def test_protein_id_attribute(tmp_path):
    gff = _gff(tmp_path, [["ctg1", ".", "CDS", 100, 500, ".", "+", "0", "ID=cds1;protein_id=Gene_123;CGC_annotation=CAZyme|GH5"]])
    genes = read_gene_regions(gff)
    assert [g.seqid for g in genes] == ["Gene_123"]


def test_id_attribute_fallback(tmp_path):
    gff = _gff(tmp_path, [["ctg1", ".", "CDS", 100, 500, ".", "+", "0", "ID=Gene_456;other=1"]])
    assert [g.seqid for g in read_gene_regions(gff)] == ["Gene_456"]


def test_cgc_gff_gene_rows_use_protein_id(tmp_path):
    gff = _gff(tmp_path, [["ctg1", ".", "gene", 100, 500, ".", "+", ".", "protein_id=P1;CGC_annotation=null"]])
    assert [g.seqid for g in read_gene_regions(gff)] == ["P1"]


def test_rows_of_one_protein_are_one_gene(tmp_path):
    gff = _gff(tmp_path, [
        ["ctg1", ".", "CDS", 100, 200, ".", "+", "0", "protein_id=P1"],
        ["ctg1", ".", "CDS", 300, 400, ".", "+", "0", "protein_id=P1"],
    ])
    genes = read_gene_regions(gff)
    assert len(genes) == 1
    g = genes[0]
    # exons are kept as separate intervals: the intron is not part of the gene length
    assert g.intervals == [("ctg1", 99, 200), ("ctg1", 299, 400)]
    assert g.length == 202
