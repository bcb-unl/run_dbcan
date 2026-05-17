from dbcan.utils.utils import GFF_record, dedupe_gff_genes, get_gene_id_from_gff_attributes


def test_protein_id_attribute():
    attr = "protein_id=Gene_123;CGC_annotation=CAZyme|GH5"
    assert get_gene_id_from_gff_attributes(attr) == "Gene_123"


def test_id_attribute_fallback():
    attr = "ID=Gene_456;other=1"
    assert get_gene_id_from_gff_attributes(attr) == "Gene_456"


def test_gff_record_uses_protein_id():
    line = ["ctg1", ".", "gene", "100", "500", ".", "+", ".", "protein_id=P1;CGC_annotation=null"]
    rec = GFF_record(line)
    assert rec.seqid == "P1"


def test_dedupe_merges_spans():
    g1 = GFF_record(["ctg1", ".", "CDS", "100", "200", ".", "+", ".", "protein_id=P1;"])
    g2 = GFF_record(["ctg1", ".", "CDS", "300", "400", ".", "+", ".", "protein_id=P1;"])
    merged = dedupe_gff_genes([g1, g2])
    assert len(merged) == 1
    assert merged[0].start == 100
    assert merged[0].end == 400
