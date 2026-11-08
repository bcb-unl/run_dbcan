"""Assembly-free abundance (dbcan_asmfree): CAZy header parsing and reads file types."""
import gzip

import click
import pytest

from dbcan.utils.diamond_unassembly import (
    PafRecord,
    _process_paf_file_streaming,
    check_read_type,
    extract_CAZy_from_Tsn,
    get_count_reads,
    seqid_from_Tsn,
)

JGI = "Krezon1_GeneCatalog_proteins_20220127.aa.fasta"


@pytest.mark.parametrize(
    "tsn, families, seqid",
    [
        ("CAD7272370.1|GH33", ["GH33"], "CAD7272370.1"),
        ("AAK49401.1|GH19", ["GH19"], "AAK49401.1"),
        ("1803517A|GH13_6|GH13_6|3.2.1.1", ["GH13_6"], "1803517A"),
        (f"CE16|390012|{JGI}", ["CE16"], f"390012|{JGI}"),
        (f"GH5_7|CBM1|390012|{JGI}", ["GH5_7", "CBM1"], f"390012|{JGI}"),
        (f"GH43_30|GH43_30|1033159|{JGI}", ["GH43_30"], f"1033159|{JGI}"),
        (f"390012|GH5|{JGI}", ["GH5"], "390012"),
        (f"GH5|390012|X.aa.faa", ["GH5"], "390012|X.aa.faa"),
        ("nothing|here", [], "nothing"),
    ],
)
def test_family_tokens_and_seqid(tsn, families, seqid):
    assert extract_CAZy_from_Tsn(tsn) == families
    assert seqid_from_Tsn(tsn) == seqid


def _paf_line(read, tsn, slen=300):
    # qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen
    return "\t".join([read, tsn, "90", "50", "0", "0", "1", "150", "1", "50", "1e-10", "80", "150", str(slen)])


def test_jgi_sequences_of_one_family_are_counted_separately(tmp_path):
    paf = tmp_path / "r.blastx"
    paf.write_text("\n".join([
        _paf_line("r1", f"CE16|390012|{JGI}"),
        _paf_line("r2", f"CE16|390013|{JGI}"),
        _paf_line("r3", f"CE16|390013|{JGI}"),
    ]) + "\n")
    res = _process_paf_file_streaming(str(paf))
    assert res["cazy2seqid"] == {"CE16": {f"390012|{JGI}", f"390013|{JGI}"}}
    assert {k: len(v) for k, v in res["seqid2readid"].items()} == {f"390012|{JGI}": 1, f"390013|{JGI}": 2}
    rec = PafRecord(_paf_line("r1", f"GH5_7|CBM1|390012|{JGI}").split("\t"))
    assert rec.CAZys == {"GH5_7", "CBM1"}


@pytest.mark.parametrize("name, kind", [
    ("x.fq", "fq"), ("x.fastq", "fq"), ("x.fq.gz", "fq"), ("x.FASTQ.GZ", "fq"),
    ("x.fa", "fa"), ("x.fasta", "fa"), ("x.fna.gz", "fa"), ("x.fasta.gz", "fa"),
])
def test_read_file_types(name, kind):
    assert check_read_type(name) == kind


def test_unsupported_read_file_type():
    with pytest.raises(click.BadParameter):
        check_read_type("x.bam")


def test_count_reads_fastq_gz_and_fasta(tmp_path):
    fq = tmp_path / "r.fastq.gz"
    with gzip.open(fq, "wt") as fh:
        fh.write("@a\nACGT\n+\nIIII\n@b\nACGT\n+\nIIII\n")
    fa = tmp_path / "r.fasta"
    fa.write_text(">a\nACGT\n>b\nACGT\n>c\nA\n")
    assert get_count_reads(str(fq)) == 2.0
    assert get_count_reads(str(fa)) == 3.0
