"""On a CGC CAZyme, labels from full-length protein databases (TCDB, prokaryotic TF DIAMOND,
SulfAtlas) are not used in cgc.gff: such a hit labels the whole protein, so it cannot be told
apart from the CAZyme domain (e.g. GH13_31 + TC 8.A.9.1.1: the rBAT/SLC3 'transport accessory'
family of TCDB is itself a GH13 homologue). Domain-level labels (STP, fungal TF HMM, MEROPS
units) stay; whether they sit on the CAZyme domain is decided later by CGCFinder."""
from dbcan.configs.base_config import GFFConfig
from dbcan.IO.gff import get_gff_processor

OVERVIEW = (
    "Gene ID\tEC#\tdbCAN_hmm\tdbCAN_sub\tDIAMOND\t#ofTools\tRecommend Results\tSubstrate\n"
    "g1\t-\tGH13_31(1-400)\tGH13_e1(1-400)\tGH13_31\t3\tGH13_31\t-\n"
    "g2\t-\tGT2(1-300)\t-\tGT2\t2\tGT2\t-\n"
    "g3\t-\tGH5(1-200)\t-\t-\t1\t-\t-\n"
)
HEAD = ("Annotate Name\tAnnotate Length\tTarget Name\tTarget Length\ti-Evalue\tAnnotate From\tAnnotate To\t"
        "Target From\tTarget To\tCoverage\tAnnotate File Name\tType\tSource\n")


def _row(name, target, typ, source):
    return f"{name}\t500\t{target}\t500\t1e-30\t1\t400\t1\t400\t0.8\t{typ}\t{typ}\t{source}\n"


def test_full_length_labels_dropped_on_cazymes_only(tmp_path):
    (tmp_path / "overview.tsv").write_text(OVERVIEW)
    (tmp_path / "total_cgc_info.tsv").write_text(
        HEAD
        + _row("8.A.9.1.1", "g1", "TC", "TC")              # full-length label on a CAZyme -> dropped
        + _row("S1_7", "g1", "Sulfatase", "Sulfatase")     # dropped
        + _row("LacI", "g1", "STP", "STP")                 # domain-level -> kept
        + _row("4.D.1.1.14", "g2", "TC", "TC")             # dropped
        + _row("M38.UNW", "g2", "Peptidase", "Peptidase")  # domain-level -> kept
        + _row("2.A.1.1.1", "g3", "TC", "TC")              # g3 is not a CGC CAZyme (#ofTools 1) -> kept
        + _row("HTH_AraC", "g4", "TF", "TF_prok")          # non-CAZyme -> kept
    )
    proc = get_gff_processor(GFFConfig(output_dir=str(tmp_path), input_gff=str(tmp_path / "x.gff"),
                                       gff_type="prodigal"))
    ann = {k: v["CGC_annotation"] for k, v in proc.load_cgc_type().items()}
    assert ann["g1"] == "CAZyme|GH13_31+STP|LacI"
    assert ann["g2"] == "CAZyme|GT2+Peptidase|M38.UNW"
    assert ann["g3"] == "TC|2.A.1.1.1"
    assert ann["g4"] == "TF|HTH_AraC"


def test_old_signature_table_without_source_keeps_all_labels(tmp_path):
    (tmp_path / "overview.tsv").write_text(OVERVIEW)
    (tmp_path / "total_cgc_info.tsv").write_text(
        HEAD.replace("\tSource", "") + _row("8.A.9.1.1", "g1", "TC", "TC").rsplit("\t", 1)[0] + "\n")
    proc = get_gff_processor(GFFConfig(output_dir=str(tmp_path), input_gff=str(tmp_path / "x.gff"),
                                       gff_type="prodigal"))
    assert proc.load_cgc_type()["g1"]["CGC_annotation"] == "CAZyme|GH13_31+TC|8.A.9.1.1"
