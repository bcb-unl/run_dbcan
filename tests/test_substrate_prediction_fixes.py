"""substrate_prediction: option mapping, PUL-CGC thresholds, IDs containing '|', optional detail tables."""
import logging
from types import SimpleNamespace
from unittest.mock import patch

import pandas as pd
import pytest
from click.testing import CliRunner

from dbcan.annotation.cgc_substrate_prediction import (
    HitParamter,
    QueryID,
    blastp_hit,
    dbCAN_substrate_prediction,
    parse_extra_pairs,
    parse_query_id,
)
from dbcan.configs.cgc_substrate_config import CGCSubstrateConfig
from dbcan.main import cli

CONTIG = "gi|49175990|ref|NC_000913.2|"
PROT = ["gi|16128008|ref|NP_414557.1|", "gi|16128009|ref|NP_414558.1|", "gi|16128010|ref|NP_414559.1|"]
CGC_OUT = (
    "CGC#\tGene Type\tContig ID\tProtein ID\tGene Start\tGene Stop\tGene Strand\tGene Annotation\n"
    f"CGC1\tCAZyme\t{CONTIG}\t{PROT[0]}\t100\t900\t+\tCAZyme|GH5_2\n"
    f"CGC1\tCAZyme\t{CONTIG}\t{PROT[1]}\t1000\t1900\t+\tCAZyme|GH43_10\n"
    f"CGC1\tTC\t{CONTIG}\t{PROT[2]}\t2000\t2900\t+\tTC|2.A.1.1.1\n"
)


def _blast_line(query, subject, bitscore=300.0, pident=60.0):
    return [query, subject, str(pident), "300", "10", "0", "1", "300", "1", "300", "1e-80", str(bitscore), "300", "300"]


@pytest.fixture
def outdir(tmp_path):
    (tmp_path / "cgc_standard_out.tsv").write_text(CGC_OUT)
    (tmp_path / "uniInput.faa").write_text("".join(f">{p}\nMKV\n" for p in PROT))
    return tmp_path


def _sp(outdir, **kw):
    return dbCAN_substrate_prediction(CGCSubstrateConfig(output_dir=str(outdir), db_dir=str(outdir), **kw))


# ---- query IDs with '|' -------------------------------------------------------------------------

def test_parse_query_id_with_pipes_in_contig_and_protein():
    q = parse_query_id(f"{CONTIG}|CGC12|{PROT[0]}|CAZyme|GH5_2")
    assert q == QueryID(CONTIG, "CGC12", PROT[0], "CAZyme", "GH5_2")
    assert q.cgc_id == f"{CONTIG}|CGC12"
    assert parse_query_id("no_cgc_here") is None


def test_parse_query_id_prefers_known_ids():
    weird = QueryID("a|CGC1|b", "CGC2", "p", "TC", "1.A.1")
    qid = "a|CGC1|b|CGC2|p|TC|1.A.1"
    assert parse_query_id(qid, {qid: weird}) is weird


def test_read_blastp_result_groups_by_full_cgc_id(outdir):
    sp = _sp(outdir)
    sp.extract_seq_in_CGC()
    queries = [rec.id for rec in sp.seqs]
    assert queries[0] == f"{CONTIG}|CGC1|{PROT[0]}|CAZyme|GH5_2"
    blast = outdir / "PUL_blast.out"
    blast.write_text("".join(
        "\t".join(_blast_line(q, f"PUL0001_{i}:PUL0001:WP_{i}:CAZyme:GH5_2")) + "\n" for i, q in enumerate(queries)
    ))
    sp.read_blastp_result(str(blast))
    assert list(sp.hitdict) == [f"{CONTIG}|CGC1"]
    assert [h.query.protein for h in sp.hitdict[f"{CONTIG}|CGC1"]] == PROT
    row = sp.hitdict[f"{CONTIG}|CGC1"][0].format_str().split("\t")
    assert row[0] == f"{PROT[0]}|CAZyme|GH5_2"
    assert row[1] == "WP_0|CAZyme|GH5_2"
    assert row[2] == f"{CONTIG}|CGC1"
    assert row[3] == "PUL0001"


# ---- -upghn / -uqcgn ----------------------------------------------------------------------------

def _hits(n_query_genes, n_pul_genes):
    hits = []
    for qi in range(n_query_genes):
        for si in range(n_pul_genes):
            h = blastp_hit(_blast_line(f"c|CGC1|q{qi}|CAZyme|GH5", f"PUL0001_{si}:PUL0001:WP_{si}:CAZyme:GH5"))
            hits.append(h)
    return hits


def test_uniq_pul_gene_and_query_gene_thresholds_are_not_swapped(outdir):
    # one CGC gene hitting two PUL genes: 2 distinct PUL genes, 1 distinct CGC gene
    hits = _hits(n_query_genes=1, n_pul_genes=2)
    common = dict(total_pair_num=1, CAZyme_pair_num=1, bitscore_cutoff=0)
    ok = _sp(outdir, uniq_pul_gene_hit_num=2, uniq_query_cgc_gene_num=1, **common)
    assert ok.Uniq_blastp_hit(hits)[0] > 0
    too_few_cgc_genes = _sp(outdir, uniq_pul_gene_hit_num=1, uniq_query_cgc_gene_num=2, **common)
    assert too_few_cgc_genes.Uniq_blastp_hit(hits)[0] == -1


def test_cazyme_pair_num_is_used(outdir):
    hits = _hits(n_query_genes=2, n_pul_genes=2)  # 4 CAZyme-CAZyme pairs
    common = dict(uniq_pul_gene_hit_num=1, uniq_query_cgc_gene_num=1, total_pair_num=1, bitscore_cutoff=0)
    assert _sp(outdir, CAZyme_pair_num=4, **common).Uniq_blastp_hit(hits)[0] > 0
    assert _sp(outdir, CAZyme_pair_num=5, **common).Uniq_blastp_hit(hits)[0] == -1


# ---- CLI mapping and errors ---------------------------------------------------------------------

@pytest.fixture
def calls():
    recorded = {}
    with patch("dbcan.core.run_dbCAN_CGCFinder_substrate", side_effect=lambda cfg: recorded.setdefault("sub", cfg)), \
         patch("dbcan.core.run_dbcan_syn_plot", side_effect=lambda cfg: recorded.setdefault("syn", cfg)):
        yield recorded


def test_cli_cpn_long_option_reaches_config(tmp_path, calls):
    result = CliRunner().invoke(cli, ["substrate_prediction", "--output_dir", str(tmp_path), "--db_dir", str(tmp_path),
                                      "--CAZyme_pair_num", "4"])
    assert result.exit_code == 0, result.output
    assert calls["sub"].CAZyme_pair_num == 4


def test_cli_cpn_short_option(tmp_path, calls):
    result = CliRunner().invoke(cli, ["substrate_prediction", "--output_dir", str(tmp_path), "--db_dir", str(tmp_path),
                                      "-cpn", "3"])
    assert result.exit_code == 0, result.output
    assert calls["sub"].CAZyme_pair_num == 3


def test_cli_extra_pair_mismatch_exits_nonzero(tmp_path, calls):
    result = CliRunner().invoke(cli, ["substrate_prediction", "--output_dir", str(tmp_path), "--db_dir", str(tmp_path),
                                      "-ept", "TC-TC,STP-STP", "-eptn", "1"])
    assert result.exit_code != 0
    assert "extra_pair_type" in result.output
    assert "sub" not in calls


def test_hitparamter_mismatch_raises():
    cfg = SimpleNamespace(uniq_pul_gene_hit_num=2, uniq_query_cgc_gene_num=2, CAZyme_pair_num=1, total_pair_num=2,
                          extra_pair_type="TC-TC,STP-STP", extra_pair_type_num="1", identity_cutoff=0,
                          coverage_cutoff=0, bitscore_cutoff=50, evalue_cutoff=0.01)
    with pytest.raises(ValueError):
        HitParamter(cfg)
    assert parse_extra_pairs("TC-TC,STP-STP", "1,2") == (["TC-TC", "STP-STP"], ["1", "2"])
    assert parse_extra_pairs(None, "0") == (None, 0)


def test_legacy_substrate_options_are_hidden_accepted_and_warned(tmp_path, calls):
    import dbcan.parameter as P
    with patch.object(P.logger, "warning") as warn:
        result = CliRunner().invoke(cli, ["substrate_prediction", "--output_dir", str(tmp_path), "--db_dir",
                                          str(tmp_path), "-o", "x.out", "-w", "/tmp", "-rerun", "true", "-env", "server"])
    assert result.exit_code == 0, result.output
    warned = {c.args[1] for c in warn.call_args_list}
    assert warned == {"--out", "--workdir", "--rerun", "--env"}
    helptext = CliRunner().invoke(cli, ["substrate_prediction", "--help"]).output
    for opt in ("--workdir", "--rerun", "--env", "--pul "):
        assert opt not in helptext
    assert "--odbcan_sub" in helptext and "--odbcanpul" in helptext


# ---- PUL table guard ----------------------------------------------------------------------------

def test_missing_pul_id_in_excel_is_a_warning(outdir, caplog):
    sp = _sp(outdir)
    sp.Puls = pd.DataFrame({"ID": ["PUL0001"], "substrate_final": ["xylan"]})
    assert sp.pul_substrate("PUL0001") == "xylan"
    with caplog.at_level(logging.WARNING):
        assert sp.pul_substrate("PUL9999") == "-"
    assert "PUL9999" in caplog.text


# ---- --odbcan_sub detail table ------------------------------------------------------------------

def test_dbcan_sub_detail_table(outdir):
    cols = ["Subfam Name", "Subfam Composition", "Subfam EC", "Substrate", "HMM Length", "Target Name",
            "Target Length", "i-Evalue", "HMM From", "HMM To", "Target From", "Target To", "Coverage", "HMM File Name"]
    rows = [
        ["GH5_e1", "GH5_2:3", "3.2.1.4:3", "cellulose", 300, PROT[0], 400, 1e-50, 1, 290, 10, 300, 0.96, "dbCAN-sub"],
        ["GH43_e2", "GH43_10:3", "3.2.1.37:2", "xylan", 300, PROT[1], 400, 1e-50, 1, 290, 10, 300, 0.96, "dbCAN-sub"],
        ["GH43_e3", "GH43_10:3", "3.2.1.37:2", "xylan", 300, PROT[1], 400, 1e-40, 1, 290, 310, 390, 0.9, "dbCAN-sub"],
    ]
    pd.DataFrame(rows, columns=cols).to_csv(outdir / "dbCANsub_hmm_results.tsv", sep="\t", index=False)
    sp = _sp(outdir, odbcan_sub=True, num_of_domains_substrate_cutoff=2, num_of_protein_substrate_cutoff=2,
             substrate_scors=2)
    sp.dbcan_sub_subfamily_substrate_prediction()
    assert sp.dbcan_sub_CGC2substrates == {f"{CONTIG}|CGC1": "xylan"}
    df = pd.read_csv(sp.dbcan_sub_intermediate_file(), sep="\t", dtype=str)
    assert len(df) == 3
    assert set(df["#cgcid"]) == {f"{CONTIG}|CGC1"}
    assert df["CGC substrate scores"].iloc[0] == "xylan:2;cellulose:1"
    assert set(df["CGC dbCAN-sub substrate"]) == {"xylan"}
    assert set(df["CGC dbCAN-sub substrate score"]) == {"2"}


# ---- multi-role annotations, row order, shared query-ID parser ------------------------------------

@pytest.mark.parametrize("annotation,gene_type,family", [
    ("CAZyme|GH32_e132+Sulfatase|S1_7", "CAZyme", "GH32_e132"),
    ("CAZyme|GH13_22|GT5_e25", "CAZyme", "GH13_22+GT5_e25"),
    ("STP|HTH_1+STP|LysR_substrate", "STP", "HTH_1+LysR_substrate"),
    ("TF|OXYR_ECOLI+STP|HTH_1+STP|LysR_substrate", "TF", "OXYR_ECOLI"),
    ("TC|1.I.1.1.3+STP|Pkinase", "TC", "1.I.1"),
    ("CAZyme|CBM50_e1400+TC|1.A.34.1.2+Peptidase|M23.UNB+STP|Peptidase_M23", "CAZyme", "CBM50_e1400"),
    ("Other|tRNA", "Other", "tRNA"),
    ("null", "null", "null"),
])
def test_gene_family_uses_the_gene_type_role(annotation, gene_type, family):
    from dbcan.annotation.cgc_substrate_prediction import Gene
    g = Gene(["CGC1", gene_type, CONTIG, PROT[0], "1", "90", "+", annotation])
    assert g.Protein_Fam == family
    qid = g.Get_CAzyID()
    assert qid == f"{CONTIG}|CGC1|{PROT[0]}|{gene_type}|{family}"
    # the query ID reads back with the shared parser, also without the known-ID map
    assert parse_query_id(qid) == QueryID(CONTIG, "CGC1", PROT[0], gene_type, family)


def test_substrate_prediction_rows_sorted_by_contig_then_cgc_number(outdir):
    from dbcan.annotation.cgc_substrate_prediction import PULhit, cgc_sort_key
    sp = _sp(outdir)
    ids = ["b|CGC10", "a|CGC2", "b|CGC9", f"{CONTIG}|CGC3", "a|CGC10", f"{CONTIG}|CGC1"]
    sp.queryCGC2hit = {c: PULhit(100.0, "PUL0001", "xylan", ["CAZyme-CAZyme"]) for c in ids[:3]}
    sp.dbcan_sub_CGC2substrates = {c: "xylan" for c in ids[3:]}
    sp.result_print_to_file()
    rows = [l.split("\t")[0] for l in open(sp.out).read().splitlines()[1:]]
    assert rows == ["a|CGC2", "a|CGC10", "b|CGC9", "b|CGC10", f"{CONTIG}|CGC1", f"{CONTIG}|CGC3"]
    assert sorted(ids, key=cgc_sort_key) == rows


@pytest.mark.parametrize("contig", [CONTIG, "scaffold|CGC7|v2"])
def test_synteny_links_use_the_substrate_query_ids(outdir, tmp_path, contig):
    """PUL_blast.out written from the CGC.faa IDs of substrate_prediction: every hit becomes a link,
    with '|' in contig and protein IDs (even a contig containing a CGC token) and a multi-role gene."""
    from dbcan.plot import cgc_viz as V
    from dbcan.plot.syntenic_plot import HomologyCutoffs, read_pul_links
    cgc = outdir / "cgc_standard_out.tsv"
    cgc.write_text(CGC_OUT.replace("CAZyme|GH43_10", "CAZyme|GH43_10+Sulfatase|S1_7").replace(CONTIG, contig))
    sp = _sp(outdir)
    sp.extract_seq_in_CGC()
    queries = [rec.id for rec in sp.seqs]
    assert queries[1] == f"{contig}|CGC1|{PROT[1]}|CAZyme|GH43_10"
    blast = outdir / "PUL_blast.out"
    blast.write_text("".join("\t".join(_blast_line(q, f"PUL0001_{i}:PUL0001:WP_{i}:CAZyme:GH5_2")) + "\n"
                             for i, q in enumerate(queries)))
    cgcs = V.read_cgc_table(cgc)
    cgcid = f"{contig}|CGC1"
    assert cgcid in cgcs
    pul_genes = [V.Gene.build(f"WP_{i}", "PUL0001", 1 + 1000 * i, 900 + 1000 * i, "+", "CAZyme|GH5_2")
                 for i in range(3)]
    links = read_pul_links(blast, cgcid, "PUL0001", cgcs[cgcid], pul_genes, HomologyCutoffs(0.01, 0, 0))
    assert sorted((l.q.pid, l.s.pid) for l in links) == [(PROT[i], f"WP_{i}") for i in range(3)]
