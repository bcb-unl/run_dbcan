"""dbCAN-sub EC / substrate parsing and the (family, EC) -> substrate rules."""
from types import SimpleNamespace

import pytest

from dbcan.process.process_dbcan_sub import DBCANSUBProcessor as P

MAPPING = (
    "Substrate_high_level\t Substrate_curated\tFamily\tName\tEC_Number\n"
    "alpha-glucan\talpha-1,3-glucan\tGH71\tmutanase\t3.2.1.59\n"
    "lignin\tlignin\tAA1\tlaccase\t1.10.3.2\n"
    "beta-glucan;beta-mannan\tx\tGH5\tendo\t3.2.1.4\n"
    "xylan\txylan\tGH10\t\"endo-1,4-xylanase\"\t\"\t3.2.1.8\"\n"
    "chitin\tchitin\tCBM5\tchitin binding\t\n"
    "lignin\tlignin\tAA8\tiron reductase\tNA\n"
    "xylan\txylan\tAA7\toxidase\t1.1.3.-\n"
)


@pytest.fixture
def proc(tmp_path):
    (tmp_path / "fam-substrate-mapping.tsv").write_text(MAPPING)
    p = P(SimpleNamespace(output_dir=str(tmp_path), db_dir=str(tmp_path)))
    return p, p.load_substrate_mapping()


@pytest.mark.parametrize(
    "name, expected_ec, expected_sub",
    [
        # no EC: the family's only substrate is assigned (GH71 -> alpha-glucan only)
        ("GH71_e0.hmm|GH71:56", "-", "alpha-glucan"),
        ("GH71_e35.hmm|CBM0:2|GH71:29", "-", "alpha-glucan"),
        # exact (family, EC) match
        ("GH71_e11.hmm|CBM24:30|GH71:83|3.2.1.59:2", "3.2.1.59:2", "alpha-glucan"),
        # EC token without ':count' is still an EC
        ("AA1_e1.hmm|AA1:57|AA1_3:52|1.10.3.2", "1.10.3.2:1", "lignin"),
        # complete EC not in the mapping -> family fallback (AA1 has a single substrate)
        ("AA1_e9.hmm|AA1:5|1.10.3.3:2", "1.10.3.3:2", "lignin"),
        # family-level mapping rows (empty / NA EC)
        ("CBM5_e1.hmm|CBM5:10", "-", "chitin"),
        ("AA8_e0.hmm|AA8:22", "-", "lignin"),
        # incomplete EC matches only its exact key
        ("AA7_e1.hmm|AA7:42|1.1.3.-:3", "1.1.3.-:3", "xylan"),
        # several ECs: one substrate entry per EC, same order as Subfam EC
        ("AA1_e2.hmm|AA1:9|1.10.3.3:1|1.10.3.2", "1.10.3.2:1;1.10.3.3:1", "lignin;lignin"),
        # incomplete EC without an exact key -> '-' (no family fallback for a.b.c.-)
        ("AA1_e3.hmm|AA1:9|1.10.3.-:1", "1.10.3.-:1", "-"),
        # a mapping value listing several substrates stays inside its own entry
        ("GH5_e3.hmm|GH5:4|3.2.1.4:3|3.2.1.8:1", "3.2.1.4:3;3.2.1.8:1",
         "beta-glucan,beta-mannan;beta-glucan,beta-mannan"),
    ],
)
def test_strict_substrate_mapping(proc, name, expected_ec, expected_sub):
    p, subs = proc
    assert p._extract_subfamily_ecs(name) == expected_ec
    assert p.get_substrates(name, subs) == expected_sub


def test_mapping_ec_whitespace_is_stripped(proc):
    p, subs = proc
    # the GH10 row has a quoted EC with a stray tab, as in the released mapping file
    assert p.get_substrates("GH10_e1.hmm|GH10:3|3.2.1.8:2", subs) == "xylan"


def test_ec_tokens_with_wildcards_and_preliminary_numbers():
    name = "GT1_e109.hmm|GT1:43|2.4.1.-:3|2.4.-.-:2|2.4.1.n1|2.4.1.17"
    assert P._parse_sorted_ec_tokens(name) == ["2.4.1.17", "2.4.1.n1", "2.4.1.-:3", "2.4.-.-:2"]
    assert P._extract_subfamily_components(name) == "GT1:43"


def test_ec_lookup_uses_the_catalytic_family_of_the_composition(tmp_path):
    """CBM1_e3 (CBM1:49|GH45:49|3.2.1.4): 3.2.1.4 is a GH45 activity, not a CBM1 one."""
    (tmp_path / "fam-substrate-mapping.tsv").write_text(
        "Substrate_high_level\tSubstrate_curated\tFamily\tName\tEC_Number\n"
        "cellulose;chitin\tx\tCBM1\tcellulose/chitin binding\t\n"
        "cellulose\tx\tGH45\tendoglucanase\t3.2.1.4\n"
    )
    p = P(SimpleNamespace(output_dir=str(tmp_path), db_dir=str(tmp_path)))
    subs = p.load_substrate_mapping()
    assert p.get_substrates("CBM1_e3.hmm|CBM1:49|GH45:49|3.2.1.4:18", subs) == "cellulose"
    # no EC: family-level row of the eHMM's own family
    assert p.get_substrates("CBM1_e9.hmm|CBM1:5", subs) == "cellulose,chitin"


def test_web_release_mapping_format_is_read(tmp_path):
    """fam-substrate-mapping-08262025.tsv: comment line, new_Substrate_high_level, type, blank cells."""
    (tmp_path / "fam-substrate-mapping.tsv").write_text(
        "#substrate mapping table manually curated (updated 8/26/2025)\t\t\t\t\t\t\n"
        "Substrate_high_level\t Substrate_curated\tFamily\tName\tEC_Number\tnew_Substrate_high_level\ttype\n"
        "lignin\tlignin\tAA1\tferroxidase\t1.10.3.-\t\tREMOVED\n"
        "\t\tAA3\tEcdysone oxidase\t1.1.3.16\tlignin\tNEW\n"
        "cellulose\tcellulose\tAA9\tlytic xylan monooxygenase\t1.14.99.-\txylan\tCHANGED\n"
        "\t\tGH3\tsomething\t2.4.1.-\t\t\n"
        "alpha-glucan\talpha-glucan\tGH71\tmutanase\t3.2.1.59\talpha-glucan\t\n"
    )
    subs = P(SimpleNamespace(output_dir=str(tmp_path), db_dir=str(tmp_path))).load_substrate_mapping()
    assert subs == {
        ("AA3", "1.1.3.16"): ["lignin"],
        ("AA9", "1.14.99.-"): ["xylan"],
        ("GH71", "3.2.1.59"): ["alpha-glucan"],
    }


def test_per_subfamily_table_from_the_web_server_takes_precedence(tmp_path):
    (tmp_path / "fam-substrate-mapping.tsv").write_text(
        "Substrate_high_level\tSubstrate_curated\tFamily\tName\tEC_Number\n"
        "cellulose;chitin\tx\tCBM1\tbinding\t\n"
    )
    (tmp_path / "dbCAN-sub-substrate.tsv").write_text(
        "Subfamily\tEC\tCount\tSubstrate\n"
        "CBM1_e1\t3.2.1.4\t3\tcellulose\n"
        "CBM1_e1\t3.2.1.132\t4\tchitosan\n"
        "CBM1_e1\t3.2.1.176\t49\tcellulose\n"
        "GT1_e9\t2.4.1.-\t2\t\n"
    )
    p = P(SimpleNamespace(output_dir=str(tmp_path), db_dir=str(tmp_path)))
    subs, table = p.load_substrate_mapping(), p.load_subfamily_substrate_table()
    name = "CBM1_e1.hmm|CBM1:217|GH7:184|3.2.1.4:3|3.2.1.132:4|3.2.1.176:49"
    assert p.get_substrates(name, subs, table) == "cellulose;chitosan;cellulose"
    assert p.get_substrates("GT1_e9.hmm|GT1:30|2.4.1.-:2", subs, table) == "-"
    # without the table: family-level fallback of CBM1 for every EC
    assert p.get_substrates(name, subs) == "cellulose,chitin;cellulose,chitin;cellulose,chitin"


def test_mapping_with_version_comment_and_extra_columns(tmp_path):
    """A versioned mapping file: a '#' provenance line without tabs and extra columns
    (Substrate_specific, Source, Change) after the 5 classic ones."""
    (tmp_path / "fam-substrate-mapping.tsv").write_text(
        "# fam-substrate-mapping version 2026-11-09 | merged from db_current and dbCAN-sub web 2025-08-26\n"
        "Substrate_high_level\tSubstrate_specific\tFamily\tName\tEC_Number\tSource\tChange\n"
        "alpha-glucan\tstarch\tGH13\tneopullulanase\t3.2.1.135\tboth\tspecific name from web\n"
        "chitin;peptidoglycan\tchitin;peptidoglycan\tCBM50\tLysM\t\tboth\t\n"
    )
    p = P(SimpleNamespace(output_dir=str(tmp_path), db_dir=str(tmp_path)))
    subs = p.load_substrate_mapping()
    assert subs[("GH13", "3.2.1.135")] == ["alpha-glucan"]
    assert p.get_substrates("GH13_e1.hmm|GH13:5|3.2.1.135:2", subs) == "alpha-glucan"
