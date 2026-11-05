"""dbCAN-sub EC / substrate parsing (GitHub issue #85, strict substrate mode)."""
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
        # #85: no EC in the eHMM -> no substrate, even though GH71 has a single substrate
        ("GH71_e0.hmm|GH71:56", "-", "-"),
        ("GH71_e35.hmm|CBM0:2|GH71:29", "-", "-"),
        # exact (family, EC) match
        ("GH71_e11.hmm|CBM24:30|GH71:83|3.2.1.59:2", "3.2.1.59:2", "alpha-glucan"),
        # EC token without ':count' is still an EC
        ("AA1_e1.hmm|AA1:57|AA1_3:52|1.10.3.2", "1.10.3.2", "lignin"),
        # complete EC not in the mapping -> '-', no family fallback
        ("AA1_e9.hmm|AA1:5|1.10.3.3:2", "1.10.3.3:2", "-"),
        # family-level mapping rows (empty / NA EC) are not used in strict mode
        ("CBM5_e1.hmm|CBM5:10", "-", "-"),
        ("AA8_e0.hmm|AA8:22", "-", "-"),
        # incomplete EC matches only its exact key
        ("AA7_e1.hmm|AA7:42|1.1.3.-:3", "1.1.3.-:3", "xylan"),
        # several ECs: one substrate entry per EC, same order as Subfam EC
        ("AA1_e2.hmm|AA1:9|1.10.3.3:1|1.10.3.2", "1.10.3.2;1.10.3.3:1", "lignin;-"),
        # a mapping value listing several substrates stays inside its own entry
        ("GH5_e3.hmm|GH5:4|3.2.1.4:3|3.2.1.8:1", "3.2.1.4:3;3.2.1.8:1", "beta-glucan,beta-mannan;-"),
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
