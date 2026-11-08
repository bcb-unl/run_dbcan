import pandas as pd

from dbcan.process.process_utils import filter_overlaps


def _overlap_fixture():
    return pd.DataFrame(
        [
            {
                "HMM Name": "GH189.hmm",
                "Target Name": "gt84_allowed",
                "Target From": 1251,
                "Target To": 1529,
                "i-Evalue": 8.7e-133,
            },
            {
                "HMM Name": "GT84.hmm",
                "Target Name": "gt84_allowed",
                "Target From": 1301,
                "Target To": 1515,
                "i-Evalue": 2.8e-104,
            },
            {
                "HMM Name": "GH192.hmm",
                "Target Name": "gt84_competes",
                "Target From": 1251,
                "Target To": 1529,
                "i-Evalue": 8.7e-133,
            },
            {
                "HMM Name": "GT84.hmm",
                "Target Name": "gt84_competes",
                "Target From": 1301,
                "Target To": 1515,
                "i-Evalue": 2.8e-104,
            },
            {
                "HMM Name": "GH13_36.hmm",
                "Target Name": "cbm58_allowed",
                "Target From": 82,
                "Target To": 540,
                "i-Evalue": 2e-109,
            },
            {
                "HMM Name": "CBM58.hmm",
                "Target Name": "cbm58_allowed",
                "Target From": 218,
                "Target To": 334,
                "i-Evalue": 4.3e-66,
            },
        ]
    )


def test_default_overlap_behavior_is_unchanged():
    result = filter_overlaps(_overlap_fixture())
    assert list(zip(result["Target Name"], result["HMM Name"])) == [
        ("cbm58_allowed", "GH13_36.hmm"),
        ("gt84_allowed", "GH189.hmm"),
        ("gt84_competes", "GH192.hmm"),
    ]


def test_pair_specific_coexistence_preserves_only_allowlisted_pairs():
    pairs = {
        frozenset({"GT84", "GH189"}),
        frozenset({"CBM58", "GH13"}),
    }
    result = filter_overlaps(_overlap_fixture(), coexistence_pairs=pairs)
    observed = set(zip(result["Target Name"], result["HMM Name"]))

    assert ("gt84_allowed", "GT84.hmm") in observed
    assert ("gt84_allowed", "GH189.hmm") in observed
    assert ("cbm58_allowed", "CBM58.hmm") in observed
    assert ("cbm58_allowed", "GH13_36.hmm") in observed
    assert ("gt84_competes", "GT84.hmm") not in observed
    assert ("gt84_competes", "GH192.hmm") in observed


def _span_fixture():
    """A long catalytic domain (1-400) and a short accessory domain (300-334) inside its span."""
    return pd.DataFrame(
        [
            {
                "HMM Name": "GH18.hmm",
                "Target Name": "gh18_cbm18",
                "Target From": 1,
                "Target To": 400,
                "i-Evalue": 1e-120,
            },
            {
                "HMM Name": "CBM18.hmm",
                "Target Name": "gh18_cbm18",
                "Target From": 300,
                "Target To": 334,
                "i-Evalue": 1e-30,
            },
        ]
    )


def test_ali_span_overlap_deletes_accessory_domain_inside_the_span():
    # Alignment-span rule (as in 5.2.9): the long domain's span 1-400 covers the
    # short domain 300-334, so the accessory domain is deleted (worse i-Evalue).
    result = filter_overlaps(_span_fixture())
    kept = set(zip(result["Target Name"], result["HMM Name"]))
    assert ("gh18_cbm18", "GH18.hmm") in kept
    assert ("gh18_cbm18", "CBM18.hmm") not in kept
