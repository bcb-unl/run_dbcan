"""The bundled --dbcan-taxon-profile fungi resolves to a shipped, loadable policy."""
import json
from pathlib import Path

import pytest

from dbcan.configs.family_cutoff import load_dbcan_family_cutoff_policy

PROFILE_DIR = Path(__file__).resolve().parent.parent / "dbcan" / "data" / "taxon_profiles"


def test_fungi_profile_is_shipped_and_wellformed():
    cutoff = PROFILE_DIR / "fungi.cutoff.json"
    assert cutoff.exists(), f"bundled fungi profile missing: {cutoff}"
    payload = json.loads(cutoff.read_text())
    assert payload["schema_version"] == 1
    assert payload["deployment_status"] in {"experimental", "approved"}
    assert isinstance(payload["families"], dict) and payload["families"]


def test_config_carries_taxon_profile():
    from dbcan.parameter import create_config
    from dbcan.configs.pyhmmer_config import PyHMMERDBCANConfig

    cfg = create_config(
        PyHMMERDBCANConfig, namespace="dbcan",
        taxon_profile_dbcan="fungi", db_dir="/tmp/db", output_dir="/tmp/out",
    )
    assert cfg.taxon_profile == "fungi"
    # bundle is resolved lazily in the processor, so the raw config is still unset here
    assert cfg.family_cutoff_file is None
