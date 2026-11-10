"""Every directory of the dbcan package is a regular package (has __init__.py).

Without __init__.py, dbcan and its subpackages are namespace packages, and a stale module left by
an old install (dbcan/constants.py from 5.1.x) shadows them: "'dbcan.constants' is not a package".
"""
from pathlib import Path

import dbcan


def test_every_package_directory_has_init():
    root = Path(dbcan.__file__).parent
    missing = [str(d.relative_to(root.parent)) for d in [root, *root.rglob("*")]
               if d.is_dir() and d.name != "__pycache__" and any(d.glob("*.py")) and not (d / "__init__.py").exists()]
    assert missing == []


def test_constants_package_wins_over_stale_module(tmp_path, monkeypatch):
    import importlib, shutil, sys
    pkg = tmp_path / "dbcan"
    shutil.copytree(Path(dbcan.__file__).parent, pkg, ignore=shutil.ignore_patterns("__pycache__"))
    (pkg / "constants.py").write_text("STALE = True\n")
    for name in [m for m in sys.modules if m == "dbcan" or m.startswith("dbcan.")]:
        monkeypatch.delitem(sys.modules, name)
    monkeypatch.syspath_prepend(str(tmp_path))
    mod = importlib.import_module("dbcan.constants.base_constants")
    assert Path(mod.__file__).parent == pkg / "constants"
