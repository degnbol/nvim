"""Fixtures shared by the pytest suites. Run: see README.md."""
from __future__ import annotations

import pathlib
import sys

import pytest


def _loaded_from(module, root):
    """Whether a module's file, or one of a package's directories, is under root."""
    paths = [getattr(module, "__file__", None), *getattr(module, "__path__", [])]
    return any(p and pathlib.Path(p).is_relative_to(root) for p in paths)


@pytest.fixture
def site(tmp_path, monkeypatch):
    """A directory on sys.path; modules imported from it are forgotten afterwards."""
    root = tmp_path / "site"
    root.mkdir()
    monkeypatch.syspath_prepend(str(root))
    yield root
    for name, module in list(sys.modules.items()):
        if _loaded_from(module, root):
            del sys.modules[name]
