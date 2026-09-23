"""Tests for the patch sequence. Run: see README.md. The rdkit test runs only
where rdkit is importable."""
from __future__ import annotations

import ast
import importlib
import math
import os
import pathlib
import re
import shutil
import subprocess
import sys

import libcst as cst
import patch_pybind_stubs
import patch_stubs
import pytest
from stub_tree import comment_keyword_targets

RUNTIME = "from math import sqrt\n\ndef py_func(x):\n    return x\n"
STUB = "def sqrt(x: float) -> float: ...\n"


@pytest.fixture
def site(tmp_path, monkeypatch):
    """A directory on sys.path; modules imported from it are forgotten afterwards."""
    root = tmp_path / "site"
    root.mkdir()
    monkeypatch.syspath_prepend(str(root))
    before = set(sys.modules)
    yield root
    for name in set(sys.modules) - before:
        del sys.modules[name]


def write(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    importlib.invalidate_caches()
    return path


def snapshot(root):
    """Bytes of every file under a directory, keyed by relative path."""
    return {p.relative_to(root): p.read_bytes() for p in root.rglob("*") if p.is_file()}


def head(path):
    return path.read_text().split("\n", 1)[0]


@pytest.fixture
def tree(site, tmp_path):
    """A `<package>-stubs` tree over an importable runtime package."""
    write(site / "fakepkg" / "__init__.py", RUNTIME)
    write(site / "fakepkg" / "other.py", "x = 1\n")
    root = tmp_path / "stubs" / "fakepkg-stubs"
    write(root / "__init__.pyi", STUB)
    write(root / "other.pyi", "x: int\n")
    return root


def test_documents_a_stubs_root(tree):
    assert patch_stubs.patch(tree) == {}
    text = (tree / "__init__.pyi").read_text()
    assert head(tree / "__init__.pyi").startswith("# patch_stubs: fakepkg unknown ")
    assert str(math.sqrt.__doc__).split("\n")[0] in text


def test_documents_a_single_module_stub(site, tmp_path):
    write(site / "fakemod.py", RUNTIME)
    stub = write(tmp_path / "stubs" / "fakemod.pyi", STUB)
    assert patch_stubs.patch(stub) == {}
    assert head(stub).startswith("# patch_stubs: fakemod ")
    assert ast.get_docstring(ast.parse(stub.read_text()).body[0]) == math.sqrt.__doc__


def test_a_package_directory_keeps_its_other_files(site):
    package = site / "fakepkg"
    write(package / "__init__.py", RUNTIME)
    write(package / "__init__.pyi", STUB)
    write(package / "native.so", "\x7fELF")
    before = snapshot(package)
    assert patch_stubs.patch(package) == {}
    after = snapshot(package)
    assert after.pop(pathlib.Path("__init__.pyi")) != before.pop(pathlib.Path("__init__.pyi"))
    assert {k: v for k, v in after.items() if k.parts[0] != "__pycache__"} == before


def test_leaves_unchanged_files_and_hard_links_alone(tree, tmp_path):
    other = tree / "other.pyi"
    inode = other.stat().st_ino
    link = tmp_path / "link.pyi"
    os.link(tree / "__init__.pyi", link)
    patch_stubs.patch(tree)
    assert other.stat().st_ino == inode
    assert link.read_text() == STUB


def test_marker_is_on_the_marker_file_only(tree):
    patch_stubs.patch(tree)
    assert head(tree / "__init__.pyi").startswith("# patch_stubs:")
    assert (tree / "other.pyi").read_text() == "x: int\n"


def test_removes_the_legacy_marker(tree):
    for path in tree.rglob("*.pyi"):
        path.write_text("# fix_pybind_stubs: fakepkg 1.0 deadbeef\n" + path.read_text())
    patch_stubs.patch(tree)
    assert not any("fix_pybind_stubs" in p.read_text() for p in tree.rglob("*.pyi"))


def test_second_run_changes_nothing(tree):
    patch_stubs.patch(tree)
    once = snapshot(tree)
    assert patch_stubs.patch(tree) is None
    assert snapshot(tree) == once


def test_nothing_to_document_writes_only_the_marker(site, tmp_path):
    write(site / "fakepkg" / "__init__.py", RUNTIME)
    root = write(tmp_path / "fakepkg-stubs" / "__init__.pyi", "def py_func(x): ...\n").parent
    assert patch_stubs.patch(root) is None
    assert (root / "__init__.pyi").read_text().split("\n", 1)[1] == "def py_func(x): ...\n"
    assert head(root / "__init__.pyi").startswith("# patch_stubs:")


def error_types(skipped):
    return {module: type(e) for module, e in skipped.items()}


def test_a_stub_that_does_not_parse_is_skipped_and_the_tree_marked(tree):
    write(tree / "other.pyi", "def f(a=1, b): ...\n")
    assert error_types(patch_stubs.patch(tree)) == {"fakepkg.other": cst.ParserSyntaxError}
    assert head(tree / "__init__.pyi").startswith("# patch_stubs:")
    assert str(math.sqrt.__doc__).split("\n")[0] in (tree / "__init__.pyi").read_text()


def test_a_failed_import_is_reported_when_nothing_else_changed(site, tmp_path):
    write(site / "fakepkg" / "__init__.py", RUNTIME)
    write(site / "fakepkg" / "broken.py", "raise ImportError('no')\n")
    root = write(tmp_path / "fakepkg-stubs" / "__init__.pyi", "def py_func(x): ...\n").parent
    write(root / "broken.pyi", "x: int\n")
    assert error_types(patch_stubs.patch(root)) == {"fakepkg.broken": ImportError}


def test_a_package_that_does_not_import_fails_unmarked(site, tree):
    write(site / "fakepkg" / "__init__.py", "raise ImportError('broken')\n")
    before = snapshot(tree)
    for _ in range(2):
        with pytest.raises(SystemExit, match="fakepkg does not import: broken"):
            patch_stubs.patch(tree)
        assert snapshot(tree) == before


def dist_info(site, version):
    """Declare an installed fakepkg distribution of a version, replacing any other."""
    for old in site.glob("fakepkg-*.dist-info"):
        shutil.rmtree(old)
    info = site / f"fakepkg-{version}.dist-info"
    write(info / "METADATA", f"Metadata-Version: 2.1\nName: fakepkg\nVersion: {version}\n")
    write(info / "top_level.txt", "fakepkg\n")


def test_a_version_change_reruns(site, tree):
    dist_info(site, "1.0")
    patch_stubs.patch(tree)
    assert head(tree / "__init__.pyi").startswith("# patch_stubs: fakepkg 1.0 ")
    dist_info(site, "2.0")
    patch_stubs.patch(tree)
    assert head(tree / "__init__.pyi").startswith("# patch_stubs: fakepkg 2.0 ")


def test_a_stub_only_distribution_is_left_alone(tmp_path):
    root = write(tmp_path / "nosuchpkg-stubs" / "__init__.pyi", STUB).parent
    assert patch_stubs.patch(root) is None
    assert (root / "__init__.pyi").read_text() == STUB


def test_a_tree_in_the_uv_archive_is_left_alone(site, tmp_path, monkeypatch):
    write(site / "fakepkg" / "__init__.py", RUNTIME)
    cache = tmp_path / "cache"
    monkeypatch.setenv("UV_CACHE_DIR", str(cache))
    root = write(cache / "archive-v0" / "abc" / "fakepkg-stubs" / "__init__.pyi", STUB).parent
    assert patch_stubs.patch(root) is None
    assert (root / "__init__.pyi").read_text() == STUB


def test_removes_stale_temp_files(tree):
    stale = write(tree / "sub" / ".patch_stubs-abc.tmp", "")
    patch_stubs.patch(tree)
    assert not stale.exists()


def test_sources_skip_tests_and_dotfiles(tmp_path):
    for name in ("requirements.txt", "a.py", ".#a.py", "test_a.py", "notes.md"):
        write(tmp_path / name, "")
    assert [p.name for p in patch_stubs.sources(tmp_path)] == ["a.py", "requirements.txt"]


def run_main(site, target):
    """patch_stubs.py run on a target, with site importable."""
    return subprocess.run([sys.executable, str(patch_stubs.__file__), target],
                          env={**os.environ, "PYTHONPATH": str(site)},
                          capture_output=True, text=True)


def test_main_exit_statuses(site, tree):
    assert run_main(site, tree).returncode == 0
    assert run_main(site, tree).returncode == patch_stubs.NOTHING_WRITTEN
    failed = run_main(site, tree.parent / "missing-stubs")
    assert failed.returncode not in (0, patch_stubs.NOTHING_WRITTEN)
    assert "missing-stubs" in failed.stderr


def test_main_prints_one_line_per_skipped_module(site, tree):
    write(site / "fakepkg" / "other.py", "raise ImportError('first\\nsecond')\n")
    write(site / "fakepkg" / "exits.py", "raise SystemExit\n")
    write(tree / "exits.pyi", "x: int\n")
    assert run_main(site, tree).stdout == (
        "fakepkg.exits (SystemExit)\nfakepkg.other (ImportError: first)\n")


def raised(parse, text):
    with pytest.raises(Exception) as caught:
        parse(text)
    return caught.value


def test_a_parser_error_is_described_by_its_message_and_line():
    syntax = raised(ast.parse, "x = 1\ndef f(a=1, b): ...\n")
    assert patch_stubs.describe(syntax) == f"SyntaxError: {syntax.msg} (line 2)"
    token = raised(comment_keyword_targets, 'x = 1\nx = """abc\n')
    assert patch_stubs.describe(token) == "TokenError: EOF in multi-line string (line 2)"
    libcst = raised(cst.parse_module, "x = 1\ndef f(a=1, b): ...\n")
    assert patch_stubs.describe(libcst) == "ParserSyntaxError: expected one of :, = (line 2)"


def test_any_other_error_is_described_by_its_first_line():
    assert patch_stubs.describe(ImportError("first\nsecond")) == "ImportError: first"
    assert patch_stubs.describe(SystemExit()) == "SystemExit"


def test_main_stdout_holds_nothing_compiled_code_prints(site, tree):
    # C's printf buffers when stdout is a pipe, so it would flush at exit.
    write(site / "fakepkg" / "other.py", "import ctypes, os\n"
          "os.write(1, b'raw\\n')\nctypes.CDLL(None).printf(b'buffered\\n')\n"
          "raise ImportError('no')\n")
    assert run_main(site, tree).stdout == "fakepkg.other (ImportError: no)\n"


_TYPED_PROPERTY = re.compile(r"@property\n    def (\w+)\(self\) -> (?:bool|int|float|str):")


def parsed(path):
    """A file's syntax tree, or None where it does not parse."""
    try:
        return ast.parse(path.read_text())
    except SyntaxError:
        return None


def typed_properties(root):
    return {(p.relative_to(root), name) for p in root.rglob("*.pyi")
            for name in _TYPED_PROPERTY.findall(p.read_text())}


def test_rdkit_sequence_matches_the_repair_alone(tmp_path):
    pytest.importorskip("rdkit")
    bundled = patch_pybind_stubs.stubs_path("rdkit")
    repaired = shutil.copytree(bundled, tmp_path / "repaired")
    patch_pybind_stubs.patch_dir(repaired, "rdkit")
    sequenced = shutil.copytree(bundled, tmp_path / "rdkit-stubs")
    patch_stubs.patch(sequenced)

    trees = {path: parsed(path) for path in sequenced.rglob("*.pyi")}
    # Some of rdkit's own stubs do not parse (a non-default parameter after a
    # default one); the sequence must not add to them.
    unparsed = {path.relative_to(sequenced) for path, tree in trees.items() if tree is None}
    assert unparsed == {path.relative_to(repaired) for path in repaired.rglob("*.pyi")
                        if parsed(path) is None}
    constructor_docs = [
        path for path, tree in trees.items() if tree
        for node in ast.walk(tree)
        if isinstance(node, ast.FunctionDef) and node.name in ("__init__", "__new__")
        and "C++ signature" in (ast.get_docstring(node) or "")]
    assert constructor_docs == []
    assert typed_properties(sequenced) >= typed_properties(repaired)
    once = snapshot(sequenced)
    assert patch_stubs.patch(sequenced) is None
    assert snapshot(sequenced) == once
