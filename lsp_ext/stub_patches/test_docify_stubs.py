"""Tests for the docstring pass. Run: see README.md. Needs Python ≥ 3.12, whose
``ast`` parses the PEP 695 fixture."""
from __future__ import annotations

import ast
import importlib
import math
import sys
import textwrap
import tokenize

import docify
import docify_stubs
import libcst as cst
import pytest
from stub_tree import side_effect_free

RUNTIME = '''\
import math

sqrt = math.sqrt
pi = 3.14


def py_func(x):
    """Python doc."""
    return x


class _Documented:
    def __init__(self, doc):
        self.__doc__ = doc


divide = _Documented("Divide element-wise.")
pattern = _Documented("Match \\\\d+ digits.")


class Array(list):
    """Array doc."""
    scale = _Documented("Scale factor.")


class Gen(list):
    pass
'''

STUB = '''\
from typing import overload

@overload
def sqrt(x: int) -> float: ...
@overload
def sqrt(x: float) -> float: ...
def py_func(x): ...
divide: _Documented
pi: float
pattern: _Documented
class _Common:
    def copy(self) -> list: ...
class Array(_Common):
    scale: _Documented
class Gen[T]:
    def copy(self) -> Gen[T]: ...
class Mode:
    None: int
    other: int
'''


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


@pytest.fixture
def documented(site, tmp_path):
    """The fixture stub after the docstring pass, and what it returned."""
    write(site / "fakepkg" / "__init__.py", RUNTIME)
    stub = write(tmp_path / "stubs" / "__init__.pyi", STUB)
    skipped = docify_stubs.document([("fakepkg", stub)], side_effect_free,
                                    {("fakepkg", "_Common"): "Array"})
    return stub.read_text(), skipped


def docstrings(source):
    """Docstring of each def and class, and the string following each annotated
    assignment, as a list per dotted name, one entry per declaration."""
    found = {}

    def visit(body, prefix):
        for i, node in enumerate(body):
            if isinstance(node, (ast.FunctionDef, ast.ClassDef)):
                found.setdefault(prefix + node.name, []).append(ast.get_docstring(node))
                if isinstance(node, ast.ClassDef):
                    visit(node.body, f"{prefix}{node.name}.")
            elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
                after = body[i + 1] if i + 1 < len(body) else None
                doc = None
                if (isinstance(after, ast.Expr) and isinstance(after.value, ast.Constant)
                        and isinstance(after.value.value, str)):
                    doc = after.value.value
                found.setdefault(prefix + node.target.id, []).append(doc)

    visit(ast.parse(source).body, "")
    return found


def test_documents_every_overload_of_a_c_function(documented):
    text, _ = documented
    assert docstrings(text)["sqrt"] == [math.sqrt.__doc__] * 2


def test_leaves_a_python_function_to_the_source_fallback(documented):
    assert docstrings(documented[0])["py_func"] == [None]


def test_documents_an_annotated_assignment_to_a_documented_instance(documented):
    assert docstrings(documented[0])["divide"] == ["Divide element-wise."]


def test_skips_a_constant_whose_doc_is_its_class(documented):
    assert docstrings(documented[0])["pi"] == [None]


def test_documents_a_class_level_assignment(documented):
    assert docstrings(documented[0])["Array.scale"] == ["Scale factor."]


def test_resolves_a_stub_only_class_through_its_alias(documented):
    assert docstrings(documented[0])["_Common.copy"] == [list.copy.__doc__]


def test_documents_a_member_of_a_pep_695_generic_class(documented):
    assert docstrings(documented[0])["Gen.copy"] == [list.copy.__doc__]


def test_backslashes_survive(documented):
    assert docstrings(documented[0])["pattern"] == ["Match \\d+ digits."]


def test_a_keyword_target_is_commented_and_its_module_documented(documented):
    text, skipped = documented
    assert "#     None: int\n" in text
    assert skipped == {}


def test_second_pass_changes_nothing(documented, tmp_path):
    stub = tmp_path / "stubs" / "__init__.pyi"
    docify_stubs.document([("fakepkg", stub)], side_effect_free,
                          {("fakepkg", "_Common"): "Array"})
    assert stub.read_text() == documented[0]


def test_a_module_failing_safe_is_not_imported(site, tmp_path):
    write(site / "fakepkg" / "__init__.py", "")
    write(site / "fakepkg" / "tests" / "__init__.py", "raise RuntimeError('imported')\n")
    stub = write(tmp_path / "tests.pyi", "def f(): ...\n")
    assert docify_stubs.document([("fakepkg.tests", stub)], side_effect_free, {}) == {}
    assert "fakepkg.tests" not in sys.modules


def error_types(skipped):
    return {module: type(e) for module, e in skipped.items()}


def test_a_module_raising_systemexit_is_returned_as_skipped(site, tmp_path):
    write(site / "fakepkg" / "__init__.py", "")
    write(site / "fakepkg" / "exits.py", "raise SystemExit(2)\n")
    stub = write(tmp_path / "exits.pyi", "def f(): ...\n")
    skipped = docify_stubs.document([("fakepkg.exits", stub)], side_effect_free, {})
    assert error_types(skipped) == {"fakepkg.exits": SystemExit}


@pytest.fixture
def two_modules(site, tmp_path):
    """Two importable modules with a C function each, and their stubs."""
    write(site / "fakepkg" / "__init__.py", "")
    for name in ("first", "second"):
        write(site / "fakepkg" / f"{name}.py", "from math import sqrt\n")
    return {name: write(tmp_path / f"{name}.pyi", "def sqrt(x: float) -> float: ...\n")
            for name in ("first", "second")}


def test_a_stub_libcst_rejects_is_returned_and_keeps_its_comments(two_modules):
    first = two_modules["first"]
    first.write_text(first.read_text() + "def f(a=1, b): ...\nclass Mode:\n    None: int\n")
    skipped = docify_stubs.document(
        [(f"fakepkg.{name}", path) for name, path in two_modules.items()], side_effect_free, {})
    assert error_types(skipped) == {"fakepkg.first": cst.ParserSyntaxError}
    assert first.read_text() == ("def sqrt(x: float) -> float: ...\n"
                                 "def f(a=1, b): ...\nclass Mode:\n#     None: int\n")
    assert docstrings(two_modules["second"].read_text())["sqrt"] == [math.sqrt.__doc__]


def test_a_stub_that_does_not_tokenize_is_returned_and_left_alone(two_modules):
    first = two_modules["first"]
    first.write_text('x = """unterminated\n')
    skipped = docify_stubs.document(
        [(f"fakepkg.{name}", path) for name, path in two_modules.items()], side_effect_free, {})
    assert error_types(skipped) == {"fakepkg.first": tokenize.TokenError}
    assert first.read_text() == 'x = """unterminated\n'
    assert docstrings(two_modules["second"].read_text())["sqrt"] == [math.sqrt.__doc__]


def test_patch_docify_is_idempotent():
    docify_stubs.patch_docify()
    docify_stubs.patch_docify()
    assert docify.get_qualname is docify_stubs.get_qualname


def test_multiline_doc_is_indented_to_its_block(site, tmp_path):
    write(site / "fakepkg" / "__init__.py", textwrap.dedent('''\
        class _Documented:
            def __init__(self, doc):
                self.__doc__ = doc
        class Array:
            scale = _Documented("First.\\n\\nSecond.")
        '''))
    stub = write(tmp_path / "__init__.pyi", "class Array:\n    scale: _Documented\n")
    docify_stubs.document([("fakepkg", stub)], side_effect_free, {})
    assert stub.read_text() == ('class Array:\n    scale: _Documented\n'
                                '    """\n    First.\n\n    Second.\n    """\n')
