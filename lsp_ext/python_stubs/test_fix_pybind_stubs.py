#!/usr/bin/env python3
"""Tests for the Boost.Python stub fixes. Run: ./test_fix_pybind_stubs.py"""
import ast
import contextlib
import pathlib
import shutil
import tempfile
import types

import fix_pybind_stubs as v

# A class block in the exact shape RDKit's bundled stubs use.
BLOCK = '''\
class SmilesParserParams(Boost.Python.instance):
    """
    Parameters controlling SMILES Parsing
    """
    def __init__(self) -> None:
        """
            C++ signature :
                void __init__(_object*)
        """
    @property
    def allowCXSMILES(*args, **kwargs):
        """
        controls whether or not the CXSMILES extensions are parsed
        """
    @allowCXSMILES.setter
    def allowCXSMILES(*args, **kwargs):
        ...
    @staticmethod
    def flush(*args, **kwargs):
        """
        flush( (Writer)self) -> None :
            Flushes the output file.

            C++ signature :
                void flush(RDKit::Writer {lvalue})
        """
    def GetFingerprint(self) -> typing.Any:
        """
            Returns the fingerprint.

            C++ signature :
                ExplicitBitVect* GetFingerprint(RDKit::ROMol)
        """
'''


def test_drops_cpp_only_constructor_docstring():
    out = v.drop_constructor_docstrings(BLOCK)
    assert "void __init__" not in out
    assert "def __init__(self) -> None:" in out
    ast.parse(out)


def test_strips_cpp_block_keeps_prose():
    out = v.strip_cpp_signatures(BLOCK)
    assert "C++ signature" not in out
    assert "Flushes the output file." in out  # prose above the block survives


def test_keeps_the_block_that_is_the_only_type_source():
    # GetFingerprint is annotated `-> typing.Any`, so its C++ return is the only
    # place the real type appears; flush's C++ return is `void`, which is not.
    out = v.strip_cpp_signatures(BLOCK, v.keepable_signature_lines(BLOCK))
    assert "ExplicitBitVect* GetFingerprint" in out
    assert "void flush" not in out


def test_keeps_no_block_for_a_fully_typed_def():
    src = ('def f(mol: Mol) -> int:\n    """\n        C++ signature :\n'
           "            int f(RDKit::ROMol)\n" '    """\n')
    assert v.keepable_signature_lines(src) == set()


def test_keepable_lines_survive_a_keyword_enum_member():
    src = ('class E(Boost.Python.enum):\n'
           '    None: typing.ClassVar[E]  # value = E.None\n'
           '    OTHER: typing.ClassVar[E]\n'
           'def f(mol) -> typing.Any:\n    """\n        C++ signature :\n'
           "            ExplicitBitVect* f(RDKit::ROMol)\n" '    """\n')
    assert v.keepable_signature_lines(src) == {6}


def test_marker_is_a_comment_naming_the_package_and_version():
    # An uninstalled package also covers the version fallback.
    assert v.marker_line("nosuchpkg") == "# fix_pybind_stubs: nosuchpkg unknown\n"


def test_comments_keyword_enum_member():
    src = "    None: typing.ClassVar[E]  # value = E.None\n"
    assert v.comment_keyword_targets(src).startswith("# ")


def test_types_property_with_default():
    class Fake:
        @property
        def allowCXSMILES(self):
            return True

    module = types.SimpleNamespace(SmilesParserParams=Fake)
    out = v.type_properties(v.clean(BLOCK), module)
    assert "def allowCXSMILES(self) -> bool:" in out
    assert "(default: True)" in out
    assert "def allowCXSMILES(self, value: bool) -> None: ..." in out  # setter kept


def test_result_is_valid_python():
    ast.parse(v.type_properties(v.clean(BLOCK), types.SimpleNamespace()))


def test_type_properties_is_idempotent():
    class Fake:
        @property
        def allowCXSMILES(self):
            return True

    module = types.SimpleNamespace(SmilesParserParams=Fake)
    once = v.type_properties(v.clean(BLOCK), module)
    assert v.type_properties(once, module) == once


@contextlib.contextmanager
def stub_tree():
    """Yield a one-file stub tree in a temp dir, shaped like the bundled stubs."""
    with tempfile.TemporaryDirectory() as tmp:
        root = pathlib.Path(tmp) / "fake-stubs"
        (root / "Chem").mkdir(parents=True)
        (root / "Chem" / "rdmolfiles.pyi").write_text(BLOCK)
        yield root


def test_patch_staged_rewrites_the_tree_in_place():
    with stub_tree() as root:
        # The fixture package is not importable, so its properties stay untyped
        # and are reported; the text cleanups need no import.
        assert v.patch_staged(root, "fake") == [
            "fake.Chem.rdmolfiles (No module named 'fake')"
        ]
        out = (root / "Chem" / "rdmolfiles.pyi").read_text()
        assert out.startswith("# fix_pybind_stubs: fake unknown\n")
        assert "void flush" not in out            # redundant block dropped
        assert "ExplicitBitVect* GetFingerprint" in out  # informative one kept
        assert not list(root.parent.glob("fake-stubs.*"))  # no staging/backup left


def test_patch_staged_is_idempotent():
    with stub_tree() as root:
        pyi = root / "Chem" / "rdmolfiles.pyi"
        v.patch_staged(root, "fake")
        once = pyi.read_text()
        v.patch_staged(root, "fake")
        assert pyi.read_text() == once


def _die(*_):
    raise RuntimeError("interrupted")


def test_patch_staged_leaves_the_original_on_failure():
    with stub_tree() as root:
        patch_dir, v.patch_dir = v.patch_dir, _die
        raised = False
        try:
            v.patch_staged(root, "fake")
        except RuntimeError:
            raised = True
        finally:
            v.patch_dir = patch_dir
        assert raised
        assert (root / "Chem" / "rdmolfiles.pyi").read_text() == BLOCK


def test_patch_staged_recovers_an_interrupted_swap():
    with stub_tree() as root:
        # State left by a kill between the two renames: the patched tree is
        # staged, the original is only reachable under .bak.
        shutil.copytree(root, root.with_name(root.name + ".new"))
        root.rename(root.with_name(root.name + ".bak"))
        v.patch_staged(root, "fake")
        assert (root / "Chem" / "rdmolfiles.pyi").read_text().startswith(v.MARKER_PREFIX)
        assert not list(root.parent.glob("fake-stubs.*"))


if __name__ == "__main__":
    tests = [f for name, f in sorted(globals().items()) if name.startswith("test_")]
    for t in tests:
        t()
        print(f"ok {t.__name__}")
    print(f"{len(tests)} passed")
