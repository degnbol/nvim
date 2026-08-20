#!/usr/bin/env python3
"""Tests for the Boost.Python stub fixes. Run: ./test_fix_pybind_stubs.py

``--slow`` additionally runs the checks whose subject is what basedpyright makes
of the vendored tree beside this file. Those need rdkit and basedpyright:
``uv run --no-project --with rdkit --with basedpyright ./test_fix_pybind_stubs.py --slow``.
"""
import ast
import contextlib
import hashlib
import json
import pathlib
import shutil
import subprocess
import sys
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


def test_marker_names_the_package_version_and_the_script_that_wrote_it():
    # An uninstalled package also covers the version fallback.
    assert v.marker_line("nosuchpkg") == (
        f"# fix_pybind_stubs: nosuchpkg unknown {v.script_fingerprint()}\n")


def test_fingerprint_is_the_first_8_of_the_script_sha256():
    # The Lua side recomputes this with vim.fn.sha256 over the same bytes.
    digest = hashlib.sha256(pathlib.Path(v.__file__).read_bytes()).hexdigest()
    assert v.script_fingerprint() == digest[:8]


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
        assert out.startswith(v.marker_line("fake"))
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


def fake_module(name, **attrs):
    """Register an importable module holding the given attributes, and return it."""
    module = types.ModuleType(name)
    module.__dict__.update(attrs)
    sys.modules[name] = module
    return module


def owned(module, name, obj):
    """Put an object on a module and make it claim that module as its owner."""
    obj.__module__ = module.__name__
    setattr(module, name, obj)
    return obj


def test_synthesises_a_lambda_built_name_with_ellipsis_defaults():
    # The descriptor modules build names in a loop, so no stub declares them
    # anywhere; the runtime default is a Mol whose repr is not even valid Python.
    module = fake_module("fakedesc")
    owned(module, "NOCount", lambda x, countUnique=True, pattern=object(): 0)
    out = v.declare_missing_names("", module)
    assert "def NOCount(x, countUnique=..., pattern=...):" in out


def test_does_not_alias_a_name_its_claimed_owner_lacks():
    # An instance inherits __module__ from its class, so the owner it names can
    # have no such attribute. Aliasing that moves the error into the stub.
    def GenKeys():
        pass

    owned(fake_module("fakeowner"), "under_another_name", GenKeys)
    out = v.declare_missing_names("", fake_module("fakeuser", GenKeys=GenKeys))
    assert "import GenKeys" not in out
    assert "def GenKeys():" in out


def test_reexports_a_name_its_owner_holds_under_the_same_name():
    def AddHs(mol):
        return mol

    owned(fake_module("fakerdmolops"), "AddHs", AddHs)
    module = fake_module("fakeallchem", AddHs=AddHs)
    assert "from fakerdmolops import AddHs as AddHs\n" in v.declare_missing_names("", module)


def test_reexports_a_submodule_instead_of_synthesising_it():
    package = fake_module("fakechem", rdchem=fake_module("fakechem.rdchem"))
    assert "from . import rdchem as rdchem\n" in v.declare_missing_names("", package)


def test_leaves_a_name_the_stub_already_imports_alone():
    # The rdkit.Chem.inchi regression: the stub types MolToInchiKey better than
    # the Python wrapper that owns it does.
    def MolToInchiKey(mol, options=""):
        pass

    owner = fake_module("fakeinchi", MolToInchiKey=MolToInchiKey)
    module = fake_module("fakechem2", MolToInchiKey=owner.MolToInchiKey)
    text = "from fakeinchi import MolToInchiKey\n"
    assert v.declare_missing_names(text, module) == text


def test_retype_drops_the_staticmethod_the_generator_emitted():
    # rdkit patches __iter__ with a function whose parameter is not called self,
    # which pybind11-stubgen reads as a staticmethod.
    block = "class V:\n    @staticmethod\n    def __iter__(vect):\n        ...\n"
    assert v._retype_methods(block, {"__iter__": "rdkit.VectIter"}) == (
        "class V:\n    def __iter__(self) -> rdkit.VectIter:\n        ...\n")


def test_retype_replaces_the_return_and_keeps_the_docstring():
    block = ('class Mol:\n    def GetAtoms(self) -> typing.Iterable[Atom]:\n'
             '        """\n        an iterator\n        """\n')
    assert v._retype_methods(block, {"GetAtoms": "rdkit.Chem._GetAtomsIterator"}) == (
        'class Mol:\n    def GetAtoms(self) -> rdkit.Chem._GetAtomsIterator:\n'
        '        """\n        an iterator\n        """\n')


def test_extension_module_test_tells_a_native_class_from_a_python_subclass():
    # A pure-Python subclass of a Boost class has only Python functions, none of
    # which overrides anything. Where the class is defined is what separates them.
    fake_module("fakeboost", __file__="/x/fakeboost.so")

    class Native:
        pass

    Native.__module__ = "fakeboost"

    class Subclass(Native):
        pass

    assert v.in_extension_module(Native)
    assert not v.in_extension_module(Subclass)


def test_keyword_commenting_leaves_docstring_prose_alone():
    # `from: <url>` inside a docstring reads as an annotated assignment.
    src = 'def f():\n    """\n    Plane of best fit\n    from: https://example.org\n    """\n'
    assert v.comment_keyword_targets(src) == src
    assert v.comment_keyword_targets("    None: int\n") == "#     None: int\n"


def test_side_effect_exclusion_matches_whole_path_components():
    assert v.side_effect_free(("DataStructs", "cDataStructs"))
    assert v.side_effect_free(("DataManip", "Metric"))
    assert not v.side_effect_free(("Chem", "Subshape", "demoCombined"))
    assert not v.side_effect_free(("sping", "PDF"))


def test_module_line_is_inserted_once_below_the_future_header():
    text = "from __future__ import annotations\nclass X: ...\n"
    once = v.ensure_module_line(text, "import typing\n")
    assert once == "from __future__ import annotations\nimport typing\nclass X: ...\n"
    assert v.ensure_module_line(once, "import typing\n") == once


def test_class_override_replaces_only_its_own_class():
    text = "class A:\n    def f(self): ...\nclass B: ...\ndef g(): ...\n"
    assert v.substitute_class(text, "A", "class A: ...\n") == (
        "class A: ...\nclass B: ...\ndef g(): ...\n")


def test_class_override_raises_when_its_target_is_gone():
    raised = False
    try:
        v.substitute_class("class B: ...\n", "A", "class A: ...\n")
    except KeyError:
        raised = True
    assert raised


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


@contextlib.contextmanager
def importable_stub_tree():
    """Yield a one-file stub tree whose module is importable and has an extra name."""
    with tempfile.TemporaryDirectory() as tmp:
        root = pathlib.Path(tmp) / "runtime-stubs"
        (root / "Chem").mkdir(parents=True)
        (root / "Chem" / "rdmolfiles.pyi").write_text(BLOCK)
        owned(fake_module("runtime.Chem.rdmolfiles"), "MolWt", lambda *x, **y: 0)
        yield root


def test_excluded_modules_are_not_introspected():
    # The exclusion is by path component, so the walk has to strip `.pyi` before
    # testing it — `conftest.pyi` is not `conftest`.
    with importable_stub_tree() as root:
        conftest = root / "conftest.pyi"
        conftest.write_text("x: int\n")
        owned(fake_module("runtime.conftest"), "computed_at_import", lambda: 0)
        v.patch_staged(root, "runtime")
        assert "computed_at_import" not in conftest.read_text()


def test_runtime_passes_leave_the_tree_byte_identical_on_a_second_run():
    # The property that the --in-place re-patching of an environment rests on.
    with importable_stub_tree() as root:
        pyi = root / "Chem" / "rdmolfiles.pyi"
        assert v.patch_staged(root, "runtime") == []
        once = pyi.read_text()
        assert "def MolWt(*x, **y):" in once
        v.patch_staged(root, "runtime")
        assert pyi.read_text() == once


PROBE = '''\
from rdkit import Chem
from rdkit.Chem import AllChem, Lipinski, PandasTools

mol = Chem.MolFromSmiles("CCO")
assert mol is not None
atoms = mol.GetAtoms()
reveal_type(atoms)
reveal_type(atoms[0])
reveal_type(len(atoms))
reveal_type([a.GetSymbol() for a in atoms])
reveal_type(Chem.InchiToInchiKey)
reveal_type(Chem.SmilesParserParams().allowCXSMILES)
AllChem.AddHs(mol)
PandasTools.LoadSDF("x.sdf")
Lipinski.NOCount()
'''
REVEALED = [
    '"atoms" is "_GetAtomsIterator"',       # the class override
    '"atoms[0]" is "Atom"',
    '"len(atoms)" is "int"',
    '"[a.GetSymbol() for a in atoms]" is "list[str]"',
    # Owned by the rdkit.Chem.inchi wrapper, but typed better where it is
    # declared: the name is bound at the stub's top level, so nothing re-exports
    # it over the top.
    '"Chem.InchiToInchiKey" is "(inchi: str) -> str"',
    '"Chem.SmilesParserParams().allowCXSMILES" is "bool"',
]


def slow_test_pyright_reads_the_vendored_tree_as_intended():
    """Type the probe against the tree beside this file, and check what comes back.

    The only check here that catches a change in pyright's own behaviour, and
    the harness the passes were measured with.
    """
    stubs = pathlib.Path(__file__).parent
    with tempfile.TemporaryDirectory() as tmp:
        work = pathlib.Path(tmp)
        (work / "pyrightconfig.json").write_text(json.dumps(
            {"typeCheckingMode": "standard", "stubPath": str(stubs),
             "reportMissingModuleSource": "none"}))
        (work / "probe.py").write_text(PROBE)
        out = subprocess.run(["basedpyright", "--outputjson", "probe.py"],
                             cwd=work, capture_output=True, text=True).stdout
    diagnostics = json.loads(out)["generalDiagnostics"]
    notes = [d["message"] for d in diagnostics if d["severity"] == "information"]
    for revealed in REVEALED:
        assert any(revealed in note for note in notes), revealed
    errors = [d["message"] for d in diagnostics if d["severity"] == "error"]
    # AllChem.AddHs and PandasTools.LoadSDF resolve, and the sole error is what a
    # synthesised declaration is meant to keep: its arity is still checked.
    assert errors == ['Argument missing for parameter "x"'], errors


if __name__ == "__main__":
    prefixes = ("test_", "slow_test_") if "--slow" in sys.argv else ("test_",)
    tests = [f for name, f in sorted(globals().items()) if name.startswith(prefixes)]
    for t in tests:
        t()
        print(f"ok {t.__name__}")
    print(f"{len(tests)} passed")
