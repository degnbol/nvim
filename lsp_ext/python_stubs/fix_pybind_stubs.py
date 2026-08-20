#!/usr/bin/env python3
"""Fix a package's bundled pybind11-stubgen stubs (RDKit).

RDKit ships official pybind11-stubgen stubs in a ``rdkit-stubs/`` sibling of
the installed package. They have real, typed function signatures (which
pybind11-stubgen alone cannot produce) but five defects, because
pybind11-stubgen cannot parse Boost.Python's own signatures and does not record
what the package's Python layer does to itself at import time.

Two modes: ``--out`` copies ``<package>-stubs`` to a destination as a *complete*
stub package (so pyright's ``stubPath`` copy wins over any ``<package>-stubs``
installed in the active environment), and ``--in-place`` rewrites the installed
``<package>-stubs`` itself. Both fix these defects:

- Strip the ``C++ signature :`` block from docstrings (keeps the readable
  Boost signature line and prose above it), except where the stub's own return
  annotation is missing or ``Any`` and the C++ line names a concrete type — there
  it is the only type information available, so it stays verbatim.
- Drop ``__init__``/``__new__`` docstrings that are only that C++ noise, so
  pyright falls back to the class docstring on constructor hover.
- Comment out lines whose assignment target is a Python keyword (e.g. an enum
  member named ``None``), which are otherwise a SyntaxError.
- Type each property from a no-arg instance of its class and record the
  default value: ``allowCXSMILES: bool`` with ``(default: True)``.
- Declare the module-level names that only exist at runtime, created by a star
  re-export or by a loop of module-level lambdas. Re-exported where the module
  owning them has a stub that declares them, synthesised from the runtime
  object otherwise.
- Retype the methods the Python layer patches over a C++ class's own, where the
  stub keeps the stale Boost signature (``Mol.GetAtoms``, and ``__iter__`` on
  the ``rdBase`` vector classes).

Every rewritten file gets a ``# fix_pybind_stubs: <package> <version> <sha>``
head comment, where ``<sha>`` is the first 8 hex digits of this script's own
sha256. Comments never reach pyright's hover, so it costs nothing to display,
and it lets a consumer detect not just that a tree was patched but whether it
was patched by *this* version of the script (a reinstall replaces the files and
drops the marker with them, so detection still self-heals).

Must run where ``<package>`` is importable: every stub outside the trees listed
in ``_SIDE_EFFECT_DIRS`` has its module imported to be introspected.
"""
import argparse
import ast
import contextlib
import fnmatch
import hashlib
import importlib
import importlib.metadata
import inspect
import io
import keyword
import os
import pathlib
import pkgutil
import re
import shutil
import sys
import textwrap
import tokenize
import types

MARKER_PREFIX = "# fix_pybind_stubs:"

_CPP_MARKER = re.compile(r"^(\s*)C\+\+ signature :\s*$")
_CPP_RETURN = re.compile(r"^\s*(.*?)\s*[\w:~]+\(.*\)\s*$")
# Boost's stand-ins for "some Python object", which say no more than the stub's
# own `Any`, plus the constructor returns that describe nothing.
_OPAQUE_CPP = frozenset({"void", "void*", "_object*", "PyObject*",
                         "boost::python::api::object", "boost::python::object"})
_WEAK_PY = frozenset({"Any", "typing.Any", "object", "Incomplete"})
_CONSTRUCTOR_DOC = re.compile(
    r'(def (?:__init__|__new__)\([^)]*\)(?:\s*->\s*[^:\n]+)?:\s*)'
    r'("""(?:[^"]|"(?!""))*""")'
)
_KEYWORD_TARGET = re.compile(r"^\s*([A-Za-z_]\w*)\s*[:=]")
_TOPLEVEL_CLASS = re.compile(r"(?m)^(?=class \w)")
_CLASS_NAME = re.compile(r"^class (\w+)")
_PROPERTY = re.compile(
    r'    @property\n'
    r'    def (?P<name>\w+)\(\*args, \*\*kwargs\):\n'
    r'(?:        """\n(?P<doc>(?:.*\n)*?)        """\n|        \.\.\.\n)'
    r'(?P<setter>    @(?P=name)\.setter\n'
    r'    def (?P=name)\(\*args, \*\*kwargs\):\n        \.\.\.\n)?'
)
_METHOD_NAME = re.compile(r"(?m)^    def (\w+)\(")
_FUTURE = "from __future__ import annotations\n"
_ADDED_HEADER = "# present at runtime, absent from the generated stub:\n"

# Trees that compute and print at import time, or need a GUI toolkit. Matched
# per path component, not as a substring of the dotted name, which would also
# drop rdkit.DataStructs and rdkit.DataManip.
_SIDE_EFFECT_DIRS = ("Contrib", "sping", "tests", "TestRunner", "conftest", "demo*")
# Dunders where the stale C++ signature is harmless, so rewriting them from the
# runtime would be churn. __iter__/__next__/__getitem__ are deliberately absent:
# typing those is what makes a patched container iterable and indexable.
_STALE_IS_HARMLESS = frozenset({"__init__", "__new__", "__reduce__",
                                "__getstate__", "__setstate__", "__ge__"})

# Classes whose shape no introspection recovers: rdkit's atom/bond iterators are
# themselves stubgen output of a pure-Python class, so every method is untyped
# and `mol.GetAtoms()[0]` is Unknown. Substituted whole, and asserted to match,
# so an rdkit rename fails loudly instead of dropping the fix.
_CLASS_OVERRIDES = {
    ("rdkit.Chem", "_GetRDKitObjIterator"): """\
class _GetRDKitObjIterator(typing.Generic[_RDKitObj]):
    def __init__(self, mol: Mol) -> None: ...
    def __iter__(self) -> _GetRDKitObjIterator[_RDKitObj]: ...
    def __next__(self) -> _RDKitObj: ...
    def __getitem__(self, i: int) -> _RDKitObj: ...
    def __len__(self) -> int: ...
""",
    ("rdkit.Chem", "_GetAtomsIterator"):
        "class _GetAtomsIterator(_GetRDKitObjIterator[Atom]): ...\n",
    ("rdkit.Chem", "_GetBondsIterator"):
        "class _GetBondsIterator(_GetRDKitObjIterator[Bond]): ...\n",
}
_OVERRIDE_PREAMBLE = {
    "rdkit.Chem": ("import typing\n", "_RDKitObj = typing.TypeVar('_RDKitObj')\n"),
}


def cpp_return_type(line):
    """Return type named by a Boost C++ signature line.

    line: a signature line such as
        ``ExplicitBitVect* GetAvalonFP(RDKit::ROMol,unsigned int)``.
    Returns the text before the function name, empty for a constructor, or None
    if the line is not a signature.
    """
    m = _CPP_RETURN.match(line)
    return m.group(1) if m else None


def keepable_signature_lines(text):
    """Line numbers of C++ signature markers whose block still adds type information.

    text: stub source.
    Returns the 1-based line numbers of ``C++ signature :`` markers inside the
    docstring of a function whose own return annotation is missing or ``Any``,
    where the C++ line names a concrete type. Elsewhere the block only restates
    the stub's own signature.
    """
    lines = text.splitlines()
    # Keyword targets are a syntax error, and line-count preserving to comment
    # out, so the numbers hold against the caller's text.
    tree = ast.parse(comment_keyword_targets(text))
    keep = set()
    for node in ast.walk(tree):
        if not isinstance(node, ast.FunctionDef):
            continue
        if node.returns is not None and ast.unparse(node.returns) not in _WEAK_PY:
            continue
        doc = node.body[0]
        if not (isinstance(doc, ast.Expr) and isinstance(doc.value, ast.Constant)
                and isinstance(doc.value.value, str)):
            continue
        # Each marker is read together with the line below it, so stop one short.
        last = min(doc.end_lineno or doc.lineno, len(lines) - 1)
        for n in range(doc.lineno, last + 1):
            if not _CPP_MARKER.match(lines[n - 1]):
                continue
            ret = cpp_return_type(lines[n])
            if ret and ret not in _OPAQUE_CPP:
                keep.add(n)
    return frozenset(keep)


def strip_cpp_signatures(text, keep=frozenset()):
    """Remove ``C++ signature :`` blocks from docstrings.

    text: stub source.
    keep: 1-based line numbers of marker lines whose block to leave in place.
    Returns the source with every other marker line and its following
    more-indented C++ type lines removed.
    """
    lines = text.splitlines(keepends=True)
    out = []
    i = 0
    while i < len(lines):
        marker = _CPP_MARKER.match(lines[i])
        if not marker or i + 1 in keep:
            out.append(lines[i])
            i += 1
            continue
        indent = len(marker.group(1))
        i += 1
        while i < len(lines):
            nxt = lines[i]
            if nxt.strip() == "":
                i += 1
                break
            if len(nxt) - len(nxt.lstrip()) > indent:
                i += 1
            else:
                break
    return "".join(out)


def drop_constructor_docstrings(text):
    """Drop ``__init__``/``__new__`` docstrings that are only C++ noise.

    text: stub source.
    Returns the source with those docstrings replaced by ``...`` so pyright
    shows the class docstring on constructor hover. Constructors carrying real
    prose (no C++ signature marker) are left untouched.
    """
    def repl(m):
        return m.group(1) + "..." if "C++ signature" in m.group(2) else m.group(0)

    return _CONSTRUCTOR_DOC.sub(repl, text)


def string_literal_lines(text):
    """1-based line numbers any string literal of a source spans.

    text: source that tokenizes, which does not require it to compile.
    Returns the set of those line numbers.
    """
    spans = set()
    for token in tokenize.generate_tokens(io.StringIO(text).readline):
        if token.type == tokenize.STRING:
            spans |= set(range(token.start[0], token.end[0] + 1))
    return spans


def comment_keyword_targets(text):
    """Comment out lines whose assignment target is a Python keyword.

    text: stub source.
    Returns the source with such lines prefixed by ``# ``. Boost.Python enums
    can have a member named ``None``, which is an invalid annotation target;
    the value stays listed in the enum's names/values dicts. Docstring interiors
    are left alone: prose reads as an assignment often enough (``from: <url>``)
    and commenting it out is both wrong and, since the comment then no longer
    matches, not idempotent.
    """
    prose = string_literal_lines(text)
    out = [
        "# " + line
        if (m := _KEYWORD_TARGET.match(line)) and keyword.iskeyword(m.group(1))
        and n not in prose
        else line
        for n, line in enumerate(text.splitlines(keepends=True), start=1)
    ]
    return "".join(out)


def property_types(cls):
    """Map each readable builtin-typed property of a class to its type and default.

    cls: a class to introspect.
    Returns a dict {property_name: (type_name, default_repr)} built from a
    no-arg instance. Empty if the class is not no-arg constructible; properties
    whose value is not a bool/int/float/str are skipped.
    """
    try:
        inst = cls()
    except Exception:
        return {}
    types = {}
    for name, attr in vars(cls).items():
        if name.startswith("_") or not isinstance(attr, property):
            continue
        try:
            value = getattr(inst, name)
        except Exception:
            continue
        if type(value) in (bool, int, float, str):
            types[name] = (type(value).__name__, repr(value))
    return types


def _type_class_properties(block, prop_map):
    """Rewrite the property stubs of one class block using an introspected map.

    block: source of a single top-level class.
    prop_map: {property_name: (type_name, default_repr)} for that class.
    Returns the block with matched properties given a typed getter (and setter
    if the original had one) and the default appended to the docstring.
    """
    def repl(m):
        name = m.group("name")
        if name not in prop_map:
            return m.group(0)
        type_name, default = prop_map[name]
        doc = (m.group("doc") or "").strip()
        doc = f"{doc} (default: {default})" if doc else f"default: {default}"
        out = f'    @property\n    def {name}(self) -> {type_name}:\n        """{doc}"""\n'
        if m.group("setter"):
            out += f"    @{name}.setter\n    def {name}(self, value: {type_name}) -> None: ...\n"
        return out

    return _PROPERTY.sub(repl, block)


def type_properties(text, module):
    """Type the property stubs in a module's source via runtime introspection.

    text: stub source of the module (after text cleanups).
    module: the imported module the stub describes.
    Returns the source with each class's builtin-typed properties rewritten to
    carry their real type and default.
    """
    blocks = _TOPLEVEL_CLASS.split(text)
    out = []
    for block in blocks:
        name = _CLASS_NAME.match(block)
        cls = getattr(module, name.group(1), None) if name else None
        if isinstance(cls, type):
            block = _type_class_properties(block, property_types(cls))
        out.append(block)
    return "".join(out)


class _Default:
    """Stand-in for a parameter default, rendering as ``...`` in a signature."""

    def __repr__(self):
        return "..."


def toplevel_names(text):
    """Names a stub binds at its own top level.

    text: stub source, parseable (keyword targets already commented out).
    Returns the set of names bound by a declaration, an assignment, an
    ``import`` or a ``from ... import``. Star imports are not followed, so a
    name reachable only through one does not count as bound.
    """
    names = set()
    for node in ast.parse(text).body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            names.add(node.name)
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
            names.add(node.target.id)
        elif isinstance(node, ast.Assign):
            names |= {t.id for t in node.targets if isinstance(t, ast.Name)}
        elif isinstance(node, ast.ImportFrom):
            names |= {a.asname or a.name for a in node.names if a.name != "*"}
        elif isinstance(node, ast.Import):
            names |= {(a.asname or a.name).split(".")[0] for a in node.names}
    return names


def ensure_module_line(text, line):
    """Add a module-level statement below the ``__future__`` header if absent.

    text: stub source.
    line: the statement, newline included.
    Returns the source unchanged when the statement is already a line of it, so
    a second run over an already-patched tree leaves the file byte-identical.
    """
    if line in text.splitlines(keepends=True):
        return text
    at = text.find(_FUTURE)
    at = at + len(_FUTURE) if at >= 0 else 0
    return text[:at] + line + text[at:]


def synthesised_def(name, obj):
    """A ``def`` for a runtime callable that no stub declares.

    name: name to declare it under.
    obj: the runtime object, read for its parameter names and docstring.
    Returns the source of a def whose every default is ``...`` — the runtime
    repr is often neither valid Python (a ``<Mol object at 0x…>``) nor stable
    between runs — falling back to ``(*args, **kwargs)`` where the signature is
    unreadable. The return stays unannotated: it is not recoverable without
    calling the object.
    """
    try:
        sig = inspect.signature(obj)
    except (TypeError, ValueError):
        params = "*args, **kwargs"
    else:
        bare = [p.replace(annotation=p.empty,
                          default=p.empty if p.default is p.empty else _Default())
                for p in sig.parameters.values()]
        params = str(sig.replace(parameters=bare,
                                 return_annotation=inspect.Signature.empty))[1:-1]
    # Raw, so a docstring's backslashes reach hover as written: an unrecognised
    # \N{...} or a truncated \u is a SyntaxError, which would take the whole run
    # down. A trailing backslash cannot be expressed raw, so drop that docstring.
    doc = inspect.getdoc(obj) or ""
    body = ('    r"""\n' + textwrap.indent(doc, "    ") + '\n    """\n'
            if doc and '"""' not in doc and not doc.endswith("\\") else "    ...\n")
    return f"def {name}({params}):\n{body}"


def declaration_for(name, obj, module):
    """A stub declaration making one runtime name of a module resolvable.

    name: the attribute name on the module.
    obj: the runtime object it holds.
    module: the imported module.
    Returns the source, or None for an object that is neither a module nor
    callable. A submodule and an object the module that owns it holds under the
    same name are re-exported, which keeps whatever type that module's stub
    gives them; anything else is synthesised here. Identity against the claimed
    owner is what makes the re-export safe: ``__module__`` alone can name a
    module that has no such attribute — an instance inherits it from its class —
    and aliasing that moves the error into the stub instead of fixing it.
    """
    if isinstance(obj, types.ModuleType):
        parent, _, leaf = obj.__name__.rpartition(".")
        if not parent:
            return f"import {leaf} as {name}\n"
        # A package holding its own submodule: relative, so the stub does not
        # import itself by its absolute name.
        where = "." if parent == module.__name__ else parent
        return f"from {where} import {leaf} as {name}\n"
    origin = getattr(obj, "__module__", None)
    if (origin and origin != module.__name__
            and getattr(sys.modules.get(origin), name, None) is obj):
        return f"from {origin} import {name} as {name}\n"
    return synthesised_def(name, obj) if callable(obj) else None


def declare_missing_names(text, module):
    """Declare the module-level names a stub leaves out but the runtime has.

    text: stub source of the module, cleaned.
    module: the imported module the stub describes.
    Returns the source with a declaration appended for each public runtime name
    not already bound at the top level. Skipping the bound ones is what keeps
    this safe: appending makes the last binding win, and a name's runtime owner
    is sometimes a worse home than the declaration the file already carries
    (``rdkit.Chem.inchi`` owns names ``rdkit.Chem``'s stub types better).
    """
    declared = toplevel_names(text)
    added = []
    for name in sorted(dir(module)):
        if name.startswith("_") or name in declared or keyword.iskeyword(name):
            continue
        obj = getattr(module, name, None)
        declaration = declaration_for(name, obj, module) if obj is not None else None
        if declaration:
            added.append(clean(declaration))
    if not added:
        return text
    return text + "\n" + _ADDED_HEADER + "".join(added)


def in_extension_module(cls):
    """Whether a class is defined in a compiled module rather than a ``.py``.

    cls: the class to test.
    Returns True for a Boost.Python class itself, False for a pure-Python
    subclass of one — where every method is a Python function overriding
    nothing, so none of them is a patch over a C++ original.
    """
    home = sys.modules.get(getattr(cls, "__module__", "") or "")
    path = getattr(home, "__file__", None)
    return bool(path) and not path.endswith(".py")


def annotation_for(cls, module):
    """How a runtime class is written as an annotation in a module's stub.

    cls: the class to name.
    module: the imported module whose stub the annotation goes in.
    Returns (annotation, import_line): a bare name for a class of this module or
    of builtins, else a dotted name and the ``import`` statement it needs.
    """
    home = cls.__module__
    if home in (module.__name__, "builtins"):
        return cls.__qualname__, None
    return f"{home}.{cls.__qualname__}", f"import {home}\n"


def boost_override_types(cls, names, module):
    """Annotate each stub-declared method Python patched over the C++ original.

    cls: the class the stub block declares, defined in a C extension module.
    names: the method names that block declares.
    module: the imported module whose stub the annotations go in.
    Returns {method_name: (annotation, import_line)} for those methods whose
    runtime attribute is a plain Python function of the instance alone, so the
    type it returns is readable off a call on a no-arg instance. Ones taking
    arguments would need a fixture value per method, and keep the stale
    annotation the generator gave them.
    """
    patched = [n for n in names if n not in _STALE_IS_HARMLESS
               and isinstance(getattr(cls, n, None), types.FunctionType)
               and len(inspect.signature(getattr(cls, n)).parameters) == 1]
    if not patched:
        return {}
    try:
        instance = cls()
    except Exception:
        return {}
    types_by_name = {}
    for name in patched:
        try:
            returned = getattr(instance, name)()
        except Exception:
            continue
        types_by_name[name] = annotation_for(type(returned), module)
    return types_by_name


def _retype_methods(block, overrides):
    """Rewrite the headers of one class block's patched methods.

    block: source of a single top-level class.
    overrides: {method_name: annotation} for the methods to rewrite.
    Returns the block with each of those declared ``(self) -> annotation``,
    dropping a ``@staticmethod`` the generator emitted for a patched method
    whose first parameter is not called ``self``. Docstrings are untouched.
    """
    for name, annotation in overrides.items():
        header = re.compile(
            rf"(?m)^(?P<indent>[ ]+)(?:@staticmethod\n(?P=indent))?"
            rf"def {re.escape(name)}\([^)]*\)(?:\s*->[^:\n]+)?:")
        block = header.sub(
            lambda m: f"{m.group('indent')}def {name}(self) -> {annotation}:", block)
    return block


def substitute_class(text, name, replacement):
    """Replace one top-level class definition with given source.

    text: stub source.
    name: the class to replace.
    replacement: its new source, newline-terminated.
    Returns the rewritten source. Raises KeyError when the class is absent, so
    an upstream rename fails loudly instead of silently dropping the override.
    """
    for node in ast.parse(text).body:
        if isinstance(node, ast.ClassDef) and node.name == name:
            lines = text.splitlines(keepends=True)
            return "".join(lines[:node.lineno - 1] + [replacement] + lines[node.end_lineno:])
    raise KeyError(f"no class {name} to override")


def type_boost_overrides(text, module):
    """Retype the methods a module patches over an extension class's own.

    text: stub source of the module (after text cleanups).
    module: the imported module the stub describes.
    Returns the source with those methods declared from a live call, and with
    the literal ``_CLASS_OVERRIDES`` entries for this module substituted in.
    """
    imports = []
    out = []
    for block in _TOPLEVEL_CLASS.split(text):
        match = _CLASS_NAME.match(block)
        cls = getattr(module, match.group(1), None) if match else None
        if isinstance(cls, type) and in_extension_module(cls):
            typed = boost_override_types(cls, _METHOD_NAME.findall(block), module)
            imports += [line for _, line in typed.values() if line]
            block = _retype_methods(block, {n: a for n, (a, _) in typed.items()})
        out.append(block)
    text = "".join(out)
    overridden = False
    for (owner, name), replacement in _CLASS_OVERRIDES.items():
        # Keyed on the runtime too, since an override describes a runtime class:
        # a package that no longer has one has nothing to describe, whereas one
        # that has it while its stub does not is the drift worth raising on.
        if owner != module.__name__ or not hasattr(module, name):
            continue
        text = substitute_class(text, name, replacement)
        overridden = True
    if overridden:
        imports += _OVERRIDE_PREAMBLE.get(module.__name__, ())
    # Each goes directly below the __future__ header, so inserting back to front
    # leaves them in the order given: the TypeVar after the import it uses.
    for line in reversed(imports):
        text = ensure_module_line(text, line)
    return text


def clean(text):
    """Apply the import-free text cleanups to one stub's source.

    text: stub source, patched or not; a marker line left by a previous run is
        dropped so the cleanups are idempotent.
    Returns the cleaned source, without a marker.
    """
    if text.startswith(MARKER_PREFIX):
        text = text.split("\n", 1)[-1]
    text = drop_constructor_docstrings(text)
    text = strip_cpp_signatures(text, keepable_signature_lines(text))
    text = comment_keyword_targets(text)
    return text


def package_version(package):
    """Version the installed distribution reports for a package.

    package: import name, which is also the distribution name for the packages
        this fixes.
    Returns the version, or ``unknown`` when no distribution of that name is
    installed — the marker is provenance, so a missing version must not stop the
    stubs being patched.
    """
    try:
        return importlib.metadata.version(package)
    except importlib.metadata.PackageNotFoundError:
        return "unknown"


def script_fingerprint():
    """Short digest of this script, identifying which version of it patched a tree.

    Returns the first 8 hex digits of the file's sha256, so a consumer holding
    the script's path can tell an outdated patched tree from a current one.
    """
    return hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()[:8]


def marker_line(package):
    """Head comment recording what rewrote a stub, and which version it described.

    package: import name whose version is recorded.
    Returns the comment line, newline included.
    """
    return f"{MARKER_PREFIX} {package} {package_version(package)} {script_fingerprint()}\n"


def module_name(pyi_path, stub_root, package):
    """Dotted module name a stub file describes.

    pyi_path: path to a ``.pyi`` file inside stub_root.
    stub_root: the stub tree's root directory.
    package: import prefix the tree describes, so the name is independent of
        stub_root's own directory name.
    Returns the dotted name, dropping a trailing ``__init__``.
    """
    parts = pyi_path.relative_to(stub_root).with_suffix("").parts
    if parts[-1] == "__init__":
        parts = parts[:-1]
    return ".".join((package, *parts))


def side_effect_free(parts):
    """Whether a module's path components steer clear of the import-time hazards.

    parts: the module's path components below the package root, as directory
        names or as the segments of its dotted name.
    Returns False for the demo, test and vendored-toolkit trees, which compute,
    print or need a GUI toolkit at import time.
    """
    return not any(fnmatch.fnmatch(part, pattern)
                   for part in parts for pattern in _SIDE_EFFECT_DIRS)


def introspection_order(stub_root):
    """The tree's stub files, ordered so a package precedes everything under it.

    stub_root: the stub tree's root directory.
    Returns the paths of every ``.pyi``, parents first. Importing a submodule
    binds it as an attribute of its package, so a package enumerated after its
    own descendants would appear to hold names that importing it alone does not
    create.
    """
    return sorted(stub_root.rglob("*.pyi"),
                  key=lambda p: (len(p.relative_to(stub_root).parts),
                                 p.name != "__init__.pyi", p))


def missing_stub_modules(package):
    """Dotted names of a package's modules that its bundled stubs do not cover.

    package: importable package name whose sibling ``<name>-stubs`` holds them.
    Returns the names in walk order, excluding the trees ``side_effect_free``
    rules out. A package that cannot be imported is reported on stderr and its
    subtree skipped, rather than dropped silently as pkgutil would.
    """
    bundled = stubs_path(package)
    root = importlib.import_module(package)
    missing = []
    with open(os.devnull, "w") as quiet, contextlib.redirect_stdout(quiet):
        for info in pkgutil.walk_packages(
                root.__path__, package + ".",
                onerror=lambda name: print(f"not walked: {name}", file=sys.stderr)):
            parts = info.name.split(".")[1:]
            relative = pathlib.Path(*parts)
            if (side_effect_free(parts)
                    and not (bundled / relative.with_suffix(".pyi")).exists()
                    and not (bundled / relative / "__init__.pyi").exists()):
                missing.append(info.name)
    return missing


def patch_dir(stub_root, package):
    """Fix the Boost.Python defects in every ``.pyi`` under a stub tree, in place.

    stub_root: directory holding the stub tree.
    package: import prefix the tree describes (e.g. ``rdkit``).
    Returns the dotted module names whose import failed, leaving them with the
    text cleanups alone.
    """
    root = stub_root.resolve()
    marker = marker_line(package)
    unimportable = []
    for p in introspection_order(root):
        # Explicit encoding: a few stubs are non-ASCII and the locale is the
        # caller's (an editor's environment), not necessarily UTF-8.
        text = clean(p.read_text(encoding="utf-8"))
        name = module_name(p, root, package)
        # Every module is introspected, so every importable one is imported.
        if side_effect_free(p.relative_to(root).with_suffix("").parts):
            try:
                module = importlib.import_module(name)
            except Exception as e:
                unimportable.append(f"{name} ({e})")
            else:
                text = type_properties(text, module)
                text = type_boost_overrides(text, module)
                text = declare_missing_names(text, module)
        p.write_text(marker + text, encoding="utf-8")
    return unimportable


def stubs_path(package):
    """Path of a package's bundled ``<package>-stubs`` sibling directory.

    package: importable package name.
    Returns the path, which need not exist.
    """
    init = importlib.import_module(package).__file__
    if init is None:
        raise SystemExit(f"{package} is a namespace package with no __file__")
    return pathlib.Path(init).parent.parent / f"{package}-stubs"


def vendor(package, dest):
    """Copy a package's bundled ``<package>-stubs`` into dest.

    package: importable package name whose sibling ``<name>-stubs`` holds the
        official stubs.
    dest: destination stub package directory. Merged into rather than replaced,
        so a caller can put stubs for the modules the bundle omits there first
        and have them fixed alongside; clearing it is that caller's business.
    Returns dest.
    """
    source = stubs_path(package)
    if not source.is_dir():
        raise SystemExit(f"no bundled stubs at {source}")
    shutil.copytree(source, dest, dirs_exist_ok=True)
    return dest


def patch_staged(stub_root, package):
    """Fix a stub tree in place, swapping a patched copy in atomically.

    stub_root: directory holding the stub tree; it is replaced by the patched
        copy built beside it, so a reader never sees a half-rewritten tree and
        an interrupted run leaves the original intact. A swap interrupted
        between its two renames is recovered on the next call.
    package: import prefix the tree describes (e.g. ``rdkit``).
    Returns the dotted module names whose import failed, leaving them with the
    text cleanups alone.
    """
    backup = stub_root.with_name(stub_root.name + ".bak")
    staging = stub_root.with_name(stub_root.name + ".new")
    if not stub_root.is_dir() and backup.is_dir():
        backup.rename(stub_root)
    if not stub_root.is_dir():
        raise SystemExit(f"no stubs at {stub_root}")
    for leftover in (backup, staging):
        if leftover.exists():
            shutil.rmtree(leftover)
    shutil.copytree(stub_root, staging)
    unimportable = patch_dir(staging, package)
    stub_root.rename(backup)
    staging.rename(stub_root)
    shutil.rmtree(backup)
    return unimportable


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("package", help="importable package name (e.g. rdkit)")
    where = ap.add_mutually_exclusive_group(required=True)
    where.add_argument("-o", "--out", type=pathlib.Path,
                       help="copy the stubs to this directory and fix the copy")
    where.add_argument("-i", "--in-place", action="store_true",
                       help="fix the stubs installed in the active environment")
    args = ap.parse_args()
    # A couple of rdkit's modules print on import; the report below goes to the
    # real stdout, which the redirect has already been lifted from.
    with open(os.devnull, "w") as quiet, contextlib.redirect_stdout(quiet):
        if args.in_place:
            unimportable = patch_staged(stubs_path(args.package), args.package)
        else:
            unimportable = patch_dir(vendor(args.package, args.out), args.package)
    if unimportable:
        print("introspection skipped (import failed): " + ", ".join(unimportable))


if __name__ == "__main__":
    main()
