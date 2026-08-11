#!/usr/bin/env python3
"""Fix a package's bundled pybind11-stubgen stubs (RDKit).

RDKit ships official pybind11-stubgen stubs in a ``rdkit-stubs/`` sibling of
the installed package. They have real, typed function signatures (which
pybind11-stubgen alone cannot produce) but three defects, because
pybind11-stubgen cannot parse Boost.Python's own signatures: the raw C++
signature is left in docstrings, enum members named after Python keywords are
emitted as invalid syntax, and every property is typed ``(*args, **kwargs) ->
Unknown``.

Two modes: ``--out`` copies ``<package>-stubs`` to a destination as a *complete*
stub package (so pyright's ``stubPath`` copy wins over any ``<package>-stubs``
installed in the active environment), and ``--in-place`` rewrites the installed
``<package>-stubs`` itself. Both fix these defects:

- Strip the ``C++ signature :`` block from docstrings (keeps the readable
  Boost signature line and prose above it).
- Drop ``__init__``/``__new__`` docstrings that are only that C++ noise, so
  pyright falls back to the class docstring on constructor hover.
- Comment out lines whose assignment target is a Python keyword (e.g. an enum
  member named ``None``), which are otherwise a SyntaxError.
- Type each property from a no-arg instance of its class and record the
  default value: ``allowCXSMILES: bool`` with ``(default: True)``.

Must run where ``<package>`` is importable (for the property introspection).
"""
import argparse
import importlib
import keyword
import pathlib
import re
import shutil

_CPP_MARKER = re.compile(r"^(\s*)C\+\+ signature :\s*$")
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


def strip_cpp_signature(text):
    """Remove ``C++ signature :`` blocks from docstrings.

    text: stub source.
    Returns the source with each marker line and its following more-indented
    C++ type lines removed.
    """
    lines = text.splitlines(keepends=True)
    out = []
    i = 0
    while i < len(lines):
        marker = _CPP_MARKER.match(lines[i])
        if not marker:
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


def comment_keyword_targets(text):
    """Comment out lines whose assignment target is a Python keyword.

    text: stub source.
    Returns the source with such lines prefixed by ``# ``. Boost.Python enums
    can have a member named ``None``, which is an invalid annotation target;
    the value stays listed in the enum's names/values dicts.
    """
    out = [
        "# " + line
        if (m := _KEYWORD_TARGET.match(line)) and keyword.iskeyword(m.group(1))
        else line
        for line in text.splitlines(keepends=True)
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


def clean(text):
    """Apply the import-free text cleanups to one stub's source.

    text: stub source.
    Returns the cleaned source.
    """
    text = drop_constructor_docstrings(text)
    text = strip_cpp_signature(text)
    text = comment_keyword_targets(text)
    return text


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


def patch_dir(stub_root, package):
    """Fix the Boost.Python defects in every ``.pyi`` under a stub tree, in place.

    stub_root: directory holding the stub tree.
    package: import prefix the tree describes (e.g. ``rdkit``).
    Returns the dotted module names whose import failed, leaving their
    properties untyped.
    """
    root = stub_root.resolve()
    unimportable = []
    for p in root.rglob("*.pyi"):
        # Explicit encoding: a few stubs are non-ASCII and the locale is the
        # caller's (an editor's environment), not necessarily UTF-8.
        text = clean(p.read_text(encoding="utf-8"))
        # Only modules with untyped Boost properties need importing; skipping
        # the rest avoids importing pure-Python demos with side effects.
        if _PROPERTY.search(text):
            name = module_name(p, root, package)
            try:
                module = importlib.import_module(name)
            except Exception as e:
                unimportable.append(f"{name} ({e})")
            else:
                text = type_properties(text, module)
        p.write_text(text, encoding="utf-8")
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
    dest: destination stub package directory (replaced if it exists).
    Returns dest.
    """
    source = stubs_path(package)
    if not source.is_dir():
        raise SystemExit(f"no bundled stubs at {source}")
    if dest.exists():
        shutil.rmtree(dest)
    shutil.copytree(source, dest)
    return dest


def patch_staged(stub_root, package):
    """Fix a stub tree in place, swapping a patched copy in atomically.

    stub_root: directory holding the stub tree; it is replaced by the patched
        copy built beside it, so a reader never sees a half-rewritten tree and
        an interrupted run leaves the original intact. A swap interrupted
        between its two renames is recovered on the next call.
    package: import prefix the tree describes (e.g. ``rdkit``).
    Returns the dotted module names whose import failed, leaving their
    properties untyped.
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
    if args.in_place:
        unimportable = patch_staged(stubs_path(args.package), args.package)
    else:
        unimportable = patch_dir(vendor(args.package, args.out), args.package)
    if unimportable:
        print("properties left untyped (import failed): " + ", ".join(unimportable))


if __name__ == "__main__":
    main()
