"""Import-free helpers for walking and cleaning a stub tree."""
from __future__ import annotations

import fnmatch
import io
import keyword
import pathlib
import re
import tokenize
from collections.abc import Sequence

_KEYWORD_TARGET = re.compile(r"^\s*([A-Za-z_]\w*)\s*[:=]")

# Trees that compute and print at import time, or need a GUI toolkit. Matched
# per path component, not as a substring of the dotted name, which would also
# drop rdkit.DataStructs and rdkit.DataManip.
_SIDE_EFFECT_DIRS = ("Contrib", "sping", "tests", "TestRunner", "conftest", "demo*",
                     "__main__")


def side_effect_free(parts: Sequence[str]) -> bool:
    """Whether a module lies outside the trees with import-time side effects.

    Args:
        parts: the module's path components below the package root, as directory
            names or as the segments of its dotted name.

    Returns:
        False for the demo, test, ``__main__`` and vendored-toolkit trees, which
        compute, print or need a GUI toolkit at import time.
    """
    return not any(fnmatch.fnmatch(part, pattern)
                   for part in parts for pattern in _SIDE_EFFECT_DIRS)


def module_name(pyi_path: pathlib.Path, stub_root: pathlib.Path, package: str) -> str:
    """Dotted module name a stub file describes.

    Args:
        pyi_path: path to a ``.pyi`` file inside stub_root.
        stub_root: the stub tree's root directory.
        package: import prefix the tree describes, so the name is independent of
            stub_root's own directory name.

    Returns:
        The dotted name, dropping a trailing ``__init__``.
    """
    parts = pyi_path.relative_to(stub_root).with_suffix("").parts
    if parts[-1] == "__init__":
        parts = parts[:-1]
    return ".".join((package, *parts))


def introspection_order(stub_root: pathlib.Path) -> list[pathlib.Path]:
    """The tree's stub files, ordered so a package precedes everything under it.

    Importing a submodule binds it as an attribute of its package, so a package
    enumerated after its own descendants would appear to hold names that
    importing it alone does not create.

    Args:
        stub_root: the stub tree's root directory.

    Returns:
        The paths of every ``.pyi``, parents first.
    """
    return sorted(stub_root.rglob("*.pyi"),
                  key=lambda p: (len(p.relative_to(stub_root).parts),
                                 p.name != "__init__.pyi", p))


def string_literal_lines(text: str) -> set[int]:
    """1-based line numbers any string literal of a source spans.

    Args:
        text: source that tokenizes, which does not require it to compile.

    Returns:
        The set of those line numbers.
    """
    spans = set()
    for token in tokenize.generate_tokens(io.StringIO(text).readline):
        if token.type == tokenize.STRING:
            spans |= set(range(token.start[0], token.end[0] + 1))
    return spans


def comment_keyword_targets(text: str) -> str:
    """Comment out lines whose assignment target is a Python keyword.

    An enum can have a member named ``None``, which is an invalid annotation
    target. Docstring interiors are left alone: prose reads as an assignment
    often enough (``from: <url>``) and commenting it out is both wrong and, since
    the comment then no longer matches, not idempotent.

    Args:
        text: stub source.

    Returns:
        The source with such lines prefixed by ``# ``.
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
