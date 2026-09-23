"""Run the patch sequence on one stub tree or single-module stub.

    uv run --no-project --python <env python> \\
        --with-requirements <dir>/requirements.txt <dir>/patch_stubs.py <stub-dir-or-.pyi>

The sequence is the docstring pass (``docify_stubs``), then the package's entry
in ``REPAIRS`` if it has one, run on a staged copy and written back file by
file, the marker file last. The head line of the marker file (``marker_file``)
then reads
``# patch_stubs: <package> <version> <digest>``, the digest covering every
patch script.

Exit status: 0 when files were patched or modules skipped, each skipped module
on its own stdout line with its error. 3 (``NOTHING_WRITTEN``) when the stub is
current, not applicable, or got only the marker. Anything else is a failure,
with the reason on stderr.
"""
from __future__ import annotations

import contextlib
import faulthandler
import fcntl
import hashlib
import importlib.metadata
import importlib.util
import os
import pathlib
import re
import shutil
import subprocess
import sys
import tempfile
import tokenize
from collections.abc import Callable, Iterator, Sequence

import patch_pybind_stubs
from stdio import discarded_stdout
from stub_tree import introspection_order, module_name, side_effect_free

REPAIRS: dict[str, Callable[[pathlib.Path, str], dict[str, BaseException]]] = {
    "rdkit": patch_pybind_stubs.patch_dir,
}
"""Package -> repair run on its documented stub tree, given the tree's root and
the package, returning the modules it skipped, each with its error."""

CLASS_ALIASES: dict[tuple[str, str], str] = {
    ("numpy", "_ArrayOrScalarCommon"): "ndarray",
}
"""(module, stub-only class) -> runtime class standing for it. Keyed on the
module, since a bare class name would also alias a same-named class elsewhere."""

NOTHING_WRITTEN = 3
WATCHDOG_S = 600

_HERE = pathlib.Path(__file__).resolve().parent
_MARKER_PREFIX = "# patch_stubs:"
_LEGACY_MARKER_PREFIX = "# fix_pybind_stubs:"
"""Written on every file by this sequence's predecessor."""
_LIBCST_ERROR_PREFIX = re.compile(r"^\w+ error: (?:error at \d+:\d+: )?")
"""Leads libcst's error messages (``parser error: error at 2:13: expected ...``),
restating what the line number says."""
_TEMP_PREFIX = ".patch_stubs-"
_TEMP_SUFFIX = ".tmp"


def sources(directory: pathlib.Path = _HERE) -> list[pathlib.Path]:
    """The patch scripts and their requirements file.

    Args:
        directory: where the patch scripts live.

    Returns:
        ``requirements.txt`` and the ``*.py`` files whose name starts with
        neither ``test_`` nor ``.``, sorted by name.
    """
    scripts = [p for p in directory.glob("*.py") if not p.name.startswith(("test_", "."))]
    return sorted([directory / "requirements.txt", *scripts], key=lambda p: p.name)


def fingerprint(paths: Sequence[pathlib.Path]) -> str:
    """Short digest of files' contents.

    Args:
        paths: the files, in the order their bytes are concatenated.

    Returns:
        The first 8 hex digits of the sha256 of that concatenation.
    """
    digest = hashlib.sha256()
    for path in paths:
        digest.update(path.read_bytes())
    return digest.hexdigest()[:8]


def distribution_version(package: str) -> str:
    """Version of the distribution providing an import name.

    Args:
        package: top-level import name, which need not be the distribution's
            name (``cv2`` is opencv-python's).

    Returns:
        The version, or ``unknown`` when no installed distribution provides it.
    """
    distributions = importlib.metadata.packages_distributions().get(package)
    return importlib.metadata.version(distributions[0]) if distributions else "unknown"


def marker_file(target: pathlib.Path) -> pathlib.Path:
    """File whose head line carries the marker.

    Args:
        target: a stub tree's directory, or a single-module ``.pyi``.

    Returns:
        The tree's top-level ``__init__.pyi``, or the file itself for a single module.
    """
    return target / "__init__.pyi" if target.is_dir() else target


def _marker(package: str) -> str:
    """Marker line for a package as currently installed and patched, newline included."""
    return (f"{_MARKER_PREFIX} {package} {distribution_version(package)} "
            f"{fingerprint(sources())}\n")


def _without_marker(text: str) -> str:
    """Stub source with a head line left by this sequence or its predecessor removed."""
    if text.startswith((_MARKER_PREFIX, _LEGACY_MARKER_PREFIX)):
        return text.split("\n", 1)[-1]
    return text


def is_current(target: pathlib.Path, package: str) -> bool:
    """Whether a stub already carries the marker the sequence would write now.

    Args:
        target: a stub tree's directory, or a single-module ``.pyi``.
        package: import name the stub describes.

    Returns:
        True when the marker names this package, its installed version and
        ``fingerprint(sources())``.
    """
    with marker_file(target).open(encoding="utf-8") as f:
        return f.readline() == _marker(package)


def applicable(target: pathlib.Path, package: str) -> bool:
    """Whether the sequence may run on a stub at all.

    Args:
        target: a stub tree's directory, or a single-module ``.pyi``.
        package: import name the stub describes.

    Returns:
        False when ``importlib.util.find_spec`` finds no runtime package (a
        stub-only distribution) or only a namespace spec (which a ``.pyi``-only
        directory yields), or when the target lies under uv's archive, which
        every environment later built from that entry would inherit.

    Raises:
        subprocess.CalledProcessError: ``uv cache dir`` failed.
    """
    spec = importlib.util.find_spec(package)
    if spec is None or spec.origin is None:
        return False
    cache = subprocess.run(["uv", "cache", "dir"], check=True, capture_output=True,
                           text=True).stdout.strip()
    return not target.resolve().is_relative_to(pathlib.Path(cache).resolve() / "archive-v0")


def _stub_root(target: pathlib.Path) -> pathlib.Path:
    """Directory a stub's files are laid out under: the tree, or a single module's parent."""
    return target if target.is_dir() else target.parent


@contextlib.contextmanager
def locked(target: pathlib.Path) -> Iterator[None]:
    """Hold an exclusive ``flock`` on a stub's directory, waiting until it is free.

    Args:
        target: a stub tree's directory, or a single-module ``.pyi``, whose
            parent is locked instead.
    """
    fd = os.open(_stub_root(target), os.O_RDONLY)
    try:
        fcntl.flock(fd, fcntl.LOCK_EX)
        yield
    finally:
        os.close(fd)


def remove_stale_temps(target: pathlib.Path) -> None:
    """Delete the temp files an interrupted ``write_back`` left beside a stub's files.

    Only safe under ``locked``, where no live run's files can match.

    Args:
        target: a stub tree's directory, or a single-module ``.pyi``.
    """
    pattern = f"{_TEMP_PREFIX}*{_TEMP_SUFFIX}"
    stale = target.rglob(pattern) if target.is_dir() else target.parent.glob(pattern)
    for path in stale:
        path.unlink()


def stage(target: pathlib.Path, dest: pathlib.Path) -> None:
    """Copy a stub's ``.pyi`` files to dest, leaving everything else behind.

    basedpyright resolves a submodule the tree omits from the environment, so
    nothing but the stubs needs to travel.

    Args:
        target: a stub tree's directory, or a single-module ``.pyi``.
        dest: directory to lay the copies out in, relative paths kept.
    """
    root = _stub_root(target)
    for path in target.rglob("*.pyi") if target.is_dir() else [target]:
        copy = dest / path.relative_to(root)
        copy.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, copy)


def write_back(staged: pathlib.Path, target: pathlib.Path) -> int:
    """Replace each changed file of a stub with its staged version.

    A file is replaced rather than written to, so its other hard links (conda's
    ``pkgs/`` cache, uv's hard-link installs) keep the old content. The marker
    file goes last, so its marker means the whole stub landed.

    Args:
        staged: directory holding the stub's ``.pyi`` files at their paths
            relative to its root.
        target: a stub tree's directory, or a single-module ``.pyi``.

    Returns:
        How many files differ in more than the marker line.
    """
    root = _stub_root(target)
    last = marker_file(target)
    changed = 0
    for copy in sorted(staged.rglob("*.pyi"), key=lambda p: root / p.relative_to(staged) == last):
        original = root / copy.relative_to(staged)
        new, old = copy.read_text(encoding="utf-8"), original.read_text(encoding="utf-8")
        if new == old:
            continue
        changed += _without_marker(new) != _without_marker(old)
        # Beside the original, so os.replace stays on one filesystem.
        fd, temp = tempfile.mkstemp(dir=original.parent, prefix=_TEMP_PREFIX, suffix=_TEMP_SUFFIX)
        try:
            with os.fdopen(fd, "w", encoding="utf-8") as f:
                f.write(new)
            shutil.copystat(original, temp)
            os.replace(temp, original)
        except BaseException:
            pathlib.Path(temp).unlink(missing_ok=True)
            raise
    return changed


def patch(target: pathlib.Path) -> dict[str, BaseException] | None:
    """Run the sequence on a stub.

    Args:
        target: a stub tree's directory, named after its package with or without
            ``-stubs``, or a single-module ``.pyi`` named after its module.

    Holding the stub's lock for more than ``WATCHDOG_S`` seconds ends the
    process, with every thread's traceback on stderr.

    Returns:
        {module: error} for the modules either pass skipped, or None when no
        module was skipped and at most the marker was written.

    Raises:
        SystemExit: the stub has no marker file, or the package does not import,
            which leaves the stub unwritten.
    """
    package = target.name.removesuffix("-stubs") if target.is_dir() else target.stem
    if not marker_file(target).is_file():
        raise SystemExit(f"no {marker_file(target)} to carry the marker")
    if not applicable(target, package):
        return None
    with locked(target):
        # Started under the lock, so waiting behind another editor's run does not
        # count. It fires from a C thread, which also ends an import stuck in C.
        faulthandler.dump_traceback_later(WATCHDOG_S, exit=True)
        try:
            if is_current(target, package):
                return None
            remove_stale_temps(target)
            with tempfile.TemporaryDirectory() as tmp:
                staged = pathlib.Path(tmp)
                skipped = _patch_staged(staged, target, package)
                changed = write_back(staged, target)
        finally:
            faulthandler.cancel_dump_traceback_later()
    return skipped if skipped or changed else None


def _patch_staged(staged: pathlib.Path, target: pathlib.Path,
                  package: str) -> dict[str, BaseException]:
    """Stage a stub, run the sequence on the copy and stamp its marker.

    Args:
        staged: empty directory to stage into.
        target: a stub tree's directory, or a single-module ``.pyi``.
        package: import name the stub describes.

    Returns:
        {module: error} for the modules either pass skipped.

    Raises:
        SystemExit: the package does not import, naming it and the error.
    """
    # Deferred: docify and libcst are only worth importing for a stub needing a run.
    import docify_stubs

    stage(target, staged)
    for path in staged.rglob("*.pyi"):
        path.write_text(_without_marker(path.read_text(encoding="utf-8")), encoding="utf-8")
    if target.is_dir():
        modules = [(module_name(p, staged, package), p) for p in introspection_order(staged)]
    else:
        modules = [(package, staged / target.name)]
    aliases = {key: runtime for key, runtime in CLASS_ALIASES.items()
               if key[0].split(".")[0] == package}
    repair = REPAIRS.get(package)
    # Some modules print on import, and stdout is where skipped modules are reported.
    with discarded_stdout():
        # Else every module would be skipped with the same error.
        try:
            importlib.import_module(package)
        except (Exception, SystemExit) as e:
            raise SystemExit(f"{package} does not import: {e}") from e
        skipped = docify_stubs.document(modules, side_effect_free, aliases)
        if repair:
            skipped |= repair(staged, package)
    marked = staged / marker_file(target).relative_to(_stub_root(target))
    marked.write_text(_marker(package) + marked.read_text(encoding="utf-8"), encoding="utf-8")
    return skipped


def main() -> None:
    skipped = patch(pathlib.Path(sys.argv[1]).absolute())
    if skipped is None:
        sys.exit(NOTHING_WRITTEN)
    for module, e in skipped.items():
        print(f"{module} ({describe(e)})")


def describe(error: BaseException) -> str:
    """One line naming an error's type and what it says.

    A parser error gives its message and the line it points at. Any other error
    gives its message's first line, since some span several.

    Args:
        error: the error.

    Returns:
        ``<type>: <message>``, or the type alone when the message is empty.
    """
    # Looked up, not imported, which would cost a run that finds the stub
    # current: a libcst error exists only once libcst is loaded.
    libcst = sys.modules.get("libcst")
    if isinstance(error, SyntaxError) and error.lineno:
        message = f"{error.msg} (line {error.lineno})"
    elif isinstance(error, tokenize.TokenError):
        text, (line, _) = error.args
        message = f"{text} (line {line})"
    elif libcst and isinstance(error, libcst.ParserSyntaxError):
        message = f"{_LIBCST_ERROR_PREFIX.sub('', error.message)} (line {error.editor_line})"
    else:
        message = str(error).partition("\n")[0]
    name = type(error).__name__
    return f"{name}: {message}" if message else name


if __name__ == "__main__":
    main()
