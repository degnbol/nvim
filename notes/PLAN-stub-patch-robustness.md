# Plan: a badly written package is reported per module, never fatal or silent
Status: spec

## Relevant files

- `lsp_ext/stub_patches/patch_pybind_stubs.py` — `patch_dir` (`:589-616`); the
  `ast.parse` calls it reaches: `keepable_signature_lines` (`:139`, via `clean`),
  `toplevel_names` (`:311`), `substitute_class` (`:509`).
- `lsp_ext/stub_patches/docify_stubs.py` — `document` (`:171-210`).
- `lsp_ext/stub_patches/patch_stubs.py` — `REPAIRS` docstring (`:39-40`),
  `write_back` (`:238`), `patch` (`:282`), `_patch_staged` (`:302-315`), `main`,
  module docstring (the stdout contract).
- `lua/autocmds/stub_patches.lua` — `finish` (`:90-104`).
- `lsp_ext/stub_patches/test_patch_pybind_stubs.py`, `test_docify_stubs.py`,
  `test_patch_stubs.py`, `README.md`.
- docify 1.2.1 `docify.py`: `run_one` (`:535-597`), `logger.exception("could not
  parse ...")` (`:563`).
- Trigger file: rdkit 2026.03.6 `rdkit-stubs/Chem/rdSynthonSpaceSearch.pyi`.

## Background

The sequence runs on whatever stubs a package ships, so a badly written stub
must cost only its own module. Today it costs the tree, or goes unreported.

**Fatal.** rdkit 2026.03.6's `Chem/rdSynthonSpaceSearch.pyi` declares a
non-default parameter after a default one. `ast.parse` raises in
`keepable_signature_lines`, `patch_dir` propagates it, and `patch_stubs.py`
exits 1: an ERROR every session in every rdkit environment at that version, and
a failed `RUNME.sh`, which does not pin rdkit. libcst rejects the file too
(`ParserSyntaxError @ 343:12`), so docify skips it as well. A stub that is not
UTF-8 fails the same way, at staging.

**Silent.** Two paths end in exit 3, state `ok`, and no message:

- docify's `run_one` catches its own parse failure and only logs it. Nothing
  configures docify's logger, so the text reaches stderr through `lastResort`,
  and the editor discards stderr on success.
- A package whose root fails to import (a missing shared library, a broken
  environment) gets no docstrings, yet the marker is stamped. `write_back`
  counts no content change, `patch` returns None, and the tree counts as
  current until the version or a script changes. Unimportable submodules are
  likewise unreported whenever nothing else in the tree changed.

The rule becomes: a module a pass cannot process keeps what earlier passes gave
it, and is reported with the reason, whether or not anything else changed. A
tree whose package root could not be imported is not marked current.

**Only faults in the stub itself are isolated**: `SyntaxError` (with its
`IndentationError`), `tokenize.TokenError`, `libcst.ParserSyntaxError` and
`UnicodeDecodeError`. Everything else still fails the run: a bug in these
scripts, a drifted `_CLASS_OVERRIDES` entry (whose `KeyError` is meant to fail
loudly), staging, the lock, write-back and the docify pin check. A defect of
the sequence reported as one line per module would read as a package problem,
and be silenced by the marker after the first session.

## Changes

### New in `lsp_ext/stub_patches/stub_tree.py`

```python
STUB_FAULTS: tuple[type[Exception], ...] = (SyntaxError, tokenize.TokenError,
                                             UnicodeDecodeError)
"""Exceptions that mean a stub file itself is malformed, as opposed to a defect
of the patch scripts. libcst's ParserSyntaxError is added by docify_stubs, since
this module imports no third-party package."""

def skip_reason(module: str, error: BaseException) -> str:
    """One report line for a module a pass skipped.

    Args:
        module: dotted module name.
        error: why it was skipped.

    Returns:
        ``<module> (<ExcType>: <first line of the message>)``, whitespace
        collapsed, so each line of the report names one module.
    """
```

Every report line in both passes goes through `skip_reason`, including the
existing import-failure lines (`f"{name} ({e})"` today), so stdout stays one
module per line: `str(ParserSyntaxError)` and some import errors span several.

### `patch_pybind_stubs.patch_dir`

Move the file read into a `try` around the read and `clean`, catching
`STUB_FAULTS`. On a fault, report `skip_reason(name, e)` and `continue`, which
leaves the staged file as it was (docified). The runtime passes stay outside
that `try`, so a bug in them or an override `KeyError` still fails the run.
The import failure keeps its own `except (Exception, SystemExit)`. The write is
already the loop's last statement, so a skipped file is never half-written.

### `docify_stubs.document`

- Around each module's body (`read_text`, `comment_keyword_targets`, alias
  binding, `run_one`), catch `STUB_FAULTS + (libcst.ParserSyntaxError,)`, report
  `skip_reason`, and continue. A fault in a module `safe` rejects is reported
  the same way, since commenting keyword targets runs on every module.
- Report docify's own parse failure. After `patch_docify`, attach a
  `logging.handlers.BufferingHandler(capacity=10**6)` at level `ERROR` to
  `docify.logger`, and remove it in a `finally`. After each `run_one`, for each
  buffered record append `skip_reason(name, record.exc_info[1])` when
  `exc_info` is set, else `f"{name} ({first line of record.getMessage()})"`,
  then `flush()`. `getMessage()` alone names only the staged temp path, which
  is gone by the time anyone reads it. Records below `ERROR` (per-attribute
  warnings that skip nothing) are no longer printed, which loses nothing: the
  editor discarded them already.

### `patch_stubs.py`

- `stage`, `write_back` and the marker stamp read and write with
  `errors="surrogateescape"`, so a stub with bytes that are not UTF-8 travels
  through staging unchanged. The passes then raise `UnicodeDecodeError` on it
  (they read strictly), which reports that module.
- Merge the two passes' lines by module: `_patch_staged` returns the lines in
  order, and `patch` keeps the first line per module name, so a file both
  parsers reject is reported once.
- `patch` reports skipped modules whether or not anything changed: it returns
  the list when it is non-empty, or when anything changed, and None only when
  neither holds.
- `_patch_staged` does not stamp the marker when the package root itself (the
  first module, `package`) is among the skipped. The docstrings still land, and
  the next session retries the tree.
- Rename `unimportable` to `skipped` throughout, including the `REPAIRS`
  docstring and `_patch_staged`'s Returns line: "the modules a pass skipped,
  one report line each".
- Module docstring: exit 0 "patched, or something to report, with the modules
  a pass skipped on stdout, one per line".

### `lua/autocmds/stub_patches.lua`

`finish`: the heading "introspection skipped (import failed):" becomes
"skipped:", since each line carries its reason. Exit 0 with nothing written but
a report still sets `patched` and notifies, which is correct: the message is
the report.

### Tests

`test_stub_tree.py`: `skip_reason` collapses a multi-line message to its first
line and names the exception type.

`test_patch_pybind_stubs.py`:
- A tree with one module that does not parse (a non-default parameter after a
  default one) and one that does. `patch_dir` returns one `SyntaxError` line for
  the first, leaves it byte-identical, and repairs the second.
- An override `KeyError` still raises out of `patch_dir`.

`test_docify_stubs.py`:
- A stub libcst cannot parse: `document` returns one line naming the module
  (not a path) with the libcst error, and documents the other modules.
  `docify.logger.handlers` is as before the call.
- A stub that does not tokenize (an unterminated triple-quoted string) is
  reported, and the next module is still documented.

`test_patch_stubs.py`:
- A tree with one unparseable module: `patch` returns exactly one line for it,
  and the marker is written.
- A package whose root raises `ImportError`: `patch` returns its line, the
  marker is not written, and a second `patch` runs again.
- A stub with a Latin-1 byte: reported, and the rest of the tree patched.
- `test_rdkit_sequence_matches_the_repair_alone` runs on whatever rdkit is
  installed. Its "does not parse" comment and assertion hold once the repair
  skips such files. Until now they have only run on 2026.3.5, which has none.
  Drop the `rdkit==2026.3.5` pin from `README.md`'s test command.

### Documentation

`lsp_ext/stub_patches/README.md`, "The sequence": a module a pass cannot parse
or decode keeps what earlier passes wrote and is reported with its reason, and
a tree whose package does not import is retried next session.

## Expected outcome

On rdkit 2026.03.6, opening a buffer patches the environment's `rdkit-stubs`
and notifies once, listing `rdkit.Chem.rdSynthonSpaceSearch (SyntaxError: ...)`.
That module keeps bundled stubs unchanged, and every other module is documented
and repaired. `./RUNME.sh` completes on current rdkit. A package that fails to
import, or a stub docify cannot parse, now appears in a notification. A drifted
override or a bug in the scripts still fails the run with an ERROR.

## Non-goals

- **Repairing the unparseable rdkit stub.** An upstream defect, to report
  upstream.
- **Retrying a skipped submodule while the root imports.** The marker covers the
  tree, so such a module stays skipped until the package version or a script
  changes. Only an unimportable root blocks the marker.
- **Reporting docify's per-attribute warnings.**

## Open questions

None.
