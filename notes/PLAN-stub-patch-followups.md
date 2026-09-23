# Plan: stub patching follow-ups from review
Status: spec

## Relevant files

- `lsp_ext/stub_patches/patch_stubs.py` — `distribution_version` (`:97`),
  `stage` docstring, `write_back` (`:244`, `copystat`), `_patch_staged`
  (`redirect_stdout`, `:310`).
- `lsp_ext/stub_patches/docify_stubs.py` — `patch_docify` (`:36-39`),
  `get_qualname` (`:47-72`).
- `lsp_ext/stub_patches/README.md`, `test_patch_stubs.py` (`site` fixture
  `:29-32`, rdkit test).
- `lua/autocmds/stub_patches.lua` — header (`:1-10`), `env_of` (`:60-66`),
  `run_next` (`:107-125`).
- `tests/plenary/stub_patches_spec.lua` — integration test (`:219`, `:239`),
  `tests/README.md` (test command).
- `lsp_ext/python_stubs/README.md` (`:4`), `notes/PLAN-rdkit-stub-typing.md`
  (`:243`, `:247`, `:274`).

## Background

Findings from reviewing the stub-docstring implementation, each small and
separate from the per-module failure handling.

- **Python floor.** `importlib.metadata.packages_distributions` is 3.10+, so a
  3.9 environment fails every session with AttributeError. 3.9 is past end of
  life (October 2025), so the floor becomes 3.10 rather than gaining a fallback.
- **mtime.** `shutil.copystat` gives rewritten content the old mtime, so caches
  keyed on mtime (parso, hence jedi-based tools) keep the pre-patch parse.
  basedpyright is unaffected (it gets `didChangeWatchedFiles`). Permissions are
  what needs copying.
- **C-level prints.** `redirect_stdout` does not catch output a C extension
  writes to file descriptor 1 on import. That text lands on stdout, which is
  the skipped-module report.
- **Interpreter guess.** A prefix can hold several Pythons (Homebrew:
  `lib/python3.12/site-packages` beside a `bin/python3` that is 3.13). The
  guessed interpreter cannot see the package, so the run exits 3 silently.
  `bin/python<X.Y>`, with `X.Y` from the site-packages path, exists in venvs,
  uv venvs, conda and Homebrew, and covers free-threaded `3.13t`.
- **Qualname scopes.** `get_qualname` steps over every non-class scope. For a
  `FunctionScope` that resolves a local name to a same-named module attribute.
  Only libcst's `AnnotationScope` needs stepping over.
- **Pin check.** `assert` vanishes under `python -O`.
- **Queue stall.** If `vim.system` fails to spawn, `running` stays true and the
  queue stops for the session.
- **Test gaps.** The integration spec calls `M.consider` itself, so LspAttach →
  `M.imports` → `locate` → deduplication is untested. The rdkit test checks
  constructors only. The `site` fixture evicts libcst and docify after every
  test, so each test re-imports them.
- **Stale prose.** Listed under Documentation.

## Changes

### `patch_stubs.py`

- `write_back`: `shutil.copymode` instead of `shutil.copystat`.
- `_patch_staged`: silence introspection at the descriptor level. A context
  manager that `os.dup`s fd 1, `os.dup2`s a devnull fd over it, and restores it
  in `finally`, replacing `redirect_stdout`:

  ```python
  @contextlib.contextmanager
  def stdout_silenced() -> Iterator[None]:
      """Send everything written to stdout, by Python or by C, to devnull."""
  ```

  Flush `sys.stdout` before and after, so buffered Python output is not
  reordered.
- `stage` docstring: only `.pyi` files are copied because only `.pyi` files are
  patched and written back. Drop the basedpyright rationale.

### `docify_stubs.py`

- `patch_docify`: raise `RuntimeError` instead of `assert` for a missing pin or
  a version mismatch. The `Raises:` section follows.
- `get_qualname`: step over `libcst.metadata.scope_provider.AnnotationScope`
  only, and raise `TypeError` for any other scope, as docify does. Drop the
  nested-`def` caveat from the docstring.

### `lua/autocmds/stub_patches.lua`

- `env_of`: capture `X.Y` from `lib/python([^/]*)/site-packages` and use
  `bin/python<X.Y>` when it exists, falling back to `bin/python` then
  `bin/python3`. `StubPatch.Env.python`'s description is unchanged.
- `run_next`: wrap `vim.system` in `pcall`. On failure set `failed`, notify
  ERROR with the message, and call `run_next` so the queue continues.
- Header: say the module patches an environment's own stubs for any project,
  since a config-less project reads them too for packages with no vendored copy.
  The config-owning case stays as the reason `stubPath` alone is not enough.

### Python floor

`lsp_ext/stub_patches/README.md` states Python ≥ 3.10. `patch_stubs.main` exits
with a message naming the floor when `sys.version_info < (3, 10)`, before any
work, so a 3.9 environment gets one clear ERROR rather than a traceback.

### Tests

- `stub_patches_spec.lua`: keep the `my.stub_patches` augroup in the
  integration test and assert the state reaches `patched` from the attach
  alone. Attach a second buffer importing the same package and assert one run.
  Add a gate-level test for a spawn failure (stub `vim.system` to error) that
  checks the next queued env still runs.
- `stub_patches_spec.lua` `env_of`: an env with `bin/python3.12` beside
  `bin/python3` resolves to `python3.12` for a `lib/python3.12` path.
- `test_docify_stubs.py`: a local annotated name inside a function body is not
  documented from a same-named module attribute.
- `test_patch_stubs.py`: the `site` fixture evicts only modules whose
  `__file__` lies under the fixture root. The rdkit test also asserts the
  count of `C++ signature` in the sequenced tree is at most the count in the
  tree repaired alone. A write-back test asserts the rewritten file's mtime is
  newer than before.

### Documentation

- `lsp_ext/python_stubs/README.md:4`: every tree but pyrosetta goes through
  `patch_stubs.py`.
- `tests/README.md`: add `-i NONE` to both plenary commands, so a test run does
  not write the user's shada.
- `notes/PLAN-rdkit-stub-typing.md`: `:243` and `:247` name
  `patch_pybind_stubs.py --out`, which no longer exists (vendoring is
  `patch_pybind_stubs.vendor` plus `patch_stubs.py`). `:274` cites
  `stub_patches.lua:120`, the removed marker check. Rewrite those passages to
  the current mechanism.

## Expected outcome

Tools caching by mtime see patched stubs. C-extension import noise no longer
appears in the report. Homebrew and other multi-Python prefixes are patched. A
3.9 environment gets one clear message. A failed spawn does not stop later
runs. The attach path is covered by a test.

## Non-goals

- A Python 3.9 fallback for `packages_distributions`.
- Patching pyrosetta's vendored tree (its size is to be measured first).

## Open questions

None.
