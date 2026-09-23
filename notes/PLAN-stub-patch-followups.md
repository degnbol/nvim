# Plan: stub patching follow-ups from review
Status: spec

## Relevant files

- `lsp_ext/stub_patches/patch_stubs.py` — `main`, `write_back` (`shutil.copystat`),
  `stage` docstring.
- `lsp_ext/stub_patches/docify_stubs.py` — `patch_docify`, `get_qualname`,
  `document` (catches only `cst.ParserSyntaxError` around `with_docstrings`).
- `lsp_ext/stub_patches/test_patch_stubs.py` — `site` fixture, rdkit sequence
  test. `test_docify_stubs.py`, `test_patch_pybind_stubs.py` (slow tests).
- `lsp_ext/stub_patches/README.md`, `lsp_ext/python_stubs/README.md`,
  `lsp_ext/python_stubs/RUNME.sh`.
- `lua/autocmds/stub_patches.lua` — header, `env_of`, `run_next`, `locate`, the
  `LspAttach` handler.
- `tests/plenary/stub_patches_spec.lua` — `with`, `fake_env`, spawn-gate tests.
- docify 1.2.1 `docify.py` — `get_obj` (returns None on `AttributeError`),
  `get_qualname` (raises `TypeError` on a non-class scope).
- Skills: `lsp` (this config's `lsp_ext/`), `neovim` → `references/lsp.md`
  (in-process servers), `references/testing.md`.

## Background

Findings from reviewing the stub-docstring implementation, each small and
independent.

- **Python floor.** `importlib.metadata.packages_distributions` is 3.10+, and
  `is_current` reaches it through `_marker`, so a 3.9 environment fails every
  session with a traceback. 3.9 is past end of life (October 2025): the floor
  becomes 3.10, stated once, instead of a fallback.
- **mtime.** `shutil.copystat` gives rewritten content the old mtime, so caches
  keyed on mtime (parso, hence jedi-based tools) keep the pre-patch parse.
  basedpyright is unaffected (it gets `didChangeWatchedFiles`). Only the mode
  needs copying; dropping BSD flags and xattrs with it is intended.
- **Interpreter guess.** A prefix can hold several Pythons (Homebrew:
  `lib/python3.12/site-packages` beside a `bin/python3` that is 3.13). The
  guessed interpreter cannot see the package, so the run exits 3 silently.
  `bin/python<X.Y>`, with `X.Y` from the site-packages path, exists in venvs,
  uv venvs, conda and Homebrew, and covers free-threaded `3.13t`.
- **Qualname scopes.** `get_qualname` steps over every scope that is neither
  global nor class, so a name declared inside a function body (`def f(): x:
  int`) gets qualname `x` and is documented from a same-named module attribute.
  Only libcst's `AnnotationScope` (PEP 695) is safe to step over. Raising for
  the rest, as docify does, would end the run: nothing between docify's
  `Transformer` and `document` catches it. Python's own answer is `__qualname__`
  syntax, `f.<locals>.x`, which `docify.get_obj` fails to resolve, returning
  None, so nothing is documented and nothing raises.
- **Pin check.** `assert` vanishes under `python -O`, and an `AssertionError`
  reaches the editor as a traceback. The run's other fatal setup errors are
  `SystemExit(message)`, one stderr line.
- **Queue stall.** `running` is set before `vim.system` spawns. If the spawn
  throws (missing cmd or cwd: ENOENT), or `finish` throws in the exit callback,
  `run_next` is never called again and the queue stops for the session.
- **Attach path untested.** The integration spec calls `M.consider` itself, so
  `LspAttach` → `M.imports` → `locate` → the `M.state` deduplication has no
  test. A real basedpyright cannot drive it without racing
  `assert_documents`'s pre-patch check; an in-process fake server answering
  `textDocument/declaration` can, deterministically.
- **Test gaps.** The rdkit sequence test does not compare `C++ signature`
  counts (98 in each tree on rdkit 2026.3.6, none higher per file in the
  sequenced one). The `site` fixture in `test_patch_stubs.py` evicts every module first
  imported during a test, including `docify_stubs` and its libcst submodules
  (loaded lazily by `_patch_staged`), so each test re-imports them.
  `test_docify_stubs.py` imports `docify_stubs` at the top and does not pay
  this. The declarations `declare_missing_names` appends (synthesised `def`s
  and re-exports) were measured never to shadow a name basedpyright already
  resolves, a property of the rdkit version rather than of the guard; nothing
  re-checks it.
- **Stale prose.** Listed under Documentation.

## Changes

### `patch_stubs.py`

- `PYTHON_FLOOR = (3, 10)`, module level. `main` begins with
  `if sys.version_info < PYTHON_FLOOR: sys.exit(...)`, before any other work
  (including reading `sys.argv`). The message names the floor and the found
  version, formatted by indexing (`sys.version_info[:2]`), so a plain tuple
  works too.
- `write_back`: `shutil.copymode` instead of `shutil.copystat`.
- `stage` docstring: only `.pyi` files are copied because only `.pyi` files are
  patched and written back. Drop the basedpyright rationale.

### `docify_stubs.py`

- `patch_docify`: raise `SystemExit(message)` instead of `assert` for a missing
  pin or a version mismatch. `Raises:` follows.
- `get_qualname(scope: meta.Scope, name: str) -> str`: `ClassScope` prepends
  `<name>.`; `FunctionScope` (a `def`'s scope) prepends `<name>.<locals>.`;
  either raises `ValueError` when nameless (a class as now; a lambda, which no
  declaration can be inside). Any other scope (`AnnotationScope`, a
  comprehension) is stepped over. This is libcst's own `<locals>` rule.
  Docstring: the result is Python's `__qualname__` for the declaration, which
  `docify.get_obj` resolves only for class and module members. Drop the
  nested-`def` caveat.

### `lua/autocmds/stub_patches.lua`

- `env_of`: capture `X.Y` with
  `^(.*)(/lib/python([^/]*)/site%-packages)/([^/]+)`. `python` is the first of
  `bin/python<X.Y>` (when `X.Y` is non-empty), `bin/python` that exists, else
  `bin/python3` (kept even when missing, so `consider` reports it).
  `StubPatch.Env.python`'s description is unchanged.
- `run_next`: `pcall` the `vim.system` call. On failure set the state `failed`,
  notify ERROR with the error, and call `run_next`. In the exit callback,
  `pcall(finish, env, out)`; on failure notify ERROR, then `run_next` either
  way. A throwing `finish` leaves whatever state it set first (every branch of
  `finish` sets the state before anything that can throw).
- Header: the module patches an environment's own stubs for any project, since
  a config-less project reads them too for a package with no vendored copy. The
  config-owning case stays as the reason `stubPath` alone is not enough.

### Tests

- `tests/plenary/stub_patches_spec.lua`:
  - `fake_env`: lay the tree out as `root/lib/python3.12/site-packages/pkg` beside
    `root/bin/python`, the layout `env_of` matches. The spawn-gate tests do not
    depend on the layout.
  - New `describe("stub_patches attach")`, placed before
    `describe("stub_patches")`, whose integration test deletes the
    `my.stub_patches` augroup (say so in a comment):
    - `it("queues one run per tree, however many buffers import it")`. Replace
      `vim.system` with a recorder (`with`) and `vim.fn.executable` with one
      returning 1. Start an in-process client named `basedpyright`
      (`vim.lsp.start` with `cmd` a function returning a fake server, and a
      shared `root_dir`) whose `initialize` reply advertises
      `declarationProvider`, whose `textDocument/declaration` reply is a
      `Location` in the `fake_env()`'s `__init__.pyi`, and which replies
      `callback(nil, nil)` to anything else (other `LspAttach` handlers may
      ask). Attach two named buffers, each `import pkg`, the second reusing the
      client. Wrap `stub_patches.env_of` to count calls; `vim.wait` until it has
      run twice (both replies handled), then assert one recorded run and the
      state `patching`. Finish the run with `code = 3`, `vim.wait` until the
      state is not `patching`, stop the client and delete both buffers, so no
      run or client outlives the test.
    - `it("keeps the queue going after a spawn failure")`: `vim.system` throws
      for the first env and records the second; enqueue both; assert the first
      is `failed`, one ERROR was notified, and the second is `patching`. Finish
      the second run with `code = 3` and wait as above.
    - `it("keeps the queue going when finishing a run throws")`: wrap
      `stub_patches.notify_changed` to throw; enqueue two envs; finish the first
      with `code = 0`; assert it is `patched`, one ERROR was notified, and the
      second is `patching`. Finish the second with `code = 3` and wait.
  - `env_of`: an env holding `bin/python3.12`, `bin/python` and `bin/python3`
    resolves to `bin/python3.12` for a `lib/python3.12` path, and one holding
    `bin/python3.13t` and `bin/python` to `bin/python3.13t` for
    `lib/python3.13t`.
- `test_docify_stubs.py`, `test_a_function_local_name_is_not_documented`: the
  stub `def f():\n    sqrt: float\n` over the existing `RUNTIME` (module-level
  `sqrt = math.sqrt`) comes back unchanged, and `document` returns `{}`. It
  fails today: the local is documented from `math.sqrt`.
- `test_patch_stubs.py`:
  - `import docify_stubs` at the top, so the `site` fixture's snapshot holds it
    and its libcst submodules.
  - `test_an_old_python_exits_naming_the_floor`: monkeypatch `sys.version_info`
    to `(3, 9, 0)`, call `patch_stubs.main`, expect `SystemExit` whose message
    names `3.10`.
  - `test_write_back_gives_rewritten_files_a_new_mtime`: backdate
    `__init__.pyi` with `os.utime(path, (0, 0))`, patch, assert its mtime is
    greater than 0 and its mode is unchanged.
  - The rdkit sequence test also asserts, per file, that the count of
    `C++ signature` in the sequenced tree is at most the count in the tree
    repaired alone.
- `test_patch_pybind_stubs.py`,
  `slow_test_appended_names_are_unresolved_without_them`: copy
  `../python_stubs/rdkit` to a temp dir with every file cut at `_ADDED_HEADER`.
  Probe `from <module> import <name>` for every name bound below that header in
  the vendored tree (`toplevel_names` of the text after it). Type the probe as
  `slow_test_pyright_reads_the_vendored_tree_as_intended` does
  (`pyrightconfig.json` with `stubPath` at the temp copy,
  `basedpyright --outputjson`), and assert every probe line has an error
  `"<name>" is unknown import symbol`.

### Documentation

- `lsp_ext/python_stubs/README.md`: `RUNME.sh` runs every tree but pyrosetta
  through `patch_stubs.py`.
- `lsp_ext/stub_patches/README.md`: Python ≥ 3.10 for the patched environment.

## Expected outcome

Tools caching by mtime see patched stubs. Homebrew and other multi-Python
prefixes are patched. A 3.9 environment, or a mismatched docify, gets one clear
ERROR line. A declaration inside a function body is not documented and does not
fail the run. A failed spawn or a throwing `finish` does not stop later runs.
The attach path and its deduplication are covered by a test.

## Non-goals

- A Python 3.9 fallback for `packages_distributions`.
- Patching pyrosetta's vendored tree (its size is to be measured first).
- Driving the real-basedpyright integration test through `LspAttach`.

## Open questions

None.
