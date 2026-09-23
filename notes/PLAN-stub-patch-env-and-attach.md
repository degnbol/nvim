# Plan: stub patching — interpreter choice and attach-path test
Status: spec

## Relevant files

- `lua/autocmds/stub_patches.lua` — header, `env_of`, `locate`, the
  `LspAttach` handler.
- `tests/plenary/stub_patches_spec.lua` — `with`, `reload`, `fake_env`,
  `describe("stub_patches.env_of")`.
- Skills: `lsp` (this config's `lsp_ext/`), `neovim` → `references/lsp.md`
  (in-process servers), `references/testing.md`.

## Background

- **Interpreter guess.** A prefix can hold several Pythons (Homebrew:
  `lib/python3.12/site-packages` beside a `bin/python3` that is 3.13). The
  guessed interpreter cannot see the package, so the run exits 3 silently.
  `bin/python<X.Y>`, with `X.Y` from the site-packages path, exists in venvs,
  uv venvs, conda and Homebrew, and covers free-threaded `3.13t`.
- **Attach path untested.** The integration spec calls `M.consider` itself, so
  `LspAttach` → `M.imports` → `locate` → the `M.state` deduplication has no
  test. A real basedpyright cannot drive it without racing
  `assert_documents`'s pre-patch check; an in-process fake server answering
  `textDocument/declaration` can, deterministically.
- **Stale header.** The header says the module exists for projects owning a
  basedpyright config. It patches an environment's own stubs for any project.

## Changes

### `lua/autocmds/stub_patches.lua`

- `env_of`: capture `X.Y` with
  `^(.*)(/lib/python([^/]*)/site%-packages)/([^/]+)`. `python` is the first of
  `bin/python<X.Y>` (when `X.Y` is non-empty), `bin/python` that exists, else
  `bin/python3` (kept even when missing, so `consider` reports it).
  `StubPatch.Env.python`'s description is unchanged.
- Header: the module patches an environment's own stubs for any project, since
  a config-less project reads them too for a package with no vendored copy. The
  config-owning case stays as the reason `stubPath` alone is not enough.

### `tests/plenary/stub_patches_spec.lua`

- `env_of`: an env holding `bin/python3.12`, `bin/python` and `bin/python3`
  resolves to `bin/python3.12` for a `lib/python3.12` path, and one holding
  `bin/python3.13t` and `bin/python` to `bin/python3.13t` for
  `lib/python3.13t`.
- `fake_env`: lay the tree out as `root/lib/python3.12/site-packages/pkg` beside
  `root/bin/python`, the layout `env_of` matches. The spawn-gate and enqueue
  tests do not depend on the layout.
- New `describe("stub_patches attach")` with `before_each(reload)`, which also
  re-creates the `my.stub_patches` augroup the integration test deletes:
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

## Expected outcome

Homebrew and other multi-Python prefixes are patched. The attach path and its
deduplication are covered by a test.

## Non-goals

- Driving the real-basedpyright integration test through `LspAttach`.

## Open questions

None.
