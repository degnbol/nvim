---
name: lsp
description: This config's LSP setup — lsp_ext/ external type sources and stubs for basedpyright, per-script clients for PEP 723 uv scripts, the R languageserver `...`-forwarding patch, and the mason-lspconfig field-override workaround. Use when editing lsp/*.lua, lua/plugins/lsp.lua, files under lsp_ext/, or configuring basedpyright/ruff/r_language_server/tinymist.
---

# LSP (this config)

General LSP knowledge (in-process servers, external `lsp/<name>.lua` config,
`vim.lsp.enable`, the "use an in-process server, not custom completefunc"
strategy) is in the global neovim skill's
`~/.claude/skills/neovim/references/lsp.md`. This file is the config-specific
setup only.

In-process LSP working examples in this config: `modules/kitty-conf.nvim`
(hover + completion), `modules/agentic.nvim/lua/agentic/completion/lsp_server.lua`
(trigger-character completion for `/` and `@`).

Inspecting a server's replies: [references/probing.md](references/probing.md).

## `lsp_ext/` — external sources and stubs

Extra type information for basedpyright, shared between neovim and the Claude
lint hook.

```
lsp_ext/
├── extraPaths/             # Source directories added to pyright's extraPaths
│   ├── kitty-source/       # Kitty source (git submodule) — provides kitty.* types
│   └── pymol_modules/      # Symlink → ../pymol-open-source/modules/
├── python_stubs/           # .pyi stub files (pyright stubPath)
├── stub_patches/           # Docstring pass + rdkit repair, run on stub trees
├── pymol-open-source/      # Full pymol source repo (git submodule)
└── r_lsp_dots.R            # R languageserver monkey-patch (see below)
```

**How it works:** `lsp/basedpyright.lua` globs `lsp_ext/extraPaths/*/` for
import resolution paths. The Claude `lib/lint-tier.sh` library uses the same
glob when generating a fallback pyright config for projects without their own
`pyrightconfig.json`.

**Adding a new source:** Drop the directory in `extraPaths/` (or symlink it
there). Both neovim and the lint hook pick it up automatically — no config
changes needed.

**`stubPath` reaches config-less projects only.** basedpyright drops the
language server's config settings (`stubPath`, `extraPaths`, …) for a project
owning a `pyrightconfig.json` (even `{}`) or a `[tool.basedpyright]` section.
A `pyproject.toml` without that section keeps them. `pythonPath` survives.
Such projects read only their environment's own stubs.
`lua/autocmds/stub_patches.lua` patches those in place without a prompt: on
`LspAttach` it asks for the declaration of every imported package and queues
`lsp_ext/stub_patches/patch_stubs.py` on the tree it lands in, then sends
`didChangeWatchedFiles`. The script decides what applies (see
`lsp_ext/stub_patches/README.md`).

## PEP 723 uv scripts — a client per script

A `# /// script` file's dependencies live in an environment of its own under
`~/.cache/uv/environments-v2/`. `lua/autocmds/uv_script_env.lua` gives such a
buffer its own basedpyright with `settings.python.pythonPath` set to that
interpreter, entered from `root_dir` in `lsp/basedpyright.lua` and stopped on
its last detach. One server per open uv script.

- `pythonPath` is not a config-file field, so unlike `stubPath`/`extraPaths`
  (above) it survives a project's own `pyrightconfig.json`.
- `reuse_client = M.same_env` must stay a **field** of the config, not only an
  argument to `vim.lsp.start`: neovim's default compares name and workspace
  folders only, so a plain `.py` buffer rooted at a script's directory would
  inherit that script's environment. Which client belongs to which environment
  is carried by a `uv_script_python` field on the config, not by
  `settings.python.pythonPath` — `LspPyrightSetPythonPath` rewrites the live
  settings in place.
- The name stays `basedpyright`, so `stub_patches` applies, once per script
  environment — and again after a dependency-list edit, which is a new
  environment.

## R language server — `...` forwarding patch

The R languageserver doesn't resolve `...` forwarding — functions like
`scale_y_log10(...)` that delegate to `scale_y_continuous(...)` only show `...`
as a parameter, with no completion for the actual arguments.

**Patch:** `lsp_ext/r_lsp_dots.R` monkey-patches the languageserver at startup
via `setHook(packageEvent("languageserver", "onLoad"), ...)`. Two patches:
1. **`get_formals`** (R6 `$set` on `PackageNamespace`) — when a function has
   only `...` as its formal, traces the body to find the target function and
   returns its formals instead.
2. **`arg_completion`** (`assignInNamespace`) — wraps the original to add
   `textEdit` (for reliable `" = "` insertion) and fix `data$funct` (so
   `completionItem/resolve` looks up docs for the underlying function, not the
   wrapper).

**Config:** `lsp/r_language_server.lua` sets a custom `cmd` that sources the
patch before `languageserver::run()`.

## mason-lspconfig field override

Mason's `automatic_enable` calls `vim.lsp.config()` which overrides fields from
`lsp/*.lua` — not just `cmd` but also `filetypes` and other fields. Fix: in
`lua/plugins/lsp.lua`, re-apply custom config after `setup()`. Already done for:
`r_language_server` (custom cmd), `basedpyright` and `ruff` (compound
filetypes), `tinymist` (`on_attach`, which pins a main so ref hover/goto-def
work — see the comment there).
