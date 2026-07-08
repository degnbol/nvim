---
name: lsp
description: This config's LSP setup — lsp_ext/ external type sources and stubs for basedpyright, the R languageserver `...`-forwarding patch, and the mason-lspconfig field-override workaround. Use when editing lsp/*.lua, lua/plugins/lsp.lua, files under lsp_ext/, or configuring basedpyright/ruff/r_language_server/tinymist.
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

## `lsp_ext/` — external sources and stubs

Extra type information for basedpyright, shared between neovim and the Claude
lint hook.

```
lsp_ext/
├── extraPaths/             # Source directories added to pyright's extraPaths
│   ├── kitty-source/       # Kitty source (git submodule) — provides kitty.* types
│   └── pymol_modules/      # Symlink → ../pymol-open-source/modules/
├── python_stubs/           # .pyi stub files (pyright stubPath)
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
