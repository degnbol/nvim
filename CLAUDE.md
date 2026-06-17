# Neovim Config

`~/nvim` = `~/.config/nvim` = `~/dotfiles/config/nvim` (symlinks). All the same directory.

## Shared utilities (`lua/utils/`)

Before reaching for `vim.api.nvim_replace_termcodes`, visual-range extraction, keymap wrappers, ColorScheme hooks, etc., check `lua/utils/init.lua`, `lua/utils/keymap.lua`, and `lua/utils/highlights.lua` — most common patterns are already wrapped (e.g. `nvim_code(s)` for `<Esc>`/`<CR>` termcode conversion, `last_visual_range()`, `get_cursor`/`set_cursor`, `press(keys)`, `hi.onColorScheme(fn)` to register a ColorScheme autocmd that also fires once on load). When you write a helper that's reusable across files, add it to `utils/` rather than duplicating inline.

## Debugging

Do not conclude that a bug is in neovim core until confirmed with `nvim --clean` - it is usually a plugin problem.

## Plugin Management (vim.pack + lz.n)

Plugins are managed by nvim built-in `vim.pack` (install/update/lockfile) and `lz.n` (lazy-loading).
For the vim.pack/lz.n interaction model and startup overhead implications, see
`~/.claude/skills/neovim/references/package-management.md`.

- **`lua/pack_specs.lua`** — `vim.pack.add()` registry, all remote plugins
- **`lua/pack_hooks.lua`** — `PackChanged` build hooks (mason, treesitter, etc.)
- **`lua/plugins/*.lua`** — lz.n specs (27 files), auto-discovered via `require("lz.n").load("plugins")`
- **Dev plugins** (modules/) — added to rtp manually in init.lua, use `load = function() end` in lz.n specs

**Disabling plugins:** Use `enabled = false` in the lz.n spec. Keep the full spec (setup config, keymaps, etc.) intact so the plugin can be re-enabled later by removing the flag. Never delete the spec body when disabling.

### lz.n gotcha: `trigger_load` requires handler registration

`trigger_load("plugin")` only finds plugins registered with a handler (event/cmd/ft/keys/colorscheme). Plugins with `lazy = true` but **no trigger field** are invisible — `trigger_load` silently skips them.

**Fixes for dep-only plugins:**
- Add a `cmd` trigger if the plugin has commands (e.g., mason.nvim → `cmd = {"Mason", ...}`)
- Remove `lazy = true` if the plugin can load eagerly (small cost)
- Inline `vim.cmd.packadd("plugin")` + setup in the parent's `before`/`after`

### blink.compat and nvim-cmp coexistence

blink.compat with `impersonate_nvim_cmp = true` replaces `package.loaded["cmp"]` with a mock. `require("cmp")` returns the mock even after clearing package.loaded (blink.compat has its own `lua/cmp/init.lua`). In `cmp.lua`, check `package.loaded["blink.cmp"]` at runtime and skip cmp setup when blink is active.

## Keymaps

Neovim can read the cmd key directly via `<D-...>` notation (e.g., `<D-v>`, `<D-\>`), as long as the terminal (kitty) isn't capturing it first. No `send_text` workaround needed.

## Zsh tree-sitter parser

A dedicated Zsh parser (`georgeharker/tree-sitter-zsh`) is now tier 1 in nvim-treesitter. It's better than the bash parser, which we don't use.

**Limitations:** No `indents.scm` (no auto-indent). `~` inside `[[ ]]` misparsed as bitwise negation ([#16](https://github.com/georgeharker/tree-sitter-zsh/issues/16)). Glob qualifier delimiters for `s::`, `n::`, `b::` must use `:`. Zsh builtins (`zparseopts`, `zstyle`, `autoload`) parse as generic commands. Still young, edge cases in complex `${(flags)name}` expansions may remain.

## Terminal cursor shape
`guicursor` with `t:ver25` doesn't affect cursor shape when actively typing in terminal buffers. The terminal process controls the cursor via escape sequences. This is a neovim limitation.

## Experimental cmdline (ui2)
Enabled in options.lua via `require('vim._core.ui2').enable()`. Eliminates "Press ENTER" prompts, highlights cmdline as you type, provides pager as a buffer. Still experimental — if issues arise, comment out the line in options.lua.

### Zsh filetype setup

Zsh files use compound filetype `"sh.zsh"` for vim regex syntax (loads sh patterns first, then zsh overrides). `vim.treesitter.language.register("zsh", "sh.zsh")` in `lua/autocmds/treesitter.lua` overrides the default first-component behaviour so treesitter uses the dedicated zsh parser.

## Fzf-lua and Ripgrep Colours

fzf-lua's grep runs ripgrep in a headless neovim subprocess (`nvim -u NONE -l spawn.lua`). This subprocess inherits `RIPGREP_CONFIG_PATH` but in practice the config file's `--colors` flags don't reliably apply — ripgrep falls back to default magenta paths. Fix: pass `--colors` explicitly in `rg_opts` (`lua/plugins/fuzzy.lua`), duplicating the relevant flags from `~/.config/ripgrep/config`.

Highlight groups for fzf-lua UI elements (prompt, header, keybinds, line numbers) are defined in the colorscheme generator (`~/dotfiles/colors/generators/colorscheme.lua`). fzf-lua sets them with `default=true`, so the colorscheme definitions take priority.

## LSP and Completion

General LSP knowledge (in-process servers, mason overrides, external config) is
in `~/.claude/skills/neovim/references/lsp.md`. Completion framework behaviour
(blink.cmp caching, `filterText`, community sources) is in
`~/.claude/skills/neovim/references/completion.md`.

Working examples: `modules/kitty-conf.nvim` (hover + completion),
`modules/agentic.nvim/lua/agentic/completion/lsp_server.lua` (trigger-character
completion for `/` and `@`).

### `lsp_ext/` — External Sources and Stubs

Extra type information for basedpyright, shared between neovim and the Claude lint hook.

```
lsp_ext/
├── extraPaths/             # Source directories added to pyright's extraPaths
│   ├── kitty-source/       # Kitty source (git submodule) — provides kitty.* types
│   └── pymol_modules/      # Symlink → ../pymol-open-source/modules/
├── python_stubs/           # .pyi stub files (pyright stubPath)
├── pymol-open-source/      # Full pymol source repo (git submodule)
└── r_lsp_dots.R            # R languageserver monkey-patch (see R Language Server section)
```

**How it works:** `lsp/basedpyright.lua` globs `lsp_ext/extraPaths/*/` for import resolution paths. The Claude `lib/lint-tier.sh` library uses the same glob when generating a fallback pyright config for projects without their own `pyrightconfig.json`.

**Adding a new source:** Drop the directory in `extraPaths/` (or symlink it there). Both neovim and the lint hook pick it up automatically — no config changes needed.

## R Language Server

The R languageserver doesn't resolve `...` forwarding — functions like `scale_y_log10(...)` that delegate to `scale_y_continuous(...)` only show `...` as a parameter, with no completion for the actual arguments.

**Patch:** `lsp_ext/r_lsp_dots.R` monkey-patches the languageserver at startup via `setHook(packageEvent("languageserver", "onLoad"), ...)`. Two patches:
1. **`get_formals`** (R6 `$set` on `PackageNamespace`) — when a function has only `...` as its formal, traces the body to find the target function and returns its formals instead.
2. **`arg_completion`** (`assignInNamespace`) — wraps the original to add `textEdit` (for reliable `" = "` insertion) and fix `data$funct` (so `completionItem/resolve` looks up docs for the underlying function, not the wrapper).

**Config:** `lsp/r_language_server.lua` sets a custom `cmd` that sources the patch before `languageserver::run()`.

**Mason-lspconfig override:** Mason's `automatic_enable` calls `vim.lsp.config()` which overrides fields from `lsp/*.lua` — not just `cmd` but also `filetypes` and other fields. Fix: in `lua/plugins/lsp.lua`, re-apply custom config after `setup()`. Already done for: `r_language_server` (custom cmd), `basedpyright` and `ruff` (compound filetypes).

## Miller DSL Highlighting

Custom tree-sitter grammar `miller` at `modules/tree-sitter-miller/` provides syntax highlighting for Miller's DSL (the language inside `put`/`filter`/`tee` verbs). Works in `*.mlr` files (nvim filetype `miller`) and is injected into zsh single-quoted strings after `put`/`filter`/`tee` verbs via `queries/zsh/injections.scm`.

Grammar name is `miller` to match nvim's built-in filetype — no `vim.treesitter.language.register()` needed. Registered in `lua/plugins/treesitter.lua` alongside pymol_select.

Queries use canonical flat structure: `modules/tree-sitter-miller/queries/highlights.scm` (not nested in a `miller/` subdir). nvim-treesitter symlinks `site/queries/miller -> modules/tree-sitter-miller/queries/` via `install_info`. Zsh injection query lives at `queries/zsh/injections.scm` (extends base zsh injections) — it's a query for the zsh parser, not the miller grammar, so it stays in the nvim config.

Regenerate after grammar changes: `cd modules/tree-sitter-miller && tree-sitter generate && cc -shared -o ~/.local/share/nvim/site/parser/miller.so -I src src/parser.c -O2`. Run `tree-sitter test` to validate.

## Miller DSL Completion

Inside `put`/`filter` DSL strings, the `blink_mlr` provider offers context-aware completions:
- After `$` → column names from referenced input files (existing)
- Otherwise → DSL builtin functions (223), keywords (40), and special variables (14)

Data generated from `mlr -F` and `mlr -K` by `lua/completion/mlr/miller_functions.json.sh` → `miller_functions.json`. Regenerate after Miller upgrades.

## PyMOL Selection Highlighting

Selection keywords (`name`, `chain`, `byres`, etc.) and representation names (`cartoon`, `sticks`, `surface`, etc.) inside Python strings are highlighted via a custom tree-sitter grammar `pymol_select` at `modules/tree-sitter-pymol-select/`. Injected into Python strings dynamically when pymol imports are detected (`ftplugin/python.lua`), scoped to function args and assignments (not docstrings).

Regenerate after grammar changes: `cd modules/tree-sitter-pymol-select && tree-sitter generate && cc -shared -o ~/.local/share/nvim/site/parser/pymol_select.so -I src src/parser.c -O2`. Restart neovim after recompiling.

**Known limitation:** Multi-part values with `+` (e.g. `chain A+B+C`) work via the `multi_value` token, but single-letter selector keywords (`b`, `q`, `x`, `y`, `z`) may be parsed as selectors instead of chain IDs in edge cases.

## PyMOL Completion

Three independent data sources, not derivations of each other:

1. **`lsp_ext/pymol-open-source/modules/pymol/`** — vendored pymol source. Added to `$PYTHONPATH` in `plugin/paths.lua` so basedpyright resolves `from pymol import cmd` (type info, hover, go-to-definition).
2. **`lua/completion/pymol/pymol.html`** — command reference generated via `cmd.write_html_ref()`. Covers `cmd.*` usage/args/examples. Not currently wired into completions.
3. **`lua/completion/pymol/pymol_settings_descriptions/*.md`** — scraped from PyMOL Wiki (separate content from the above). Used by the custom blink.cmp `pymol_settings` provider (`lua/completion/pymol/blink_pymol_settings.lua`), which offers setting completions inside `set()`/`cmd.set()` calls.

4. **`lua/completion/pymol/blink_pymol_select.lua`** — blink.cmp provider for selection keywords, builtins, operators, and representation names inside Python strings. Activates only where the `pymol_select` treesitter injection is active (uses `parser:language_for_range():lang()` to check), so no completions in docstrings or non-pymol strings.

Loading is conditional: `ftplugin/python.lua` scans the first 10 lines for `import.*pymol`, or use `<localleader>+` manually. Sets `vim.g.loaded_pymol` which gates the blink providers.

## File Templates

A new files get initial text from templates in `lua/autocmds/templates.lua` based on filetype.

## WGSL / Bevy Shader Highlighting

See [notes/wgsl.md](notes/wgsl.md).

## Testing

Read `tests/README.md`.

