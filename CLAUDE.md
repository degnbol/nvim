# Neovim Config

`~/nvim` = `~/.config/nvim` = `~/dotfiles/config/nvim` (symlinks). All the same directory.

## Keymap Notes

### Cmd key (`<D-...>`) in Neovim
Neovim can read the cmd key directly via `<D-...>` notation (e.g., `<D-v>`, `<D-\>`), as long as the terminal (kitty) isn't capturing it first. No `send_text` workaround needed.

### Key patterns
- `<CR>` / `<S-CR>` — REPL run/paste (with motion or visual)
- `<M-CR>` — Agent send (with motion or visual)
- `<D-]>` / `<D-[>` — Switch to REPL window (kitty level)
- `<D-\>` — Toggle agent window

### AI/Agent keymaps (`<leader>i` = intelligence)
Shared keymaps work with whichever plugin is enabled (claudecode.nvim or agentic.nvim).
Toggle `enabled` in `lua/plugins/agents.lua` to switch.

| Key | Action |
|-----|--------|
| `<leader>ii` | Toggle agent |
| `<leader>if` | Focus agent |
| `<leader>iq` | Close agent |
| `<leader>in` | New session |
| `<leader>ix` | Stop generation |
| `<leader>ib` | Add current buffer |
| `<leader>is` | Send selection (visual) |

**claudecode.nvim only:**
| `<leader>ir` | Resume session |
| `<leader>ic` | Continue session |
| `<leader>im` | Select model |
| `<leader>ia` | Accept diff |
| `<leader>id` | Deny diff |

## Treesitter Commands

nvim-treesitter removed `TSBufToggle`, `TSToggle`, `TSInstallInfo`, and other `TS*` commands. Treesitter highlighting is now built into neovim core — use `vim.treesitter.start()`/`vim.treesitter.stop()` to control it. For install info, use `:checkhealth nvim-treesitter`.

### Zsh tree-sitter parser (2025-02)

A dedicated Zsh parser (`georgeharker/tree-sitter-zsh`) is now tier 1 in nvim-treesitter. Install with `:TSInstall zsh` — no more `vim.treesitter.language.register("bash", "zsh")` hack. Actively maintained, handles most Zsh-specific syntax (parameter expansion flags, glob qualifiers, short-form loops, anonymous functions).

**Limitations:** No `indents.scm` (no auto-indent). `~` inside `[[ ]]` misparsed as bitwise negation ([#16](https://github.com/georgeharker/tree-sitter-zsh/issues/16)). Glob qualifier delimiters for `s::`, `n::`, `b::` must use `:`. Zsh builtins (`zparseopts`, `zstyle`, `autoload`) parse as generic commands. Still young (~5 months), edge cases in complex `${(flags)name}` expansions may remain.

## Known Issues

### Terminal cursor shape
`guicursor` with `t:ver25` doesn't affect cursor shape when actively typing in terminal buffers. The terminal process controls the cursor via escape sequences. This is a neovim limitation.

### Experimental cmdline (`vim._extui`)
Disabled in options.lua due to treesitter query error with "tab" node after nvim update. Re-enable when fixed upstream.

### Markdown treesitter crash (nvim 0.12)
The bundled markdown parser in nvim 0.12-dev crashes when `vim.treesitter.start()` is called during initial buffer load with `foldmethod=expr` and treesitter foldexpr. Workaround in `lua/autocmds/treesitter.lua` uses `vim.schedule()` to delay treesitter start for markdown files. Remove workaround when fixed upstream.

### TSV column hiding and undo
The TSV ftplugin's column hiding (`zc`/`zo`/`za`) works by actually removing text from the buffer and storing it in `vim.b`. The `modified` flag is preserved so hiding doesn't mark clean buffers as dirty.

**Undo behaviour**: Hide operations are in the undo tree. When undo restores hidden text, a `TextChanged` autocmd clears the stale hidden state to keep `vim.b.tsv_hidden` in sync. Press `zc` to re-hide after undo if needed.

**Why not use concealment?** Concealment (`conceal` extmark option) was attempted but has fundamental issues:
- Cursor still navigates through concealed text (confusing)
- Tab alignment breaks because tabs expand based on buffer position, not visual position
- Would require reimplementing entire tab/column system with virtual text

### Zsh filetype setup

Zsh files use compound filetype `"sh.zsh"` for vim regex syntax (loads sh patterns first, then zsh overrides). `vim.treesitter.language.register("zsh", "sh.zsh")` in `lua/autocmds/treesitter.lua` overrides the default first-component behaviour so treesitter uses the dedicated zsh parser.

## Large File Handling

Files >50MB are handled specially to avoid freezing. See `lua/largefile.lua` and `lua/autocmds/largefile.lua`.

### Event sequence for command line files

1. `init.lua` runs — detect large file, clear argv, wipe buffer, set eventignore
2. Plugins load (events blocked)
3. `UIEnter` fires
4. Double-scheduled callback opens file with `edit`
5. `BufAdd` in largefile.lua intercepts (for session files, not needed for command line)

## In-process LSP Servers

Neovim can host LSP servers inside its own process — no external binary needed. See `kitty-conf.nvim` for a working example.

- **`lsp/*.lua`**: Return `{ cmd = function() ... end, filetypes = {...} }`. The `cmd` function returns `request`, `notify`, `is_closing`, `terminate` handlers.
- **blink.cmp community sources**: Convention is `lua/<plugin-name>/init.lua`, referenced as `module = "<plugin-name>"` in the provider config.
- **blink.cmp CompletionItem fields**: `labelDetails.description` shows inline in the menu (right of label). `labelDetails.detail` appends directly after the label (no gap). `documentation` shows in the hover popup when an item is selected. Top-level `detail` only shows in the documentation popup, not inline in the menu.

## R Language Server

The R languageserver doesn't resolve `...` forwarding — functions like `scale_y_log10(...)` that delegate to `scale_y_continuous(...)` only show `...` as a parameter, with no completion for the actual arguments.

**Patch:** `lsp_ext/r_lsp_dots.R` monkey-patches the languageserver at startup via `setHook(packageEvent("languageserver", "onLoad"), ...)`. Two patches:
1. **`get_formals`** (R6 `$set` on `PackageNamespace`) — when a function has only `...` as its formal, traces the body to find the target function and returns its formals instead.
2. **`arg_completion`** (`assignInNamespace`) — wraps the original to add `textEdit` (for reliable `" = "` insertion) and fix `data$funct` (so `completionItem/resolve` looks up docs for the underlying function, not the wrapper).

**Config:** `lsp/r_language_server.lua` sets a custom `cmd` that sources the patch before `languageserver::run()`.

**Mason-lspconfig override:** Mason's `automatic_enable` calls `vim.lsp.config()` with `cmd = { "r-languageserver" }` (Mason's wrapper), overriding the custom cmd from `lsp/*.lua`. Fix: in `lua/plugins/lsp.lua`, the mason-lspconfig config function calls `vim.lsp.config('r_language_server', { cmd = ... })` after `setup()` to re-apply our custom cmd.

## PyMOL Selection Highlighting

Selection keywords (`name`, `chain`, `byres`, etc.) and representation names (`cartoon`, `sticks`, `surface`, etc.) inside Python strings are highlighted via a custom tree-sitter grammar `pymol_select` at `tree-sitter-pymol-select/`. Injected into Python strings dynamically when pymol imports are detected (`ftplugin/python.lua`), scoped to function args and assignments (not docstrings).

Regenerate after grammar changes: `cd tree-sitter-pymol-select && tree-sitter generate && cc -shared -o ~/.local/share/nvim/site/parser/pymol_select.so -I src src/parser.c -O2`. Restart neovim after recompiling.

**Known limitation:** Multi-part values with `+` (e.g. `chain A+B+C`) work via the `multi_value` token, but single-letter selector keywords (`b`, `q`, `x`, `y`, `z`) may be parsed as selectors instead of chain IDs in edge cases.

## PyMOL Completion

Three independent data sources, not derivations of each other:

1. **`lsp_ext/pymol-open-source/modules/pymol/`** — vendored pymol source. Added to `$PYTHONPATH` in `plugin/paths.lua` so basedpyright resolves `from pymol import cmd` (type info, hover, go-to-definition).
2. **`lua/completion/pymol/pymol.html`** — command reference generated via `cmd.write_html_ref()`. Covers `cmd.*` usage/args/examples. Not currently wired into completions.
3. **`lua/completion/pymol/pymol_settings_descriptions/*.md`** — scraped from PyMOL Wiki (separate content from the above). Used by the custom blink.cmp `pymol_settings` provider (`lua/completion/pymol/blink_pymol_settings.lua`), which offers setting completions inside `set()`/`cmd.set()` calls.

4. **`lua/completion/pymol/blink_pymol_select.lua`** — blink.cmp provider for selection keywords, builtins, operators, and representation names inside Python strings. Activates only where the `pymol_select` treesitter injection is active (uses `parser:language_for_range():lang()` to check), so no completions in docstrings or non-pymol strings.

Loading is conditional: `ftplugin/python.lua` scans the first 10 lines for `import.*pymol`, or use `<localleader>+` manually. Sets `vim.g.loaded_pymol` which gates the blink providers.

## Testing

Unit tests use plenary.nvim. See `tests/README.md`.
