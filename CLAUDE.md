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

### vim.NIL in buffer variables
When storing sparse Lua tables (e.g., `{[3] = 10}`) in `vim.b`, Vim pads missing indices with `vim.NIL`. This userdata value:
- Is NOT equal to Lua `nil` (`vim.NIL == nil` is `false`)
- Is truthy in conditions
- **Cannot be indexed** — causes "attempt to index a userdata value"

Pattern `t[k] and t[k][j]` fails when `t[k]` is `vim.NIL`. Use a helper to convert `vim.NIL` back to `nil` when reading from `vim.b`. See `denilify()` in `ftplugin/tsv.lua`.

### TSV column hiding and undo
The TSV ftplugin's column hiding (`zc`/`zo`/`za`) works by actually removing text from the buffer and storing it in `vim.b`. The `modified` flag is preserved so hiding doesn't mark clean buffers as dirty.

**Undo behaviour**: Hide operations are in the undo tree. When undo restores hidden text, a `TextChanged` autocmd clears the stale hidden state to keep `vim.b.tsv_hidden` in sync. Press `zc` to re-hide after undo if needed.

**Why not use concealment?** Concealment (`conceal` extmark option) was attempted but has fundamental issues:
- Cursor still navigates through concealed text (confusing)
- Tab alignment breaks because tabs expand based on buffer position, not visual position
- Would require reimplementing entire tab/column system with virtual text

## Treesitter Query `#match?` Escaping

Tree-sitter query strings strip backslashes from unknown escape sequences (`\~` → `~`, `\.` → `.`). Nvim's `check_magic` then prepends `\v` (very magic) to `#match?` patterns. In `\v` mode, `~` means "last substitute string" (E33 crash if no prior `:s`), and `.` means "any char".

**Double-escape** Vim-special characters in `#match?` patterns: `\\~` → `\~` (literal tilde), `\\.` → `\.` (literal dot). Bare `?`, `+`, `*`, `{n,m}` work as quantifiers in `\v` mode without escaping.

`#lua-match?` uses Lua patterns (no Vim regex) and avoids this issue entirely.

## Compound Filetypes and Runtime Loading

Neovim splits compound filetypes (e.g., `"sh.zsh"`) on `.` and iterates left-to-right. For each component, `runtime!` searches **all** rtp entries (including `after/` dirs). This means `after/syntax/a.vim` loads *before* `syntax/b.vim`.

Load order for filetype `"a.b"`:
1. `syntax/a.{vim,lua}` across all rtp (user → bundled → after/)
2. `syntax/b.{vim,lua}` across all rtp (user → bundled → after/)

Same pattern for `ftplugin/` and `indent/`. Last-write-wins for options, so `b` overrides `a`.

**Treesitter uses only the first component** for parser and query selection (`vim.split(filetype, '.')[1]`). Override with `vim.treesitter.language.register("lang", "a.b")` to use a different parser for the compound filetype. Queries load from `queries/{parser_lang}/`, not from the second filetype component.

**`; inherits`** in treesitter queries loads queries from another language's directory (e.g., `; inherits sh` loads `queries/sh/highlights.scm`). Use this to share queries between parsers.

**`b:current_syntax` guard**: The first syntax component sets it. Subsequent components with the standard `if exists("b:current_syntax") | finish` guard will be skipped — only addon-style syntax files (without the guard) load for the second component.

### Zsh filetype setup

Zsh files use compound filetype `"sh.zsh"` for vim regex syntax (loads sh patterns first, then zsh overrides). `vim.treesitter.language.register("zsh", "sh.zsh")` in `lua/autocmds/treesitter.lua` overrides the default first-component behaviour so treesitter uses the dedicated zsh parser.

## Large File Handling

Files >50MB are handled specially to avoid freezing. See `lua/largefile.lua` and `lua/autocmds/largefile.lua`.

### Key learnings

**Command line files load before plugin/ autocmds register.** Autocmds in `plugin/` can't intercept files specified on the command line — they're processed before those files load. Solution: check `vim.fn.argv()` in `init.lua` before lazy.nvim.

**`eventignore` alone doesn't help.** Even `eventignore=all` doesn't prevent freezes if plugins do expensive work during their setup (not via autocmds).

**Must remove file from argv AND wipe buffer.** Use `argdelete *` + `bwipeout!` to prevent neovim from auto-opening the file. Just `argdelete` isn't enough — the buffer may already exist.

**`noautocmd edit` still triggers expensive operations.** The freeze wasn't from autocmds but from something during the edit command itself. Using regular `edit` after blocking events works better.

**Double `vim.schedule` after UIEnter.** Single schedule runs before startup plugins (dashboard, etc.) finish. Double schedule ensures file opens after everything else completes.

**`nvim --clean` works because no plugins load.** Useful baseline for debugging — if `--clean` works but normal startup doesn't, the issue is plugin-related.

### Event sequence for command line files

1. `init.lua` runs — detect large file, clear argv, wipe buffer, set eventignore
2. Plugins load (events blocked)
3. `UIEnter` fires
4. Double-scheduled callback opens file with `edit`
5. `BufAdd` in largefile.lua intercepts (for session files, not needed for command line)

## Macro Replay and RecordingLeave

`RecordingLeave` fires before the register is populated (in headless tests at least). With `vim.schedule`, the register is available in interactive use.

**`norm @b` doesn't work in `vim.schedule` after `RecordingLeave`.** The `@` command silently fails — no error, but no text changes. Use `vim.fn.feedkeys(reg, "nx")` with the raw register content instead. The `"n"` flag prevents remapping, `"x"` executes immediately.

**Don't pass register content through `nvim_replace_termcodes`.** Register content is already in vim's internal key encoding (`<80><fd>...` sequences for modifiers). `replace_termcodes` double-encodes it, producing garbage.

See `lua/keymaps/blockim.lua` for working example.

## Dynamic Esc Mappings

Avoid `vim.keymap.set`/`vim.keymap.del` for temporary `<Esc>` overrides — they conflict with any global Esc mapping and with other plugins that may also remap Esc. Prefer `RecordingLeave` or `ModeChanged` autocmds, or use a flag variable checked by a single global Esc mapping.

## Debugging Tips

**`nvim --clean --headless`** is useful for automated tests but has limitations — macro recording/register timing differs from interactive use. Some things that work interactively fail in headless mode and vice versa.

**`pcall` inside `vim.schedule`** — errors in scheduled callbacks can be silent. Wrap in `pcall` and print errors to surface them.

**`vim.fn.getreg("x")` returns `""` not `nil` for empty registers.** The pattern `not vim.fn.getreg("x")` is always false (empty string is truthy in Lua). Use `vim.fn.getreg("x") == ""` instead.

## In-process LSP Servers

Neovim can host LSP servers inside its own process — no external binary needed. See `kitty-conf.nvim` for a working example.

- **`lsp/*.lua`**: Return `{ cmd = function() ... end, filetypes = {...} }`. The `cmd` function returns `request`, `notify`, `is_closing`, `terminate` handlers.
- **`plugin/*.lua`**: Call `vim.lsp.enable()` here, not in `ftplugin/`. `enable` registers a `FileType` autocmd internally, so calling it from ftplugin is circular and too late. The `filetypes` in the `lsp/` config gates which buffers actually start the LSP.
- **blink.cmp community sources**: Convention is `lua/<plugin-name>/init.lua`, referenced as `module = "<plugin-name>"` in the provider config.
- **blink.cmp CompletionItem fields**: `labelDetails.description` shows inline in the menu (right of label). `labelDetails.detail` appends directly after the label (no gap). `documentation` shows in the hover popup when an item is selected. Top-level `detail` only shows in the documentation popup, not inline in the menu.

## Testing

Unit tests use plenary.nvim. See `tests/README.md`.
