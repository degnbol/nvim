# Neovim Config

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

## Testing

Unit tests use plenary.nvim. See `tests/README.md`.
