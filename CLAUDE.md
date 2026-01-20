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
