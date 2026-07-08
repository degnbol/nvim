# Neovim Config

`~/nvim` = `~/.config/nvim` = `~/dotfiles/config/nvim` (symlinks). All the same directory.

Config-specific detail lives in local skills (auto-loaded by description), not here:
`plugins`, `blink`, `snippets`, `lsp`, `treesitter-config`, `math-images`,
`version-audit`. General neovim knowledge is in the global `neovim` skill and its
`references/`. Add new detail to the relevant skill rather than growing this file.

## Shared utilities (`lua/utils/`)

Before reaching for `vim.api.nvim_replace_termcodes`, visual-range extraction, keymap wrappers, ColorScheme hooks, etc., check `lua/utils/init.lua`, `lua/utils/keymap.lua`, and `lua/utils/highlights.lua` — most common patterns are already wrapped (e.g. `nvim_code(s)` for `<Esc>`/`<CR>` termcode conversion, `last_visual_range()`, `get_cursor`/`set_cursor`, `press(keys)`, `hi.onColorScheme(fn)` to register a ColorScheme autocmd that also fires once on load). When you write a helper that's reusable across files, add it to `utils/` rather than duplicating inline.

## Experimental cmdline (ui2)

Enabled in `lua/options.lua` via `require('vim._core.ui2').enable()`. If cmdline behaviour breaks, comment out that line first. (What ui2 does: neovim skill → `references/ui.md`.)

## File Templates

New files get initial text from templates in `lua/autocmds/templates.lua` based on filetype.

## WGSL / Bevy Shader Highlighting

See [notes/wgsl.md](notes/wgsl.md).

## Testing

Read `tests/README.md`.
