# Add-Snippet Workflow (`lua/luasnippets/add.lua`)

Creates new LuaSnip snippet definitions from visual selections.

## Functions

- `add_from_visual()` — visual-select text, then call this. It:
  1. Reads the selected lines
  2. Determines filetype and opens `luasnippets/<ft>.lua`
  3. Inserts a blank line before the closing `}` of the return table
  4. Expands a LuaSnip meta-snippet at that position with:
     - `i(1)` — trigger name
     - `i(2)` — description
     - `c(3, ...)` — choice between `t""` body and `fmta()` body

- `edit()` — opens `luasnippets/<ft>.lua` for the current filetype.

## Implementation details

- Uses `ls.snip_expand()` (not `vim.snippet.expand`) to expand the meta-snippet.
- The snippet object is built dynamically via `ls.s("", {...})` with an empty trigger.
- Text is escaped for Lua double-quoted strings (`escape_lua`). No LSP snippet escaping needed.
- For multi-line fmta body, uses `[[...]]` long strings (no Lua escaping needed inside).
- File creation: new files get the standard header (`---@diagnostic disable: undefined-global\nreturn {\n}\n`).
- **Blink keymap workaround:** Before expanding, fires `nvim_exec_autocmds('InsertEnter')` to force blink.cmp to apply its buffer-local keymaps (snippet jump keys) to the newly opened buffer. Without this, `<C-.>`/`<C-,>` don't work until the user manually enters insert mode, because blink only applies keymaps on `InsertEnter`.

## Conditional snippets

`lua/luasnippets/lazy.lua` and `python_pymol.lua` use `ls.add_snippets()` inside a
conditional guard (e.g. checking buffer name). These load programmatically rather than
via the `return {}` pattern.
