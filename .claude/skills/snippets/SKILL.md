---
name: snippets
description: >
  Neovim snippet development with LuaSnip and vim.snippet. Use when editing
  luasnippets/*.lua files, lua/luasnippets/ helpers, writing new snippets,
  or discussing snippet nodes (insert, text, choice, function, dynamic, fmt/fmta).
---

# Snippets in Neovim

Two snippet systems are active:

1. **`vim.snippet`** — built-in LSP snippet expansion (`$1`, `${1:placeholder}`, `$0`). Used for file templates. Documented in the neovim skill.
2. **LuaSnip** — node-based snippet engine with choice nodes, function nodes, dynamic nodes, regex triggers. Handles all interactive snippets.

## Directory layout

```
luasnippets/<ft>.lua              # Per-filetype snippet definitions (LuaSnip auto-loads)
lua/luasnippets/add.lua           # Create snippets from visual selection
lua/luasnippets/lazy.lua          # Conditional: lazy.nvim plugin spec snippets
lua/luasnippets/python_pymol.lua  # Conditional: pymol snippets
lua/utils/luasnip.lua             # Re-exports + custom utilities
```

## Snippet file boilerplate

```lua
---@diagnostic disable: unused-local
local ls = require "utils/luasnip"
local s, t, i, c, f, d, sn, fmta, conds, rep, ms = ls.s, ls.t, ls.i, ls.c, ls.f, ls.d, ls.sn, ls.fmta, ls.conds, ls.rep, ls.ms
return {
-- snippets here
}
```

## Completion integration

LuaSnip snippets appear in blink.cmp via the built-in `luasnip` source. Snippet jump keymaps (`<C-.>`/`<C-,>`) are shared between blink and LuaSnip. For blink internals (keyword detection, trigger flow, fuzzy matching), see the [blink skill](.claude/skills/blink/SKILL.md).

## References

- [LuaSnip node types, choice nodes, triggers](references/luasnip-nodes.md)
- [Custom utilities in `lua/utils/luasnip.lua`](references/utils-luasnip.md)
- [Add-snippet workflow (`lua/luasnippets/add.lua`)](references/add-snippet.md)
