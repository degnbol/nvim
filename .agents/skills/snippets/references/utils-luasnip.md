# Custom LuaSnip Utilities (`lua/utils/luasnip.lua`)

Re-exports all standard LuaSnip symbols so snippet files can destructure from one require:

```lua
local ls = require "utils/luasnip"
local s, t, i, c, f, d, sn, fmta, conds, rep, ms = ls.s, ls.t, ls.i, ls.c, ls.f, ls.d, ls.sn, ls.fmta, ls.conds, ls.rep, ls.ms
```

## Custom functions

- `get_visual(args, parent)` — returns `sn` with visual selection text (`SELECT_RAW`) or empty insert node. Use with `d()` node.
- `re(i)` — function node returning regex capture group `i` from the trigger match.
- `virt(text)` — node_ext_opts table adding virtual text hint (e.g. for choice nodes: `c(1, {t("", virt("^l for alt")), t"alt"})`).
- `match_ahead(n)` — custom trigEngine that requires the trigger pattern to match `n` characters ahead of the cursor position. The ahead characters are not replaced.
- `upper(i)` — function node that uppercases the text of node `i`.
- `rep_jump(N, ref)` — like `rep(ref)` but creates a dynamic node at position `N` that is jumpable (editable copy of node `ref`).
- `upper_jump(N, ref)` — jumpable uppercase copy of node `ref`.
