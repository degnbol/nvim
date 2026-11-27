# LuaSnip Node Reference

## Node types

| Node | Constructor | Purpose |
|------|-------------|---------|
| Text | `t"text"` or `t{"line1", "line2"}` | Static text |
| Insert | `i(N, "default")` | Tab stop with optional default |
| Choice | `c(N, {opt1, opt2})` | Cycle options with `<C-n>`/`<C-p>` |
| Function | `f(fn, argnode_refs)` | Computed text from other nodes |
| Dynamic | `d(N, fn, argnode_refs)` | Returns a snippet_node dynamically |
| Snippet node | `sn(N, nodes)` | Group of nodes as a single node |
| `fmta` | `fmta("template <>", {i(1)})` | Format with `<>` placeholders |
| `rep` | `rep(N)` | Mirror another node's text |

## Choice nodes

- Options can be any node type: `t`, `i`, `sn`, `fmta(...)`, etc.
- `t` nodes work directly as choice options — no `sn()` wrapper needed.
- Cycle with `<C-n>` / `<C-p>` (or however LuaSnip choice keys are mapped).

## `fmta` vs `fmt`

- `fmta` uses `<>` as placeholder delimiters (angle brackets).
- `fmt` uses `{}` as placeholder delimiters (curly braces).
- Escape delimiters by doubling: `<<>>` in fmta, `{{}}` in fmt.
- Both return a list of nodes, usable directly as snippet body or inside `sn()`.

## Multi-line text nodes

`t{"line1", "line2"}` produces two lines. Each table entry is one line. Empty string `""` produces a blank line.

## Snippet triggers

| Key | Values |
|-----|--------|
| `trig` | Trigger string or regex |
| `dscr` | Description shown in completion menu |
| `snippetType` | `"snippet"` (default) or `"autosnippet"` |
| `wordTrig` | `true` (default) — only trigger at word boundary |
| `trigEngine` | Custom function for advanced matching |
| `condition` | Function returning bool, checked before expansion |
| `show_condition` | Function for completion menu visibility |
| `priority` | Higher overrides lower (default 1000) |
