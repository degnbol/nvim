---
name: treesitter-config
description: This config's tree-sitter setup — the custom miller and pymol_select grammars, the dedicated zsh parser and sh.zsh filetype wiring, and their highlight/injection queries. Use when editing modules/tree-sitter-*, lua/plugins/treesitter.lua, queries under those grammars, ftplugin/python.lua pymol injection, or zsh filetype/parser behaviour. For grammar authoring (grammar.js, lexer conflicts) see the global treesitter skill; for nvim treesitter internals see the neovim skill's treesitter-integration.md.
---

# Tree-sitter (this config)

Custom grammars and parser wiring specific to this nvim config. Grammar
authoring itself is in the global `treesitter` skill; neovim treesitter
integration (injections, `language_for_range`, compilation) is in
`~/.claude/skills/neovim/references/treesitter-integration.md`.

## Zsh parser and filetype

A dedicated Zsh parser (`georgeharker/tree-sitter-zsh`) is tier 1 in
nvim-treesitter — better than the bash parser, which this config doesn't use.

Zsh files use compound filetype `"sh.zsh"` for vim regex syntax (loads sh
patterns first, then zsh overrides). `vim.treesitter.language.register("zsh",
"sh.zsh")` in `lua/autocmds/treesitter.lua` overrides the default
first-component behaviour so treesitter uses the dedicated zsh parser.

**Parser limitations:** No `indents.scm` (no auto-indent). `~` inside `[[ ]]`
misparsed as bitwise negation
([#16](https://github.com/georgeharker/tree-sitter-zsh/issues/16)). Glob
qualifier delimiters for `s::`, `n::`, `b::` must use `:`. Zsh builtins
(`zparseopts`, `zstyle`, `autoload`) parse as generic commands. Still young;
edge cases in complex `${(flags)name}` expansions may remain.

## Miller DSL grammar

Custom grammar `miller` at `modules/tree-sitter-miller/` highlights Miller's DSL
(the language inside `put`/`filter`/`tee` verbs). Works in `*.mlr` files (nvim
filetype `miller`) and is injected into zsh single-quoted strings after
`put`/`filter`/`tee` verbs via `queries/zsh/injections.scm`.

Grammar name is `miller` to match nvim's built-in filetype — no
`vim.treesitter.language.register()` needed. Registered in
`lua/plugins/treesitter.lua` alongside `pymol_select`.

Queries use canonical flat structure:
`modules/tree-sitter-miller/queries/highlights.scm` (not nested in a `miller/`
subdir). nvim-treesitter symlinks `site/queries/miller ->
modules/tree-sitter-miller/queries/` via `install_info`. The zsh injection query
lives at `queries/zsh/injections.scm` (extends base zsh injections) — it's a
query for the zsh parser, not the miller grammar, so it stays in the nvim
config.

Regenerate after grammar changes: `cd modules/tree-sitter-miller && tree-sitter
generate && cc -shared -o ~/.local/share/nvim/site/parser/miller.so -I src
src/parser.c -O2`. Run `tree-sitter test` to validate.

## PyMOL selection grammar

Selection keywords (`name`, `chain`, `byres`, etc.) and representation names
(`cartoon`, `sticks`, `surface`, etc.) inside Python strings are highlighted via
grammar `pymol_select` at `modules/tree-sitter-pymol-select/`. Injected into
Python strings dynamically when pymol imports are detected
(`ftplugin/python.lua`), scoped to function args and assignments (not
docstrings).

Regenerate after grammar changes: `cd modules/tree-sitter-pymol-select &&
tree-sitter generate && cc -shared -o
~/.local/share/nvim/site/parser/pymol_select.so -I src src/parser.c -O2`.
Restart neovim after recompiling.

**Known limitation:** Multi-part values with `+` (e.g. `chain A+B+C`) work via
the `multi_value` token, but single-letter selector keywords (`b`, `q`, `x`,
`y`, `z`) may be parsed as selectors instead of chain IDs in edge cases.
