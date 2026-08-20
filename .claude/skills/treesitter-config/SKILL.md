---
name: treesitter-config
description: This config's tree-sitter setup — the custom miller, pymol_select and smarts grammars, the dedicated zsh parser and sh.zsh filetype wiring, and their highlight/injection queries. Use when editing modules/tree-sitter-*, plugin/treesitter.lua, lua/plugins/treesitter.lua, lua/chem/*, queries under those grammars, ftplugin/python.lua pymol injection, or zsh filetype/parser behaviour. For grammar authoring (grammar.js, lexer conflicts) see the global treesitter skill; for nvim treesitter internals see the neovim skill's treesitter-integration.md.
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
"sh.zsh")` in `plugin/treesitter.lua` overrides the default
first-component behaviour so treesitter uses the dedicated zsh parser.

**Parser limitations:** No `indents.scm` (no auto-indent). `~` inside `[[ ]]`
misparsed as bitwise negation
([#16](https://github.com/georgeharker/tree-sitter-zsh/issues/16)). Glob
qualifier delimiters for `s::`, `n::`, `b::` must use `:`. Zsh builtins
(`zparseopts`, `zstyle`, `autoload`) parse as generic commands. Still young;
edge cases in complex `${(flags)name}` expansions may remain.

## Local grammars

Three grammars live in this repo, all wired the same way: registered in
`lua/plugins/treesitter.lua` through `install_info` with a local `path` and no
`queries` field. `init.lua` puts every `modules/*` directory on the runtimepath
at startup (via `lua/dev_plugins.lua`, shared with `tests/minimal_init.lua`), so
each grammar's queries reach nvim by the standard runtime layout —
`modules/tree-sitter-<name>/queries/<lang>/highlights.scm`.

Setting `install_info.queries` instead makes nvim-treesitter symlink
`~/.local/share/nvim/site/queries/<lang>` to that path (`install.lua`,
`do_link_queries`) — the field names the directory holding the `.scm` files, not
a runtime-shaped `queries/` root. Avoided here: the symlink goes stale whenever
the grammar moves, and it duplicates what the runtimepath already provides.

A grammar named after an existing nvim filetype needs nothing further. Otherwise
`vim.treesitter.language.register()` in `plugin/treesitter.lua` maps filetype to
parser, and `plugin/ftdetect.lua` maps extension to filetype.

Regenerate and rebuild after any grammar change, then restart neovim:

```zsh
cd modules/tree-sitter-<name>
tree-sitter generate
cc -shared -o ~/.local/share/nvim/site/parser/<name>.so -I src src/parser.c -O2
tree-sitter test
```

### Miller DSL

Grammar `miller` highlights Miller's DSL — the language inside `put`/`filter`/
`tee` verbs — in `*.mlr` files. The name matches nvim's built-in `miller`
filetype, so no `register()` call is needed.

Also injected into zsh single-quoted strings following those verbs, by
`queries/zsh/injections.scm`. That query drives the *zsh* parser, so it lives in
the nvim config rather than the miller module.

### PyMOL selection

Grammar `pymol_select` highlights selection keywords (`name`, `chain`, `byres`)
and representation names (`cartoon`, `sticks`, `surface`) inside Python strings.
Injected dynamically when `ftplugin/python.lua` detects a pymol import, scoped to
function arguments and assignments, not docstrings.

**Known limitation:** multi-part values with `+` (e.g. `chain A+B+C`) work via the
`multi_value` token, but single-letter selector keywords (`b`, `q`, `x`, `y`, `z`)
may parse as selectors instead of chain IDs.

### SMARTS / SMILES

Grammar `smarts` serves two filetypes. `smarts` (`.smarts`, `.sma`, `.smirks`)
resolves by name; `smiles` (`.smi`, `.smiles`, `.cxsmiles`, `.rsmi`) needs
`language.register("smarts", "smiles")`, which also carries the parser into
markdown fences tagged `smiles`.

The accepted alphabets are themselves generated: `symbols.js.py` probes RDKit to
write `symbols.js`, read by both `grammar.js` and the colour generator below.
Only the organic subset (`B C N O F P S Cl Br I`) parses unbracketed, so a bare
`Mn` is an ERROR node rather than manganese — and error recovery emits it as a
*sibling* of the record, leaving the valid prefix in a record that reports no
error.

Two highlight layers, concatenated because the second opens `; extends`
(`:h treesitter-query` — a base query plus every extending file on the
runtimepath): the grammar's own
`modules/tree-sitter-smarts/queries/smarts/highlights.scm`, and
`after/queries/smarts/highlights.scm` for bonds, the reaction arrow, the
recursion anchor and hydrogen (which has no symbol node of its own).

**Element identity is not a capture.** The highlighter derives a group from the
capture name and language alone, so N colours would need N patterns — and every
bracket atom is `element_symbol` or `atomic_number`, so it would evaluate all ~97
of them (102 ms per redrawn window, measured). Instead `queries/smarts/atoms.scm`
— a base query, no `; extends`, deliberately out of the highlights chain — finds
the atom nodes, and `lua/chem/highlight.lua` looks their text up in a
`text → group` hash (`Fe`, `fe` and `#26` all iron) and stores the result as
extmarks at priority 101, above the grammar's `@constant`. Text-driven, so it is
computed on `on_changedtree` rather than per redraw; see
`notes/PLAN-smarts-element-perf.md` for why not a decoration provider.

`lua/chem/elements.lua.py` (uv script, `ase` + `coloraide`) writes
`lua/chem/elements.lua`: a Jmol CPK colour per element, lightened or darkened
until it clears 2:1 WCAG contrast against the light and dark colorschemes'
`Normal` background, plus which symbols have a legal aromatic spelling and the
suspect symbols with no isolable compound. Rerun it after `symbols.js` changes;
`tests/plenary/smarts_spec.lua` fails when the two files disagree.

`lua/chem/highlight.lua` both defines the `@chem.*` groups and paints the marks.
`plugin/treesitter.lua` calls its `setup()` — a plugin rather than an ftplugin, so
fenced structures with no filetype of their own still find the groups defined —
and its `attach()` beside `vim.treesitter.start`, gated on `vim.b.ts_highlight` so
the `largefile` and disabled-filetype policies are inherited rather than
restated.
