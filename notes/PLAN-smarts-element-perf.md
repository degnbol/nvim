# PLAN: element colours without a pattern per element

The generated query layer costs 102 ms per redrawn window on bracket-heavy input.
Highlighting is per-line, so the buffer visibly paints one line at a time. Element
identity has to stop being expressed as query patterns.

## The regression

Per-line `iter_captures` over a 50-line window, the shape the highlighter uses:

| buffer | with the generated layer | without it | lua, for scale |
| --- | --- | --- | --- |
| organic (`CC(=O)Oc1ccccc1C(=O)O`) | 6.05 ms | 0.31 ms | 0.62 ms |
| bracket-heavy (`[#6;H1][C;R1](=[O;X1])[N;!H0]`) | 102.14 ms | 1.35 ms | — |

The earlier cost analysis measured the wrong axis: it compared per-element
node-type lists against one uniform union and concluded the lists were
load-bearing. They only ever helped `organic_symbol`. Every bracket atom is
`element_symbol` or `atomic_number`, and those appear in **all 97** patterns, so a
bracket atom evaluates ~97 `#any-of?` predicates either way — which is why the
bracket column barely moved between the two variants (45.9 vs 49.9 ms) while the
organic one did.

## Why the query cannot do this

`highlighter.lua:38-49` — `get_hl_from_capture` builds the group from the capture
name and the language, nothing else. `on_range_impl` reads only `priority`,
`conceal`, `conceal_lines` and `url` out of metadata, and derives `spell` from the
capture name; `metadata.range`/`offset!` move a range but never the group. A
`#set!` directive computing `hl_group` from node text works in some query
consumers and is ignored here.

So N colours require N capture names, and N capture names require N patterns.
Tree-sitter indexes patterns by the start node's type, which bounds the damage per
*node type* — but the periodic table is one node type. The design is O(elements)
per atom by construction, not by mistuning.

## The fix

Find atoms with a query, colour them in Lua, and store the result as extmarks that
survive between redraws. Element identity depends only on buffer text, so it must
be computed when the text changes, not when the screen repaints.

`queries/smarts/atoms.scm` — a base query, no `; extends`, so it does not join the
highlights chain. On disk rather than built in Lua for the reason the previous plan
gave for the element query: `:InspectTree` and `ts_captures.sh` can be pointed at
it with no Lua having run.

```scheme
[(organic_symbol) (aromatic_symbol) (element_symbol) (atomic_number)] @atom
```

A `text → highlight group` hash built once from `chem/elements` — `Fe`, `fe` and
`#26` all mapping to `@chem.element.iron`, the twenty suspect symbols to
`@chem.suspect`. Lookup is exact, which is what keeps `[#0]`, `[#06]` and `[#999]`
out: all three parse, none is an element. That rationale currently lives in the
generated query's header and needs a home in the lookup's construction.

Measured with a decoration provider — an upper bound on the persistent design,
which pays this only on edit rather than on every redraw:

| buffer | now | finder + lookup |
| --- | --- | --- |
| organic | 6.13 ms | 0.98 ms |
| bracket-heavy | 101.98 ms | 2.17 ms |

## Wiring

`LanguageTree:register_cbs({ on_changedtree = … }, true)` — recursive, so injected
trees are covered — paired with `on_child_added`/`on_child_removed`, as
`highlighter.lua:114-137` does. Not `on_lines`: it fires *before* the reparse, so
the tree at that moment does not describe the new text. This is the one place the
design departs from mini.hipatterns, which matches regexes and needs no tree.

The changed `Range6`s bound the work. Clear the namespace over each range before
re-adding. **No debounce** — `on_changedtree` already fires only after a parse, and
0.12's parses are async and range-limited, so a timer on top adds lag to an
existing deferral.

Attach on the same signal as `vim.treesitter.start`, gated on
`vim.b[buf].ts_highlight`. That inherits two policies rather than restating them:
`plugin/treesitter.lua:63` returns before `start` for a `largefile` buffer, and
line 48 disables treesitter for `latex`/`plaintex`/`tex`. In both, element colours
correctly do not appear. Do **not** gate with `vim.treesitter.get_parser` — it
creates and caches a parser for any buffer whose filetype maps to an installed
language (verified: on a fresh lua buffer it returns a parser while `ts_highlight`
stays `nil`), so using it as a question instantiates parsers everywhere.

**Priority 101**, and the constraint is beating the grammar's own query, which
gives every element node `@constant` at the default 100 — not the `after/` layer.
Disjoint attributes merge across namespaces regardless of priority, so bonds
(`bold`, no colour) and anchors (`underline`, no colour) still let the element
colour through. Conflicting attributes at *equal* priority resolve by set order,
which is why the element mark must sit above 100 rather than at it.

## What each file becomes

**`queries/smarts/highlights.scm`** — deleted. nvim goes back to concatenating the
grammar's query and the `after/` one.

**`lua/chem/elements.lua.py`** — stops emitting a query; `spellings()` and
`pattern()` go. It still reads `symbols.js`, to know which symbols exist and which
have a legal aromatic spelling, so only those get a lowercase key. Hydrogen keeps
its exclusion from the symbol keys — the grammar cannot produce `H` as an
`element_symbol` — but the reason should survive as a comment now that it is no
longer implicit in a generated pattern.

**`lua/chem/elements.lua`** — gains the suspect symbols as a **separate export**,
which today exist only inside the query being deleted. Every `ChemElement` keeps
mandatory `light` and `dark`, so "each element has a colour clearing its
background's floor" stays a total invariant instead of one each loop must
remember.

**`lua/chem/highlight.lua`** — `lua/chem/palette.lua` renamed and extended: it
defines the groups and now also applies them. One `setup()`, since splitting
"defines" from "applies" would give two modules one caller each and put two setup
calls in `plugin/treesitter.lua` for one feature. `colors_spec.lua:277` requires
`chem/palette` and calls `setup()`, so the rename is a two-line test edit.

The painting logic factors out as `element_highlights(buf, root, srow, erow)`
returning `{srow, scol, erow, ecol, group}` — named for what it returns, and taking
a row range because the unit of work is a `Range6`, not a line.

**`after/queries/smarts/highlights.scm`** — keeps bonds, the reaction arrow, the
recursion anchor and the four structural hydrogen patterns. Two comments point at
the deleted layer: the header, and the line saying "the generated layer reaches the
`[#1]` spelling". Hydrogen stays split — `[H]` matched structurally in the query,
`[#1]` falling out of the lookup, both resolving to `@chem.element.hydrogen`.

## Cost of the persistent design

The first `on_changedtree` after attach arrives with the sentinel whole-buffer
range (`highlighter.lua:261`), so the initial pass is proportional to file size,
not window size. That is the real trade against a per-redraw design, and the reason
the `largefile` gate matters rather than being a nicety.

## Tests

Persistent extmarks are readable with `nvim_buf_get_extmarks`, so element colour is
directly assertable — the pure-function factoring is for clarity, not for reaching
otherwise-untestable state.

`tests/plenary/smarts_spec.lua`:

- `chem_groups` reads capture names and stops reporting element identity. Replace
  with extmark reads, keeping the existing cases: `[Fe]`/`[#26]` agreeing, `Sn`
  being sulfur plus nitrogen, bracket aromatics, suspects, `[H]` versus `[CH3]`.
- "gives every glyph at most one element group" goes **vacuous** — it would assert
  over a set containing no `chem.element.*` at all. Move it to the extmark path or
  delete it; keeping it as-is is worse than either.
- The drift check reads the suspect set from `chem/elements` instead of parsing the
  query file.
- Keep the markdown-fence case: it is the reason the wiring is recursive.

`tests/plenary/colors_spec.lua` — unaffected under the separate-export choice; its
three loops keep seeing only coloured elements.

One case worth a test rather than an assumption: a markdown buffer with no
`smiles` fence, and so no smarts tree at all — the ordinary state of most markdown
files, not an edge case.

## Rejected

**A decoration provider with ephemeral extmarks**, which an earlier draft of this
plan chose. It recomputes per redraw something that depends only on buffer text —
wrong in kind, not merely expensive, since making it 2 ms instead of 102 ms
optimises work that should not repeat. mini.hipatterns solves the same problem
(dynamic groups from buffer content, both unbounded for hex colours and a fixed
text→group map for TODO/NOTE) with persistent extmarks and no decoration provider.
Ephemeral marks are also invisible to `vim.inspect_pos`/`:Inspect`, which the
persistent design keeps.

**A node type per element in the grammar** (`(iron) @chem.element.iron`, no
predicate) — the escalation the previous plan named. O(matches) and no Lua, but it
puts chemistry in the parser, doubles the token set to cover `#26` as well as `Fe`,
inflates the parse tables, and changes the tree shape every existing query and test
is written against.

## Stale after this

`.claude/skills/treesitter-config/SKILL.md` § SMARTS / SMILES: the three-layer
table, the claim that the generator "writes both the middle layer and
`lua/chem/elements.lua`", the drift sentence naming the generated query, and the
warning not to collapse the per-element node-type lists with its 7× figure.

`lua/chem/palette.lua:1-10` opens "Colours for the `@chem.*` captures of the smarts
queries"; `plugin/treesitter.lua:16-18` justifies the setup call in terms of
element captures reaching a fence. Both describe the mechanism being removed.

`notes/PLAN-smarts-element-generator.md` reaches the wrong conclusion in its Costs
section. Retire it, or correct that section and let it stand as the record of how
the colours are computed — still accurate, and not restated here.
