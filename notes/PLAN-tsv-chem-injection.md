# PLAN: inject the smarts parser into a TSV chemistry column

## Goal

Highlight the cells of a TSV column whose header names a chemical line notation
(`smiles`, `cxsmiles`, `smarts`, …), leaving every other column alone. The
`smarts` grammar and both of its highlight layers already work, and the `tsv`
parser is installed; what is missing is the injection and its column gate.

## Why a host parser, and not a cell-at-a-time string parse

`lua/chem/highlight.lua`'s `attach()` walks the buffer's parser for `smarts`
language trees, registers `on_changedtree` on each, and catches later ones
through `on_child_added`. A decoration provider parsing each visible cell with
`get_string_parser` would have to reimplement that, so the host parser is the
cheaper route — but see § The element layer needs one change.

## What the parser already does (measured, config loaded)

| | |
| --- | --- |
| `ts_highlight` | `true` — the highlighter is active on every TSV buffer |
| `syntax` | `''` — `syntax/tsv.vim` no longer applies |
| `vartabstop`, `vim.b.tsv_widths` | set as before — `ftplugin/tsv.lua` is unaffected |
| `foldlevel` | `0` on every line — no `folds.scm` for `tsv`, so no folds appear |
| captures | `string/tsv` on every text cell, plus `number`, `number.float`, `boolean` |

Two riders worth recording rather than acting on. `indentexpr` is now
nvim-treesitter's for `tsv` (`lua/options.lua:34`) where it previously returned
`-1` for want of a parser; neither `tsv` nor `smarts` ships `indents.scm`, and
TSV has no leading whitespace, so this is inert. And losing `syntax/tsv.vim`
costs nothing visible: a `#` line already resolved to csv.vim's `csvCol0`, which
is undefined, so the ftplugin's `Comment` match was dead before. Do **not** add
`tsv` to `additional_vim_regex_highlighting` to get it back; `syntax/tsv.lua`
(which only clears `csvCol1..11`) is now dead code and can go separately.

Add `tsv` to the install list in `lua/plugins/treesitter.lua` so this is
reproducible rather than resting on a manual `:TSInstall`, and have the new spec
skip when `vim.treesitter.language.add("tsv")` fails, as the smarts specs do.

## Grammar shape (verified against the pinned `f6bf6e3`)

```
document → row → field → choice(text, number, float, boolean)
```

- `field` has exactly one named child with an **identical extent**.
- An empty cell **does not node at all** (measured, parser as installed):
  `x<TAB><TAB>y` is two fields, a trailing `<TAB>` adds none, and a row whose
  first cell is empty merges into the row above (its `\n` lands inside the
  previous row's last field). So there are no zero-length cells to reject, and
  sibling counting cannot give a column index — see § Column index.
- `text` has a quoted alternative, `seq('"', repeat(choice(/[^"]/, '""')), '"')`.
- There is no header and no comment concept: both are ordinary `row`s.

## Suppressing only the host's string capture

Keep the numeric captures, drop `(text) @string`. A non-extending
`queries/tsv/highlights.scm` in the config root shadows the shipped query:

```scm
; Deliberately no (text) capture: a text cell's colour is the injected
; language's business, or nothing.
(number) @number
(float) @number.float
(boolean) @boolean
```

**Not** under `after/` — measured, with a scratch rtp entry against nvim's own
lua queries:

| where the non-extending file sits | resulting query |
| --- | --- |
| rtp entry *before* the shipped one | the file is the whole query |
| rtp entry *after* it (incl. `after/`) | the file is ignored |

The rule is *first* non-extending file in runtimepath order wins as the base and
every later one is silently dropped — not "last wins". What the config root beats
is not nvim-treesitter's own directory, which is not on the runtimepath at all:
installing a parser symlinks `~/.local/share/nvim/site/queries/<lang>` to the
plugin's `runtime/queries/<lang>`, and the site directory (index 2) sits below
`~/.config/nvim` (index 1). Same position as `queries/smarts/atoms.scm`,
`queries/zsh/`, `queries/pymol_select/`, `queries/pymol/` and
`queries/karabiner/`.

Two non-options: an `; extends` file cannot cancel a capture, and linking
`@string.tsv` to `@none` only changes what the group looks like.

Keep a comment line in the file. A 0-byte query makes
`vim.treesitter.query.get` return `nil` rather than an empty query — both are
handled, but an apparently-empty file invites deletion.

Side effect to note: `runtime/queries/csv/highlights.scm` is `; inherits: tsv`
plus one `","` pattern, and the inherited base resolves through the same rtp
search — so a `csv` buffer keeps the delimiter and the numeric captures and
loses `@string` on its text cells too. Not observable yet: the csv parser is not
installed, and `query.get("csv", "highlights")` is `nil` without it.

## The column gate

The injected language is always `smarts` — one parser serves both notations — so
this is a **predicate** plus the built-in `#set!`, in the style of
`#command-is?` and `#arg-after?`, not an `inject-*` directive:

```scm
((field (text) @injection.content)
 (#chem-column? @injection.content)
 (#set! injection.language "smarts"))
```

Capture the `(text)` child, **not** `(field)`. `get_node_ranges`
(`languagetree.lua:900`) masks named children out unless
`injection.include-children` is set, and `field`'s single child has an identical
extent, so `(field) @injection.content` yields **zero** ranges;
`add_injection` then returns early on `#ranges == 0` with no diagnostic — a
silent no-op. Capturing `text` also drops `number`/`float`/`boolean` cells for
free, and having no named children means `metadata[capture].range` is honoured
without `include-children`.

A failing predicate drops the whole match, so this never relies on the
omit-`injection.language` path (which does work — `languagetree.lua:1128` logs
and skips — but is the weaker signal).

`#chem-column?` answers one question: is this cell in a chemical column? It:

1. Skips the header row and any comment row — the grammar has no notion of
   either, so `# generated by rdkit` is a `row` whose first field is `(text)`
   and would otherwise inject.
2. Takes the column index from the node's preceding `field` siblings.
3. Accepts the index if it is marked by hand, or if the header names it (below).
4. Rejects a cell holding `""`, which is one escaped quote that no range can
   undouble.

Deliberately **not** in the predicate: quote stripping. That is `#unquote!`, a
one-job directive beside `#trim!`/`#head!`/`#tail!` in `lua/utils/treesitter.lua`
— `#trim! @cap 1 1` conditioned on the quotes being there, so one pattern serves
both cell shapes rather than two gated on `#lua-match?`/`#not-lua-match? "^\""`.

## Column index: tabs in the line prefix, not `field` siblings

The two diverge on an **empty cell**, which the grammar drops (above) and which
occurs in any table with missing values: every column right of it would shift.
So the index is the tab count of the line up to the cell's start byte — the same
measure `ftplugin/tsv.lua:286` takes, and the same tab-free-cell assumption it
already declares, but no shared code: nothing here calls the ftplugin, so there
is no contract to agree on and no dependency to sequence.

Rows still come from the tree — which row is the header, and which rows are
comments — because that is what a parse is good for.

Rows spanning lines are likewise ruled out by practice, so the header can be
identified by line number.

Guard against a header shorter than its data rows — R's `write.table` omits the
rowname column's name — by comparing the header's field count with the first
data row's and refusing to inject on a mismatch. Otherwise every index is off by
one.

## Header detection

The header is the first non-comment row, read from the tree. `.bed` also resolves
to filetype `tsv` and is headerless; the grammar's own `(number)`/`(float)`
nodes make a numeric field in row 1 good evidence of that, using a parse that has
already happened.

Which names count as chemical belongs in `lua/chem/notation.lua` returning
`"smiles"|"smarts"|nil`, not a boolean: `notes/PLAN-chem_draw.md:81-83` needs the
notation itself so that a `reaction_smarts` column selects SMARTS
automatically. The predicate only tests for non-nil.

An **exact set**, not a pattern: `smiles`, `cxsmiles`, `smarts`, `smirks` and the
`reaction_`/`product_`/`reactant_` spellings worth naming outright. A pattern
would catch the rest for free but eventually matches a column that is not one,
and the hand-marking below covers the misses at no risk.

## Marking a column by hand

A keymap that marks the column under the cursor as chemical, for a header the set
does not name. This is the same on-demand philosophy as
`notes/PLAN-chem_draw.md` — the header set is a convenience, not the mechanism.

The mark is per buffer and per column index, and the cursor's column index comes
from the `field` node under it, so the keymap and the predicate share one
definition and neither touches the ftplugin. Store it densely (a list of
booleans, or a table keyed by bufnr in the module) — **not** a sparse table in
`vim.b`, which reads back with truthy `vim.NIL` in the holes.

Marking has to force a reparse for the injection query to run again:
`parser:invalidate(true)`, as `ftplugin/python.lua`'s `Load_pymol()` already does
after changing what a query matches.

The mark says "this column holds a chemical structure" and takes no language
argument. Both notations parse with the one `smarts` parser, so the injected
language is constant; a "highlight this column as `<language>`" keymap would be a
different, currently unused feature, and it would force the gate back into a
directive resolving a variable language.

## The element layer needs one change

It does not come entirely for free. `lua/chem/highlight.lua` cleared with
`nvim_buf_clear_namespace(buf, ns, srow, erow)` — row-granular — over whole-row
tree bounds. Two chem columns on one row therefore wipe each other's marks, and
which survives depends on `pairs()` order over the trees. So `paint` clears the
**range** (`nvim_buf_get_extmarks` over it, then `nvim_buf_del_extmark`) and
`repaint` clamps range-wise. This is an owned change, not a consequence.

Two things the row-granular version got for free and the range one has to state:

- The repainted part of a region is the change's **rows** intersected with the
  region's **columns**. A mark moves with the buffer edit, not with the bytes
  tree-sitter calls changed: `nvim_buf_set_lines` over a line collapses every
  mark on it to the start of the next line — past the last line when the line
  replaced was the last, where only a sweep to the buffer's end finds them
  again. Columns-from-the-region is what keeps a sibling cell's marks.
- `on_child_removed` has to clear off the removed tree's `included_regions()`,
  not its trees: `invalidate(true)` — which the mark keymap calls — empties
  `_trees` before the reparse that drops the language, so the regions are the
  only record left of where the marks were.

## Interaction with the column hiding

`ftplugin/tsv.lua` hides a column by removing its text from the buffer and
stashing it in `vim.b`, so a hidden chem column would inject over truncated text
and fill with ERROR nodes. The predicate has to skip it — the one place this work
touches ftplugin state.

Two traps there. `getCommentChar`, `isComment` and the `vim.b.tsv_*` accessors
are file-locals, so they cannot be called from a predicate; and a sparse table
round-tripped through `vim.b` comes back with `vim.NIL` in the holes, which is
**truthy** — `if hidden[col] then return end` would skip every column left of the
hidden one. That is what `denilify` exists for. Since the ftplugin is heading for
a shareable plugin under `modules/`, the accessor belongs in that plugin's public
API; until then this needs an explicit, denilified read.

`hide()` does one `nvim_buf_set_text` per row, and `BufWritePre`/`BufWritePost`
unhide then rehide the whole buffer on every save — each an injection re-run.

## Cost

`notes/PLAN-smarts-element-perf.md` is the standing warning: a query layer that
looked cheap cost 102 ms per redrawn window. Here the dominant cost is not the
predicate but **region churn**: when the region count changes,
`set_included_regions` (`languagetree.lua:846`) discards every tree and fires
`changedtree` for all of them, so each scroll reparses and repaints the whole
viewport's regions. Measured on a comparable host, ~160 callbacks per viewport
parse in ~2 ms, ~4 ms with the element layer attached.

`injection.combined` is not a mitigation: besides concatenating cells into one
bogus record, it forces a full-document scan on every injection pass
(`languagetree.lua:1102`).

## Tests

`tests/plenary/tsv_spec.lua` exists but currently only covers `vim.b.tsv_*`
round-tripping. Add:

- header → index mapping, and a `#` comment row above the header
- a comment row at the chem column's index does not inject
- a headerless `.bed` first row is not read as a header
- injection fires on the named column only
- a quoted cell injects without its quotes
- two chem columns in one row both inject, and both keep their element marks
- a hidden chem column does not inject, and its element marks are cleared
- a header with fewer fields than its data rows refuses to inject
- text cells are no longer `@string`, while numbers and floats still capture
- hand-marking a column whose header is not in the set makes it inject, and the
  mark survives an edit to an unrelated row

## Later

- rdkit call sites in Python — already noted in `notes/TODO.md`.
- `.csv`: the `csv` parser's `requires` pulls in `tsv`, and the predicate would
  need to stop assuming the separator.
- The ftplugin's promotion to a shareable plugin under `modules/`. This work is a
  consumer of its eventual API (the hidden-columns accessor), not a driver of it.

## Accepted limitations

A headerless SMILES TSV (`CCO<TAB>ethanol`, no numeric field anywhere in row 1)
defeats the `.bed` heuristic: row 1 is read as a header, `CCO` becomes a column
name, and nothing injects. Accepted — hand-marking the column is the workaround,
and a `.bed`-style headerless file with a structure column is not a shape that
occurs here yet.

Empty cells are **not** among them: with the index and the header both counted in
tabs off the line, a missing value costs nothing — `ethanol<TAB><TAB>CCO` injects
`CCO` at index 3, and a row whose first cell is empty still injects its own
cells at their own indices, even though the grammar folds it into the row above.

## Note

`largefile` (1 MB, `lua/largefile.lua:13`) gates treesitter off entirely, so
large TSVs get none of this. That is existing cross-filetype policy and not this
plan's to change.
