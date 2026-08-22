# PLAN: scrolling a buffer whose injections are per-row

Scrolling `rxns.tsv.gz` (4568 rows, 8.3 MB, a `cxsmiles` column injected into the
smarts parser) costs 25–28 ms per line. One row enters the window and every
injected cell on screen is reparsed and repainted.

## The cost

Latency measured in a real TUI: nvim in a pty, one `<C-e>` written, timed to the
last byte of the frame that came back. 220x60, warmed over the same rows first,
median of 20–25 keys.

| | median | p90 |
| --- | --- | --- |
| as-is | 25.4 ms | 31.4 |
| `chem/highlight` stubbed out at startup — no element marks | 20.1 ms | 27.9 |
| `queries/tsv/injections.scm` emptied — no smarts child at all | 1.8 ms | 2.5 |

So the smarts injection is ~18 ms and the element extmarks ~5 ms of every
scrolled line. Extmark *drawing* is not involved: with the marks in place but
every `@chem.element.*` group emptied, the frame is 298 474 bytes against 302 935,
and the redraw time does not move.

Counting extmark IDs across one `<C-e>`, after warming:

```
warm             marks=4129  regions=61  trees=60
after scroll 1   marks=4254  kept=327   new=3927
after scroll 2   marks=4316  kept=382   new=3934
```

Nearly the whole visible set is destroyed and re-created per line.

## Why

`runtime/lua/vim/treesitter/languagetree.lua:846`, `set_included_regions`:

```lua
if #self:included_regions() ~= #new_regions then
  -- TODO(lewis6991): inefficient; invalidate trees incrementally
  for _, t in pairs(self._trees) do
    self:_do_callback('changedtree', t:included_ranges(true), t)
  end
  self._trees = {}
  self:invalidate()
else
  self:_iter_regions(function(i, region) return vim.deep_equal(new_regions[i], region) end)
end
```

A region is identified by its index in the list. The highlighter parses the
visible range, so the list *is* the injected cells currently on screen: scrolling
one line changes its length, or shifts every entry by one. The first branch drops
every tree and reports every one of them as changed. The second invalidates each
shifted region and hands `parser:parse` the tree of the region that *used to* sit
at that index (`languagetree.lua:436`), so the incremental diff runs between two
unrelated cells and the whole region comes back changed anyway.

Unchanged on neovim master as of 2026-08-22 (same branch, refactored into
`_do_changedtree_callbacks()`; checked through the GitHub API), and 0.12.4 is
current stable, so there is nothing to upgrade to.

## Process injections over a block, not the viewport

What stops the churn is not a complete region list but
`_processed_injection_region` *containing* the requested range
(`languagetree.lua:680-687`), which skips `_get_injections` and
`set_included_regions` outright. Any range gives that. Measured on the fixture,
headless with `vim.g._ts_force_sync_parsing`:

| block | one-off | marks | 10 scrolls inside it | edit + parse |
| --- | --- | --- | --- | --- |
| 300 lines | 171 ms | 19 744 | 0.18 ms | 36 ms |
| 600 lines | 326 ms | 40 247 | 0.12 ms | 35 ms |
| 1200 lines | 690 ms | 89 801 | 0.31 ms | 49 ms |
| whole file | 2630 ms | 306 058 | 0.0 ms | **2384 ms** |

The whole-file column is why the block is bounded rather than "parse everything
once": a full parse makes the *next keystroke* reprocess injections over the
visible range, `set_included_regions(41)` against 4567, drop every tree, and run
4567 repaints. It also exceeds the async budget — `_async_parse` caps *cumulative*
parse time at `'redrawtime'` (2000 ms, `languagetree.lua:596-604`), and the full
parse costs 2630 ms with the painter attached, so it aborts with `err=TIMEOUT`
having possibly done nothing.

Where it goes: how much of a buffer's injections to process is buffer-level
treesitter policy — it changes what the highlighter, folds and `:InspectTree` see,
not just element colours — so the trigger and its gate belong in
`plugin/treesitter.lua` beside `start_treesitter` and the `largefile` gate, and
the "process injections over a range around the viewport" helper in
`lua/utils/treesitter.lua`.

Open:

- **Block size.** 300 lines buys the whole win on this fixture. What the right
  number is depends on the cost of crossing an edge (164 ms for a 300-line block),
  which sets the floor.
- **The re-block trigger.** `WinScrolled` within N lines of the edge, scheduled so
  it does not land inside the scroll it is reacting to. Untested.
- **Whether painting follows the parse.** Processing a block paints every atom in
  it, displayed or not — 19 744 marks for 300 lines. Painting only what is visible
  would need its own paint path on scroll, and the marks are the cheap half.
- **p90 stays at 27.6 ms** after a whole-file parse while the median drops to 6.0.
  If the marktree size is the cause, the block bounds it; one measurement with the
  painter stubbed settles it.

## Repaint only what differs

`paint` clears a range and re-adds. When a region is reported changed but its text
did not change — every scroll — the marks it removes and the marks it adds are the
same set. Compare instead of clearing, and return when they match.

Two requirements, both from what `clear` currently does beyond the range:

- Any fetched mark whose `details.invalid` is set forces the repaint. `clear` also
  removes invalid marks past the range end, which is how a mark survives
  `nvim_buf_set_lines` collapsing a row; a comparison built from
  `element_highlights` alone cannot see them.
- Compare sorted multisets of `(srow, scol, erow, ecol, hl_group)`.
  `nvim_buf_get_extmarks` returns traversal order and `iter_captures` match order,
  and identical signatures can repeat in one range.

Comparing rather than caching a hash of the text: a cache would need invalidating
everywhere marks are cleared, which is the second source of truth that produced
the mark leak this module was already fixed for once.

Worth the delete-and-re-add half of `paint` — measured over 10 scroll parses on
the fixture, `del_extmark` 2.2 ms and `set_extmark` 5.0 ms per scroll. The
`get_extmarks` call stays and grows a little (it needs `details`), and the
`atoms.scm` query stays.

## What each fixes

Blocking stops the invalidation; comparing stops paying for it when it happens
anyway — on an edit, which reprocesses injections over the changed range and
changes the region count again, and on `toggle_column`'s `invalidate(true)`.
Neither substitutes for the other: comparing cannot touch the reparse that is 73%
of the cost, and blocking does not reach the edit path.

## Tests

The comparison is assertable in `tests/plenary/smarts_spec.lua` beside the
existing mark cases: force a redundant reparse over unchanged text
(`parser:invalidate(true)` then `parse()`) and assert the extmark **IDs** are
unchanged — the existing cases assert positions, which pass either way. That path
now reaches a clean starting state: `clear`'s sentinel-row clamp and
`toggle_column`'s stale marks were both fixed first, with their own regressions.

Blocking is a latency property and the suite cannot time a redraw, but the
mechanism is assertable: hook `set_included_regions` and assert zero calls across
a scrolled range parse inside the block, or assert `parser:is_valid(nil, range)`.
Set `vim.g._ts_force_sync_parsing` so the parse finishes rather than racing the
test.

## Known gap

A region set that shrinks while the language survives strands the dropped
region's marks — nothing reports its rows again. `toggle_column` now handles this
by dropping the buffer's marks and letting the reparse repaint; an edit that makes
a cell stop qualifying (a doubled quote appears, the header is renamed) does not.
