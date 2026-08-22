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

## Processing injections over a block: tried, reverted

What stops the churn is not a complete region list but
`_processed_injection_region` *containing* the requested range
(`languagetree.lua:680-687`), which skips `_get_injections` and
`set_included_regions` outright. Asking the parser for a block of lines around
the viewport, on a grid so that every viewport inside one block asks for the same
range, gives that: the highlighter's later parses of the visible rows find the
injections already processed and touch nothing.

It works, and it is not worth having. Both halves of that are measured.

Working: `<C-e>` in the TUI from row 1000, 120 timed keys, median 15 ms against
27 and `set_included_regions` 1 against 100. Not worth having, for two reasons.

**The gain does not reach the input.** Scrolling is paced by the events, not by
the work. 300 wheel notches at a 16 ms cadence move the same 296 lines in the
same 4.8 s either way, because nvim coalesces redraws while input is pending and
absorbs a costlier frame by drawing fewer of them. What the block buys there is
~2× the frames drawn and a tenth of the injection processing — smoother and
cheaper, not shorter, and nothing a hand on a wheel can tell apart.

**It moves the cost rather than removing it, into a shape that is worse.**
Entering a block reparses every region of the range asked for, so the amortised
cost is flat in the block size and the size only sets how much lands at once —
and it lands whole if the scroll stops there. Scrolling 700 lines from row 1000
in a 60-row window, headless:

| block | per scrolled line | per block entered |
| --- | --- | --- |
| 60 lines | 0.72 ms | 40 ms |
| 100 | 0.72 | 73 |
| 150 | 0.59 | 100 |
| 300 | 0.75 | 318 |
| 600 | 0.88 | 612 |

At 150 lines, 200 `<C-e>` at a 16 ms cadence left the screen 0.54 s behind the
last key, against 0.07 s with no block at all — the transition landing after the
input stopped. Sizing the block to the viewport fixes that (0.00 s) and moves the
same cost into occasional 80–90 ms frames at the transitions instead. A failure
mode that depends on where a scroll stops relative to an invisible grid is worse
than a uniform cost, and it took an outside review to find at all.

The whole-file end of that table is why no block is large enough to escape this:
2630 ms one-off and 2384 ms for the next edit, which also exceeds the async
budget — `_async_parse` caps *cumulative* parse time at `'redrawtime'` (2000 ms,
`languagetree.lua:596-604`), so it aborts with `err=TIMEOUT` having possibly done
nothing.

## Where the time actually goes

The residual ~15 ms per drawn frame is not the reparse, not the element marks
(dropping them: 16.2 ms against 16.0) and not the smarts captures (emptying that
query: 16.7 ms). Against the 1.8 ms with `queries/tsv/injections.scm` emptied,
what is left is nvim's per-redraw cost of the ~60 injected trees a window holds
— by elimination, and independent of what any of them highlights.

So the only lever on scrolling this file is **fewer injected regions per screen**:
inject the cell under the cursor and no other, or drop the per-cell injection and
let the element marks carry the column, losing the bond, bracket and charge
captures with the trees. Both are design changes, not constants to tune, and both
are unmeasured.

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

It buys the delete-and-re-add half of `paint`: the parse that follows an edit
costs 66 ms clearing and re-adding against 58 ms comparing, over the same marks
(34 742 either way). That is what is left to win once the reparse itself is out of
reach — it does not touch the reparse, which is 73% of a scrolled line. The
`get_extmarks` call stays and grows a little (it needs `details`), and the
`atoms.scm` query stays.

The regression is in `tests/plenary/smarts_spec.lua` beside the existing mark
cases: a redundant reparse over unchanged text (`parser:invalidate(true)` then
`parse()`) leaves the extmark **IDs** unchanged — the existing cases assert
positions, which pass either way.

## Known gaps

A region set that shrinks while the language survives strands the dropped
region's marks — nothing reports its rows again. `toggle_column` now handles this
by dropping the buffer's marks and letting the reparse repaint; an edit that makes
a cell stop qualifying (a doubled quote appears, the header is renamed) does not.

Deleting the *last* atom of a structure leaves that atom's mark behind with
`invalid` set: no region reports those rows again, so nothing runs the comparison
that would notice. Deleting a middle atom and replacing a whole row both come out
right. Predates the comparison — `clear` alone never saw those rows either.
