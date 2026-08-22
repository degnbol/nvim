---
name: math-images
description: Inline LaTeX math rendering in neovim (LaTeX→PNG via a terminal-graphics image backend) and its conceallevel-tiered sizing/baseline gotchas. Use when working on inline `$…$` math images — width, size, baseline, or conceal behaviour in markdown/tex buffers. The active backend is snacks.image, configured in lua/plugins/pickers.lua.
---

# Inline Math Images

Inline LaTeX `$…$` / `$$…$$` renders by compiling the expression to a PNG and
placing it in the buffer through the terminal graphics protocol (kitty unicode
placeholders). The active backend is **snacks.image**, but most of the hard
parts are inherent to *any* such backend in neovim. The gotchas below are split
accordingly: general ones first, then the snacks-specific section.

## Design intent: tier the render by `conceallevel`

This is the level-1-vs-2 "conceals harder" distinction, made meaningful for
math. It is a design choice, not a backend feature — any backend should honour
it:

| `conceallevel` | Source `$…$` | Render | Layout intent |
|---|---|---|---|
| **0** | shown literally | **none** | editing — see the raw LaTeX |
| **1** | one residual cell | fills the source footprint | **no reflow**: text never moves as the cursor toggles conceal on the line |
| **2+** | fully hidden | **true glyph size**, no surrounding whitespace | reflows around the real footprint |

Intended consequences (not bugs):

- At **cl=1** a wide expression (`$k_{cat}$`) is stretched to fill its source
  width, so it looks "fat". That is the cost of zero reflow; cl=2 trades it for
  true size at the cost of reflow.
- At **cl=2+** every glyph renders at a consistent size regardless of the
  expression around it (standalone `$k$` matches the `k` in `$k_{cat}$`), with
  no padding on either side.

## Backend-independent gotchas (any inline-image renderer)

- **Cell-box stretch-fill.** The image fills an integer `width × height` cell
  box exactly (kitty unicode placeholders). The width you pick controls *size*,
  not padding — there is no "render at native size, left-aligned" mode. Native
  size = the box width *equals* the rounded true cell width.
- **Compute true aspect from raw pixels, not pre-ceiled cell sizes.** A helper
  that converts px→cells typically `ceil`s *both* axes first; taking a ratio of
  those double-rounds and squashes wide expressions. Go back to the raw PNG px
  (`img.info.size`) and the terminal cell dimensions:
  `round(png_w/png_h · cell_h/cell_w)`.
- **No baseline awareness.** A backend that stretches the PNG to fill a cell box
  maps the box *bottom* to the line *bottom* — it knows nothing about the text
  baseline. A descender-less glyph then hovers above the line. Fix on the LaTeX
  side with an asymmetric strut whose depth fraction matches the editor font's
  descender (~0.15). Two knobs, both in `\baselineskip` units: the **ratio**
  `depth : (depth+height)` is the baseline split (0.15 is the floor — shallower
  and a subscript drops below the strut and grows the box, so `$k$` and
  `$k_{cat}$` stop sharing a height); the **sum** is the size knob (glyph ink is
  fixed px, so a taller box shrinks the glyph once the box collapses to one
  row). Scale both together to resize without moving the baseline. Exact
  alignment only holds for simple glyphs; fractions/limits extend both sides.
- **Trim crops the strut.** Once you pad with a transparent strut, image-convert
  trimming (e.g. `-trim`) crops it straight back off. Trim must stay off, so
  PNGs carry the strut's vertical padding by design.
- **conceallevel changes the residual footprint.** At cl=1 the concealed
  `$...$` leaves one residual cell; at cl=2 it leaves none. Width therefore
  depends on the level, and most backends only recompute placements on
  scroll/edit/enter — so a level change needs an explicit re-render trigger.

## snacks.image specifics

snacks runs every image type through one fit-into-a-cell-box pipeline and, for
math, **discards the `$`/`$$`/`\[` delimiter and always wraps display `\[…\]`** —
so inline vs block is indistinguishable downstream and display glue makes even
one glyph's box ~3 cells tall. Config alone can't fix it; the work is three
overrides on exported snacks tables, applied after `setup` so they survive
`vim.pack` updates. All live in `pickers.lua` with the full *why* in comments —
read those when touching this:

- `doc.transforms.latex` — inspects raw `img.content` *before* snacks strips the
  delimiter; rewrites inline `$…$`/`\(…\)` to `\begin{math}<strut>…\end{math}`
  (drops display glue, floors to one line height). Block math passes through.
- `doc.find_visible` — filters out `type == "math"` matches when the window is
  at cl=0, so the inline manager closes them and the source shows.
- `placement.state` — sets `loc.width` per the tier table. cl=1 uses source
  display width minus the residual cell; cl=2+ computes the true cell width from
  raw px. Don't reuse snacks' own `ceil(w/h)+2` (placement.lua) — its
  `pixels_to_cells` already ceils both axes, the pre-ceiled bug above.
- An `OptionSet conceallevel` autocmd re-fires the inline manager's own
  `BufWinEnter` handler (scoped to its `snacks.image.inline.<buf>` augroup) to
  re-render when the level changes.
