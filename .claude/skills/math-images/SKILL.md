---
name: math-images
description: Rendered LaTeX math images in neovim markdown/tex buffers (`$…$`, `$$…$$`, snacks.image). Use when working on their size, baseline, conceal or layout.
---

# Math Images

LaTeX math (`$…$`, `$$…$$`) renders by compiling the expression to a PNG and
placing it in the buffer through the terminal graphics protocol (kitty unicode
placeholders). The active backend is **snacks.image**, but most of the hard
parts are inherent to *any* such backend in neovim.

## Design intent: tier the render by `conceallevel`

A design choice, not a backend feature — any backend should honour it:

| `conceallevel` | Source `$…$` | Render | Layout intent |
|---|---|---|---|
| **0** | shown literally | **none** | editing — see the raw LaTeX |
| **1** | keeps its footprint (image overlaid) | inline: centred in the source footprint; display block: native size | **no reflow**: text never moves as the cursor enters/leaves the line |
| **2+** | fully hidden | inline: **true glyph size**, no surrounding whitespace; display block: snacks' default | reflows around the real footprint |

Intended consequences (not bugs):

- At **cl=1** an image wider than its source is fit to the source width, so it
  is drawn smaller and off the baseline. That is the cost of zero reflow; cl=2
  trades it for true size at the cost of reflow.
- At **cl=1** a display block keeps native size instead. Blank cells cover the
  rest of the source; a taller image continues in virtual lines, so nothing
  reflows.
- At **cl=2+** every inline glyph renders at a consistent size regardless of the
  expression around it (standalone `$k$` matches the `k` in `$k_{cat}$`), with
  no padding on either side.

## Backend-independent gotchas (any inline-image renderer)

- **Cell-box fit.** Kitty fits a unicode-placeholder image into its integer
  `width × height` cell box with the aspect ratio kept: to the box height,
  centred horizontally, or, if the image is wider, to the box width, centred
  vertically (kitty `graphics.c` `grman_put_cell_image`). A box narrower than
  the image makes it smaller; a wider box adds centred whitespace. Native size =
  box width equals the rounded-up true cell width.
- **Compute true aspect from raw pixels, not pre-ceiled cell sizes.** A helper
  that converts px→cells typically `ceil`s *both* axes first; taking a ratio of
  those double-rounds and shrinks wide expressions. Go back to the raw PNG px
  (`img.info.size`) and the terminal cell dimensions:
  `ceil(png_w/png_h · cell_h/cell_w)`.
- **No baseline awareness.** A backend that fits the PNG to the cell height
  maps the box *bottom* to the line *bottom* — it knows nothing about the text
  baseline. A descender-less glyph then hovers above the line. Fix in the math
  document with a strut: its height sets the scale (box = one cell) and its depth
  fraction the baseline. Only fractions of the cell matter, so the fit is
  zoom-independent. `latex.cell_fit_tex` and `typst.cell_fit_typ` compute it
  from the terminal's cell layout (`kitty.cell_fractions`): strut depth = the
  cell's space below the baseline; font size so that math x-height, cap-height
  and descender match the prose with the smallest worst-case error. An expression that extends beyond
  the strut grows its own PNG and is drawn smaller (subscripts shrink ~3-11%,
  as the cell has little space below the baseline). unicode-math picks script-style
  glyphs by pt size relative to the size at `\setmathfont`, so `cell_fit_tex`
  reruns each `\setmathfont` at the fitted size (~9.5 pt).
- **Trim crops the strut.** Once you pad with a transparent strut, image-convert
  trimming (e.g. `-trim`) crops it straight back off. Trim must stay off, so
  PNGs carry the strut's vertical padding by design.
- **Conceal overrides line highlights.** At cl=1 a concealed `$…$` leaves one
  cell drawn with the `Conceal` attr *replacing* the line attr (sign `linehl`,
  `line_hl_group`, diff) — neovim/neovim#31555. Instead leave the source
  unconcealed and cover it with overlay cells (image, plus blank padding for
  blocks) with `hl_mode="combine"` (`"replace"` paints `Normal` bg). Overlays
  clip at the window edge, so a span crossing a screen-line wrap must fall back
  to inline+conceal. Layout thus depends on conceallevel *and* wrap geometry;
  re-render when either changes.
- **virt_lines get no line attr.** Image rows continuing in `virt_lines` lose
  the source line's `line_hl_group`: combine it into every chunk and end each
  row with a chunk in that group reaching the window edge. Re-render when the
  line's hl changes (e.g. gitsigns re-diff).

## snacks.image specifics

snacks runs every image type through one fit-into-a-cell-box pipeline and, for
math, **discards the `$`/`$$`/`\[` delimiter and always wraps display `\[…\]`** —
so inline vs block is indistinguishable downstream and display glue makes even
one glyph's box ~3 cells tall. Config alone can't fix it; the work is
overrides on exported snacks tables, applied after `setup` so they survive
`vim.pack` updates. All live in `pickers.lua` with the full *why* in comments —
read those when touching this:

- PNG size: snacks shows `px / dpi · 96 · scale` px, `scale = cell_w / 8`
  (image/util.lua `fit`), so a PNG `aspect / 12` in tall is one cell at any zoom
  (`image_placement.row_height_in`). One `math_look` table sets the fonts of
  both engines (`fonts_tex`/`fonts_typ`, before `${header}`). The cell fit goes
  after `${header}` in `math.latex.tpl` and `math.typst.tpl`, set on the first
  math render. Typst: the fit drops snacks' 2 pt margin and starts each
  equation with a strut box; snacks' `bounds` text edges then make an inline
  PNG the union of strut and ink (`utils/typst.lua`).

- `doc.transforms.latex` — inspects raw `img.content` *before* snacks strips the
  delimiter. Except in `.tex` buffers and buffers with a snacks header, it
  rewrites the math to typst via mitex (`typst.mitex_typ`, `img.ext =
  "math.typ"`), ~0.06 s against tectonic's ~2 s. On tectonic it rewrites
  inline `$…$`/`\(…\)` to `\begin{math}<strut>…\end{math}` (drops display
  glue, floors to one line height); display math likewise in `\displaystyle`,
  unless the body starts with `\begin`. Tags `img.inline`, which
  `inline.update` passes on as the `display` placement opt.
- `doc.find_visible` — filters out `type == "math"` matches when the window is
  at cl=0, so the inline manager closes them and the source shows.
- `placement.state` — sets `loc.width` per the tier table: cl=1 source display
  width + custom `loc.overlay` flag (width−1, inline, when `utils.crosses_wrap`
  reports a wrap crossing); cl=2+ true cell width from raw px. Don't reuse
  snacks' `ceil(w/h)+2` (placement.lua) — its `pixels_to_cells` pre-ceils both
  axes (see above). `display` math at cl=1 with a clean rectangular footprint
  (`footprint_width`) sets `loc.box_width`, capped to fit the text area.
- `placement._render` — on `loc.overlay`, turns the inline extmark into an
  unconcealed `overlay`, `hl_mode="combine"`. On `loc.box_width`, `fill_box`
  (`utils/image_placement.lua`) does the same to every image row, pads it to
  the box and swaps conceal for blank overlay rows. Any placement's virt_lines
  go through `with_line_hl`. A decoration-provider `on_end` re-renders
  placements whose line hl went stale.
- An `OptionSet` autocmd (conceallevel + wrap-geometry options) re-fires the
  inline manager's own `BufWinEnter` handler (its `snacks.image.inline.<buf>`
  augroup) to re-render.
