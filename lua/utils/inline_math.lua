local M = {}

--- True cell width of a collapsed inline-math image at conceallevel 2+.
--- The PNG stretch-fills an integer `width × 1` cell box, so the only knob is
--- the whole-cell width. Round UP, never to nearest: a box narrower than the
--- glyph's native width compresses it (`$T$` is ~1.2 cells wide, and
--- `floor(native+0.5)` squashed it to 1). Over-width is just whitespace;
--- compression distorts the glyph, so the box must never be smaller than
--- native. `$k$` (~0.9) and `$k_{cat}$` (~2.7) don't sit in the round-down
--- zone, which is why only `$T$` looked wrong.
--- @param png_w number PNG pixel width
--- @param png_h number PNG pixel height
--- @param cell_w number terminal cell pixel width
--- @param cell_h number terminal cell pixel height
--- @return number cells
function M.cell_width(png_w, png_h, cell_w, cell_h)
    local native = png_w / png_h * (cell_h / cell_w)
    return math.max(1, math.ceil(native))
end

return M
