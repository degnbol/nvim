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

--- Copy of an image's extmarks, reshaped to fill a `width`-cell box without
--- conceal. Each image row (an `overlay` or `inline` mark) becomes an overlay,
--- loses `conceal` and is padded with blank cells to `width`. A `conceal_lines` mark becomes one blank
--- overlay row per line it covered, at the first overlay's `virt_text_win_col`.
--- All overlay rows get `hl_mode = "combine"`, so the line bg shows under image
--- and padding. Other extmarks (e.g. `virt_lines`) pass through unchanged.
--- @param extmarks table[] extmark specs as passed to nvim_buf_set_extmark, plus `row`/`col`
--- @param width number box width in cells
--- @return table[] extmarks
function M.fill_box(extmarks, width)
    local win_col
    local boxed = {}
    for _, e in ipairs(vim.deepcopy(extmarks)) do
        if e.virt_text_pos == "overlay" or e.virt_text_pos == "inline" then
            win_col = win_col or e.virt_text_win_col
            e.virt_text_pos = "overlay"
            e.conceal = nil
            e.hl_mode = "combine"
            local image_width = 0
            for _, chunk in ipairs(e.virt_text) do
                image_width = image_width + vim.fn.strdisplaywidth(chunk[1])
            end
            if width > image_width then
                table.insert(e.virt_text, { (" "):rep(width - image_width) })
            end
            boxed[#boxed + 1] = e
        elseif e.conceal_lines then
            for row = e.row, e.end_row do
                boxed[#boxed + 1] = {
                    row = row,
                    col = 0,
                    virt_text_pos = "overlay",
                    virt_text = { { (" "):rep(width) } },
                    hl_mode = "combine",
                    virt_text_hide = false,
                    virt_text_win_col = win_col,
                }
            end
        else
            boxed[#boxed + 1] = e
        end
    end
    return boxed
end

--- Widest display width of a span's lines right of its start column, or
--- nil when other text shares those lines: non-blank text left of the start
--- column on any line, or after the end column on the last line. The start
--- column is a byte column, and it is used as a display column on every line.
--- @param lines string[] buffer lines the span covers
--- @param start_col number 0-indexed byte column of the span's start on its first line
--- @param end_col number 0-indexed exclusive byte column of the span's end on its last line
--- @return number|nil width in cells
function M.footprint_width(lines, start_col, end_col)
    if lines[#lines]:sub(end_col + 1):find("%S") then
        return nil
    end
    local text_width = 0
    for _, line in ipairs(lines) do
        if line:sub(1, start_col):find("%S") then
            return nil
        end
        text_width = math.max(text_width, vim.fn.strdisplaywidth(line))
    end
    return text_width - start_col
end

--- Copy of `virt_lines` drawn in the line highlight `hl`: each chunk's group
--- becomes `{ hl, group }` (`hl` for a chunk without one), and each row ends
--- with a chunk of `pad` spaces in `hl`.
--- @param virt_lines table[] rows of `{ text, group }` chunks, as for nvim_buf_set_extmark
--- @param hl string highlight group
--- @param pad number cells of the trailing chunk. With 0, neovim 0.13+ extends
--- its hl to the window edge (neovim/neovim#41289).
--- @return table[] virt_lines
function M.with_line_hl(virt_lines, hl, pad)
    return vim.tbl_map(function(line)
        local chunks = vim.tbl_map(function(chunk)
            return { chunk[1], chunk[2] and { hl, chunk[2] } or hl }
        end, line)
        chunks[#chunks + 1] = { (" "):rep(pad), hl }
        return chunks
    end, virt_lines)
end

return M
