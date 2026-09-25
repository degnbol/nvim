---@diagnostic disable: undefined-global
-- PNG dimensions below were measured by compiling each expression through the
-- same pipeline pickers.lua uses (strut-wrapped \begin{math}, tectonic, magick
-- -density 192, no -trim). Cell size 9×18 is snacks' headless fallback
-- (terminal.lua M.size). Regenerate if the strut or font_size changes.
local image_placement = require "utils/image_placement"
local cell_width = image_placement.cell_width
local fill_box = image_placement.fill_box
local footprint_width = image_placement.footprint_width
local with_line_hl = image_placement.with_line_hl
local CELL_W, CELL_H = 9, 18

-- (png_w, png_h) → native cell width = png_w/png_h * cell_h/cell_w
local glyphs = {
    ["$k$"] = { 17, 39 }, -- native ~0.87
    ["$T$"] = { 23, 39 }, -- native ~1.18 — round-to-nearest squashed this to 1
    ["$k_{cat}$"] = { 53, 40 }, -- native ~2.65
}

describe("inline math cell width", function()
    -- Conceallevel-2 sizing rule: the box must never be narrower than the
    -- glyph's native width, or the stretch-fill compresses the glyph.
    for src, png in pairs(glyphs) do
        it("never compresses " .. src, function()
            local native = png[1] / png[2] * (CELL_H / CELL_W)
            assert.is_true(cell_width(png[1], png[2], CELL_W, CELL_H) >= native)
        end)
    end

    it("rounds $T$ up to 2 cells (regression: was 1)", function()
        assert.are.equal(2, cell_width(23, 39, CELL_W, CELL_H))
    end)
end)

describe("fill_box", function()
    -- A 3×2-cell image over a 4-line block, as snacks' render_grid emits it:
    -- two overlay rows (the first concealing the whole range) and a
    -- conceal_lines mark over the two lines below the image. "x" stands in for
    -- one placeholder cell.
    local function block_marks()
        return {
            { row = 4, col = 2, end_row = 7, end_col = 4, conceal = "", virt_text_pos = "overlay",
              virt_text = { { "xxx", "Img" } }, virt_text_hide = false, virt_text_win_col = 2 },
            { row = 5, col = 0, virt_text_pos = "overlay",
              virt_text = { { "xxx", "Img" } }, virt_text_hide = false, virt_text_win_col = 2 },
            { row = 6, end_row = 7, col = 0, conceal_lines = "", virt_text_hide = false },
        }
    end

    it("pads each overlay row to the box width and drops conceal", function()
        local marks = fill_box(block_marks(), 5)
        for i = 1, 2 do
            assert.is_nil(marks[i].conceal)
            assert.are.same({ { "xxx", "Img" }, { "  " } }, marks[i].virt_text)
        end
    end)

    it("sets hl_mode combine", function()
        for _, m in ipairs(fill_box(block_marks(), 5)) do
            assert.are.equal("combine", m.hl_mode)
        end
    end)

    it("replaces conceal_lines with blank overlay rows at the image column", function()
        local marks = fill_box(block_marks(), 5)
        assert.are.equal(4, #marks)
        for i, row in ipairs { 6, 7 } do
            local m = marks[2 + i]
            assert.is_nil(m.conceal_lines)
            assert.are.equal(row, m.row)
            assert.are.equal(0, m.col)
            assert.are.equal("overlay", m.virt_text_pos)
            assert.are.same({ { "     " } }, m.virt_text)
            assert.are.equal(2, m.virt_text_win_col)
            assert.is_false(m.virt_text_hide)
        end
    end)

    it("adds no rows when the image is as tall as the box", function()
        local marks = block_marks()
        marks[3] = nil
        assert.are.equal(2, #fill_box(marks, 5))
    end)

    it("adds no padding when the image is as wide as the box", function()
        assert.are.same({ { "xxx", "Img" } }, fill_box(block_marks(), 3)[1].virt_text)
    end)

    it("returns virt_lines marks unchanged", function()
        local vl = { row = 8, col = 0, virt_lines = { { { "  " }, { "xxx", "Img" } } }, virt_text_hide = false }
        assert.are.same(vl, fill_box({ vl }, 5)[1])
    end)

    it("turns a single inline mark into a padded, unconcealed overlay", function()
        -- render_grid's shape for a one-line source and a one-row image.
        local marks = fill_box({
            { row = 4, col = 0, end_row = 4, end_col = 12, conceal = "", virt_text_pos = "inline",
              virt_text = { { "xxx", "Img" } }, virt_text_hide = true },
        }, 12)
        assert.are.equal(1, #marks)
        assert.are.equal("overlay", marks[1].virt_text_pos)
        assert.is_nil(marks[1].conceal)
        assert.are.equal("combine", marks[1].hl_mode)
        assert.are.same({ { "xxx", "Img" }, { (" "):rep(9) } }, marks[1].virt_text)
        assert.is_nil(marks[1].virt_text_win_col)
    end)

    it("does not modify the input", function()
        local marks = block_marks()
        fill_box(marks, 5)
        assert.are.same(block_marks(), marks)
    end)
end)

describe("footprint_width", function()
    it("measures the widest line right of the start column", function()
        assert.are.equal(9, footprint_width({ "  $$", "  a + bcdef", "  $$" }, 2, 4))
    end)

    it("works on a one-line span", function()
        assert.are.equal(14, footprint_width({ "  $$ E = mc^2 $$" }, 2, 16))
    end)

    it("is nil with text left of the start column on the first line", function()
        assert.is_nil(footprint_width({ "x $$", "  a", "  $$" }, 2, 4))
    end)

    it("is nil with text left of the start column on a later line", function()
        assert.is_nil(footprint_width({ "  $$", "x a", "  $$" }, 2, 4))
    end)

    it("is nil with text after the end column on the last line", function()
        assert.is_nil(footprint_width({ "$$", "a", "$$ x" }, 0, 2))
    end)
end)

describe("with_line_hl", function()
    local function rows()
        return { { { "  " }, { "xxx", "Img" } }, { { "  " }, { "yyy", "Img" } } }
    end

    it("combines the hl into each chunk and gives it to chunks without a group", function()
        local lines = with_line_hl(rows(), "Ln", 0)
        assert.are.same({ "  ", "Ln" }, lines[1][1])
        assert.are.same({ "xxx", { "Ln", "Img" } }, lines[1][2])
        assert.are.same({ "yyy", { "Ln", "Img" } }, lines[2][2])
    end)

    it("appends a pad chunk in the hl to each row", function()
        for _, line in ipairs(with_line_hl(rows(), "Ln", 4)) do
            assert.are.same({ "    ", "Ln" }, line[#line])
            assert.are.equal(3, #line)
        end
    end)

    it("appends an empty chunk for pad 0", function()
        assert.are.same({ "", "Ln" }, with_line_hl(rows(), "Ln", 0)[1][3])
    end)

    it("does not modify the input", function()
        local lines = rows()
        with_line_hl(lines, "Ln", 4)
        assert.are.same(rows(), lines)
    end)
end)
