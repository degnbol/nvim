---@diagnostic disable: undefined-global
-- PNG dimensions below were measured by compiling each expression through the
-- same pipeline pickers.lua uses (strut-wrapped \begin{math}, tectonic, magick
-- -density 192, no -trim). Cell size 9×18 is snacks' headless fallback
-- (terminal.lua M.size). Regenerate if the strut or font_size changes.
local cell_width = require("utils.inline_math").cell_width
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
