---@diagnostic disable: undefined-global
local scroll_target = require("utils/hscroll").scroll_target

-- widths index i = buffer line (top + i - 1); need = leftcol + count.
describe("hscroll.scroll_target", function()
    it("returns the cursor line when it qualifies", function()
        local widths = { 5, 100, 5 }
        assert.are.equal(2, scroll_target(widths, 1, 2, 10, 1, 3))
    end)

    it("stays on the cursor line even if a nearer line is shorter", function()
        -- cursor line (2) qualifies; band offers no better move.
        local widths = { 3, 20, 3 }
        assert.are.equal(2, scroll_target(widths, 1, 2, 10, 1, 3))
    end)

    it("picks the nearest qualifying line when the cursor line is short", function()
        -- cursor on line 3 (short); line 4 near+long, line 1 far+long.
        local widths = { 100, 3, 3, 100 }
        assert.are.equal(4, scroll_target(widths, 1, 3, 10, 1, 4))
    end)

    it("breaks distance ties toward the screen centre", function()
        -- cursor 2; qualifying lines 1 and 3 are equidistant (dist 1),
        -- centre = 3, so line 3 (nearer centre) wins.
        local widths = { 100, 3, 100, 3, 3 }
        assert.are.equal(3, scroll_target(widths, 1, 2, 10, 1, 5))
    end)

    it("returns nil when nothing in the band qualifies", function()
        local widths = { 3, 3, 3 }
        assert.is_nil(scroll_target(widths, 1, 2, 10, 1, 3))
    end)

    it("respects the band, ignoring qualifying lines outside it", function()
        -- line 5 is long but outside the band [2,4]; cursor line 3 short.
        local widths = { 3, 3, 3, 3, 100 }
        assert.is_nil(scroll_target(widths, 1, 3, 10, 2, 4))
    end)

    it("uses a strict width > need boundary", function()
        -- width exactly equal to need does not qualify.
        assert.is_nil(scroll_target({ 10 }, 1, 1, 10, 1, 1))
        assert.are.equal(1, scroll_target({ 11 }, 1, 1, 10, 1, 1))
    end)
end)
