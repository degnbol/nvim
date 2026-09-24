---@diagnostic disable: undefined-global
-- Only the span scan is tested: headless has no redraw to drive the
-- decoration provider.

local textcolor = require("tex.textcolor")

describe("tex.textcolor.spans", function()
    it("covers the body of a known colour, braces excluded", function()
        assert.are.same(
            { { col = 18, end_col = 21, group = "texTextcolorRed" } },
            textcolor.spans([[x \textcolor{red}{abc} y]]))
    end)

    it("finds several on one line, including the *text variants", function()
        local spans = textcolor.spans([[\textcolor{blue}{a}\textcolor{yellowtext}{bc}]])
        assert.are.same({
            { col = 17, end_col = 18, group = "texTextcolorBlue" },
            { col = 42, end_col = 44, group = "texTextcolorYellowtext" },
        }, spans)
    end)

    it("skips unknown colours and an empty body", function()
        assert.are.same({}, textcolor.spans([[\textcolor{chartreuse}{a} \textcolor{red}{}]]))
    end)
end)
