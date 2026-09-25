---@diagnostic disable: undefined-global
local latex = require "utils/latex"

describe("is_inline", function()
    for _, src in ipairs { "$x$", "$`x`$", "\\(x\\)", " $x$ " } do
        it("is true for " .. src, function()
            assert.is_true(latex.is_inline(src))
        end)
    end

    for _, src in ipairs { "$$x$$", "\\[x\\]" } do
        it("is false for " .. src, function()
            assert.is_false(latex.is_inline(src))
        end)
    end
end)

describe("math_body", function()
    it("strips $", function()
        assert.are.equal("x", latex.math_body("$x$"))
    end)

    it("keeps whitespace inside $$", function()
        assert.are.equal(" x ", latex.math_body("$$ x $$"))
    end)

    it("strips \\[ \\]", function()
        assert.are.equal("x", latex.math_body("\\[x\\]"))
    end)

    it("keeps the newline before an environment", function()
        local body = latex.math_body("$$\n\\begin{aligned}x\\end{aligned}\n$$")
        assert.is_truthy(body:find("^\n\\begin{aligned}"))
    end)

    it("starts with \\begin when the environment follows the delimiter", function()
        local body = latex.math_body("$$\\begin{align}x\\end{align}$$")
        assert.is_truthy(body:find("^\\begin"))
    end)

    it("strips a backtick inside $", function()
        assert.are.equal("x", latex.math_body("$`x`$"))
    end)
end)
