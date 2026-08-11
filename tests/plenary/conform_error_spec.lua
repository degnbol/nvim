---@diagnostic disable: undefined-global
local parse = require("conform_error").parse

describe("conform_error.parse", function()
    it("extracts line and column from a shfmt error", function()
        local msg = "Formatter 'shfmt' error: /Users/x/step_to_scheme.tsv.sh:10:13: `>` must be followed by a word"
        local line, col, text = parse(msg)
        assert.are.equal(10, line)
        assert.are.equal(13, col)
        assert.are.equal("`>` must be followed by a word", text)
    end)

    it("is agnostic to the filename token", function()
        local line, col = parse("Formatter 'x' error: <standard input>:5:2: oops")
        assert.are.equal(5, line)
        assert.are.equal(2, col)
    end)

    it("defaults column to 1 when only a line is given", function()
        local line, col, text = parse("tool: file.py:7: bad thing")
        assert.are.equal(7, line)
        assert.are.equal(1, col)
        assert.are.equal("bad thing", text)
    end)

    it("trims trailing whitespace from the message", function()
        local _, _, text = parse("f:3:1: trailing  \n")
        assert.are.equal("trailing", text)
    end)

    it("returns nil when there is no position", function()
        assert.is_nil(parse("Formatter 'shfmt' timeout"))
    end)
end)
