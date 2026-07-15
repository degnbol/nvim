---@diagnostic disable: undefined-global
-- Tests for lua/julia_fndef.lua

local jf = require("julia_fndef")

--- Set up a julia buffer with deterministic indent (4-space, expandtab) and a
--- parsed tree, place the cursor, and return the buffer.
local function make_buf(lines, row, col)
    local buf = vim.api.nvim_create_buf(true, false)
    vim.api.nvim_set_current_buf(buf)
    vim.api.nvim_buf_set_lines(buf, 0, -1, false, lines)
    vim.bo[buf].filetype = "julia"
    vim.bo[buf].shiftwidth = 4
    vim.bo[buf].expandtab = true
    vim.treesitter.get_parser(buf, "julia"):parse()
    vim.api.nvim_win_set_cursor(0, { row or 1, col or 0 })
    return buf
end

local function lines(buf)
    return vim.api.nvim_buf_get_lines(buf, 0, -1, false)
end

describe("julia_fndef.toggle", function()
    it("short → long", function()
        local buf = make_buf({ "f(x) = x^2" })
        jf.toggle()
        assert.are.same({ "function f(x)", "    x^2", "end" }, lines(buf))
    end)

    it("long → short", function()
        local buf = make_buf({ "function f(x)", "    x^2", "end" })
        jf.toggle()
        assert.are.same({ "f(x) = x^2" }, lines(buf))
    end)

    it("keeps a return-type annotation on the signature", function()
        local buf = make_buf({ "f(x)::Int = x" })
        jf.toggle()
        assert.are.same({ "function f(x)::Int", "    x", "end" }, lines(buf))
    end)

    it("keeps a where clause on the signature", function()
        local buf = make_buf({ "f(x) where T = x" })
        jf.toggle()
        assert.are.same({ "function f(x) where T", "    x", "end" }, lines(buf))
    end)

    it("preserves a macro prefix", function()
        local buf = make_buf({ "@inline f(x) = x" }, 1, 8)
        jf.toggle()
        assert.are.same({ "@inline function f(x)", "    x", "end" }, lines(buf))
    end)

    it("preserves a preceding docstring", function()
        local buf = make_buf({ '"""doc"""', "function f(x)", "    x", "end" }, 2, 0)
        jf.toggle()
        assert.are.same({ '"""doc"""', "f(x) = x" }, lines(buf))
    end)

    it("respects the def's own indentation", function()
        local buf = make_buf({ "module M", "    f(x) = x", "end" }, 2, 4)
        jf.toggle()
        assert.are.same(
            { "module M", "    function f(x)", "        x", "    end", "end" },
            lines(buf))
    end)

    it("refuses to inline a multi-statement body", function()
        local src = { "function f(x)", "    a = 1", "    a", "end" }
        local buf = make_buf(src)
        jf.toggle()
        assert.are.same(src, lines(buf))
    end)

    it("ignores a plain variable assignment", function()
        local buf = make_buf({ "a = 1" })
        jf.toggle()
        assert.are.same({ "a = 1" }, lines(buf))
    end)
end)

describe("julia_fndef.conform", function()
    it("conforms all same-name methods to the cursor's long form", function()
        local buf = make_buf({
            "f(x) = x",
            "f(x, y) = x + y",
            "function f(x, y, z)",
            "    x",
            "end",
        }, 3, 0)
        jf.conform()
        assert.are.same({
            "function f(x)",
            "    x",
            "end",
            "function f(x, y)",
            "    x + y",
            "end",
            "function f(x, y, z)",
            "    x",
            "end",
        }, lines(buf))
    end)

    it("conforms to short form, leaving other names untouched", function()
        local buf = make_buf({
            "function f(x)", "    x", "end",
            "g(x) = x",
            "f(y) = y",
        }, 5, 0)
        jf.conform()
        assert.are.same({ "f(x) = x", "g(x) = x", "f(y) = y" }, lines(buf))
    end)

    it("skips multi-statement bodies when conforming to short", function()
        local src = {
            "f(x) = x",
            "function f(y)", "    a = 1", "    a", "end",
        }
        local buf = make_buf(src, 1, 0)
        jf.conform()
        assert.are.same(src, lines(buf))
    end)
end)
