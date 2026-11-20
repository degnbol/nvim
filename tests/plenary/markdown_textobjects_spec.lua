---@diagnostic disable: undefined-global
-- Tests for queries/markdown_inline/textobjects.scm: @emphasis.inner/outer
-- selected through mini.ai with the `*` spec ftplugin/markdown.lua sets, in a
-- markdown buffer, where markdown_inline is an injected language.

vim.cmd.packadd "mini.nvim"
vim.cmd.packadd "nvim-treesitter-textobjects"
vim.cmd "filetype plugin off"
require("mini.ai").setup()

local LINE = "a *b c* **d e** ***f*** _g_ *i **j** k* z"

--- Yank emphasis textobject `ai_type` with the cursor on the first occurrence
--- of `target`.
--- @param target string text in LINE to place the cursor on
--- @param ai_type string "i" or "a", with an optional count before it, e.g. "2a"
--- @return string yanked
local function yank(target, ai_type)
    local buf = vim.api.nvim_create_buf(false, false)
    vim.api.nvim_set_current_buf(buf)
    vim.api.nvim_buf_set_lines(buf, 0, -1, false, { LINE })
    vim.bo[buf].filetype = "markdown"
    vim.cmd.source "ftplugin/markdown.lua"
    vim.treesitter.get_parser(buf):parse(true)
    vim.api.nvim_win_set_cursor(0, { 1, assert(LINE:find(target, 1, true)) - 1 })
    vim.fn.setreg('"', "")
    vim.cmd.normal("y" .. ai_type .. "*")
    vim.api.nvim_buf_delete(buf, { force = true })
    return vim.fn.getreg('"')
end

describe("markdown emphasis textobject", function()
    it("inner of *emphasis* excludes delimiters", function()
        assert.are.equal("b c", yank("b c", "i"))
    end)
    it("inner of **strong** excludes both delimiter pairs", function()
        assert.are.equal("d e", yank("d e", "i"))
    end)
    it("inner of ***both*** is the text", function()
        assert.are.equal("f", yank("f*", "i"))
    end)
    it("inner of _underscore_ emphasis", function()
        assert.are.equal("g", yank("g_", "i"))
    end)
    it("inner picks the innermost span around the cursor", function()
        assert.are.equal("j", yank("j", "i"))
        assert.are.equal("i **j** k", yank("i *", "i"))
    end)
    it("outer includes delimiters on both sides", function()
        assert.are.equal("*b c*", yank("b c", "a"))
        assert.are.equal("**d e**", yank("d e", "a"))
    end)
    it("count selects the enclosing span", function()
        assert.are.equal("*i **j** k*", yank("j", "2a"))
    end)
end)
