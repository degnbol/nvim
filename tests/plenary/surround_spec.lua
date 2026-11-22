---@diagnostic disable: undefined-global
-- Tests for lua/surroundings.lua through mini.surround with its default
-- mappings: `sa` add, `sd` delete, `sr` replace.

vim.cmd.packadd "mini.nvim"
vim.cmd.packadd "nvim-treesitter-textobjects"
-- Provides the julia injection query that parses docstrings as markdown.
vim.cmd.packadd "nvim-treesitter"
-- ftplugin/markdown.lua requires plugins the harness does not set up.
vim.cmd "filetype plugin off"
require("mini.surround").setup { custom_surroundings = require "surroundings" }

--- Lines of a scratch buffer of `lines` after normal-mode `keys`, with the
--- cursor on the first byte of the first `target`. The buffer is wiped after.
--- @param lines string[]
--- @param filetype string parsed with its tree-sitter parser, injections included, when one exists
--- @param target string
--- @param keys string e.g. "sdF"
--- @return string[]
local function lines_after(lines, filetype, target, keys)
    local buf = vim.api.nvim_create_buf(false, true)
    vim.api.nvim_set_current_buf(buf)
    vim.api.nvim_buf_set_lines(buf, 0, -1, false, lines)
    vim.bo[buf].filetype = filetype
    local parser = vim.treesitter.get_parser(buf, nil, { error = false })
    if parser then parser:parse(true) end
    local row = assert(vim.iter(ipairs(lines)):find(function(_, l) return l:find(target, 1, true) end))
    vim.api.nvim_win_set_cursor(0, { row, lines[row]:find(target, 1, true) - 1 })
    vim.api.nvim_feedkeys(keys, "mx", false)
    local result = vim.api.nvim_buf_get_lines(buf, 0, -1, true)
    vim.api.nvim_buf_delete(buf, { force = true })
    return result
end

describe("t triple quotes", function()
    it("deletes '''", function()
        assert.are.same({ "x = x" }, lines_after({ "x = '''x'''" }, "python", "x'''", "sdt"))
    end)
    it('deletes """', function()
        assert.are.same({ "x = x" }, lines_after({ 'x = """x"""' }, "python", 'x"""', "sdt"))
    end)
end)

describe("F code fence", function()
    it("deletes both fence lines", function()
        local lines = { "a", "", "```py", "x = 1", "y = 2", "```", "", "b" }
        assert.are.same({ "a", "", "x = 1", "y = 2", "", "b" }, lines_after(lines, "markdown", "x =", "sdF"))
    end)
    it("deletes a fence that ends on the last line", function()
        assert.are.same({ "a", "x" }, lines_after({ "a", "```", "x", "```" }, "markdown", "x", "sdF"))
    end)
    it("deletes indented fence lines in a list item", function()
        local lines = { "- a", "  ```py", "  x", "  ```", "- b" }
        assert.are.same({ "- a", "  x", "- b" }, lines_after(lines, "markdown", "x", "sdF"))
    end)
    it("deletes indented fence lines in a blockquote", function()
        local lines = { "> ```", "> x", "> ```", "" }
        assert.are.same({ "> x", "" }, lines_after(lines, "markdown", "x", "sdF"))
    end)
    it("leaves a fence without a closing delimiter", function()
        local lines = { "a", "```", "x", "y" }
        assert.are.same(lines, lines_after(lines, "markdown", "x", "sdF"))
    end)
    it("deletes a ~~~ fence", function()
        assert.are.same({ "x", "" }, lines_after({ "~~~", "x", "~~~", "" }, "markdown", "x", "sdF"))
    end)
    it("deletes the 4-backtick fence around 3-backtick lines", function()
        local lines = { "````md", "```", "x", "```", "````", "" }
        assert.are.same({ "```", "x", "```", "" }, lines_after(lines, "markdown", "x", "sdF"))
    end)
    it("deletes a fence in a julia docstring", function()
        local lines = { '"""', "Doc.", "", "```julia", "f(1)", "```", '"""', "f(x) = x" }
        local expected = { '"""', "Doc.", "", "f(1)", '"""', "f(x) = x" }
        assert.are.same(expected, lines_after(lines, "julia", "f(1)", "sdF"))
    end)
    it("replaces the fence lines, keeping the next line", function()
        local lines = { "```py", "x", "```", "b" }
        assert.are.same({ '"x"', "b" }, lines_after(lines, "markdown", "x", 'srF"'))
    end)
    it("adds fence lines around a paragraph", function()
        local lines = { "a", "b", "", "c" }
        assert.are.same({ "```", "a", "b", "```", "", "c" }, lines_after(lines, "markdown", "a", "saipF"))
    end)
end)
