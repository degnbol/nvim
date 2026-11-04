---@diagnostic disable: undefined-global
-- Tests for lua/utils/markdown_table.lua

local md = require "utils.markdown_table"

-- ftplugin/markdown.lua requires mini.hipatterns, which the harness lacks.
vim.cmd("filetype plugin off")

local buf

--- @param lines string[]
local function set(lines)
    vim.api.nvim_buf_set_lines(buf, 0, -1, false, lines)
end

--- @return string[]
local function get()
    return vim.api.nvim_buf_get_lines(buf, 0, -1, false)
end

--- Align rows [r1, r2] (default: whole buffer) and return the buffer lines.
--- @param opts utils.markdown_table.Opts|nil
--- @param r1 integer|nil
--- @param r2 integer|nil
--- @return string[]
local function align(opts, r1, r2)
    md.align(buf, r1 or 0, r2 or vim.api.nvim_buf_line_count(buf) - 1, opts)
    return get()
end

describe("markdown_table.align", function()
    before_each(function()
        buf = vim.api.nvim_create_buf(true, false)
        vim.bo[buf].filetype = "markdown"
    end)

    after_each(function()
        pcall(vim.api.nvim_buf_delete, buf, { force = true })
    end)

    it("aligns at raw width", function()
        set({ "| a | b |", "|---|---|", "| **x** | yy |" })
        assert.are.same({
            "| a     | b   |",
            "| ----- | --- |",
            "| **x** | yy  |",
        }, align())
    end)

    it("lets a cell over max_width overflow while other rows stay aligned", function()
        set({
            "| id | note |",
            "| --- | --- |",
            "| 1 | short |",
            "| 2 | much longer than five |",
        })
        assert.are.same({
            "| id  | note  |",
            "| --- | ----- |",
            "| 1   | short |",
            "| 2   | much longer than five |",
        }, align({ max_width = 5 }))
    end)

    it("keeps delimiter colons and pads all cells left by default", function()
        set({ "| L | C | R |", "|:--|:-:|--:|", "| a | b | cc |" })
        assert.are.same({
            "| L   | C   | R   |",
            "| :-- | :-: | --: |",
            "| a   | b   | cc  |",
        }, align())
    end)

    it("pads right and centre columns with align_cells", function()
        set({ "| L | C | R |", "|:--|:-:|--:|", "| a | b | cc |", "| a | bb | c |" })
        assert.are.same({
            "| L   |  C  |   R |",
            "| :-- | :-: | --: |",
            "| a   |  b  |  cc |",
            "| a   | bb  |   c |",
        }, align({ align_cells = true }))
    end)

    it("keeps an escaped pipe inside its cell", function()
        set({ "| a | b |", "|---|---|", "| x \\| y | z |" })
        assert.are.same({
            "| a      | b   |",
            "| ------ | --- |",
            "| x \\| y | z   |",
        }, align())
    end)

    it("leaves pipes in a code fence untouched and formats a markdown fence", function()
        local lines = {
            "```",
            "|a|b|",
            "|-|-|",
            "```",
            "",
            "```markdown",
            "|a|b|",
            "|---|---|",
            "```",
        }
        set(lines)
        local out = align()
        assert.are.same(vim.list_slice(lines, 1, 6), vim.list_slice(out, 1, 6))
        assert.are.same({ "| a   | b   |", "| --- | --- |", "```" }, vim.list_slice(out, 7, 9))
    end)

    it("adds missing leading and trailing pipes", function()
        set({ "a | b", "--|--", "1 | 2" })
        assert.are.same({
            "| a   | b   |",
            "| --- | --- |",
            "| 1   | 2   |",
        }, align())
    end)

    it("keeps the prefix of a table in a blockquote and in a list item", function()
        set({ "> |a|b|", "> |---|---|", "", "- item", "", "  |a|b|", "  |---|---|" })
        assert.are.same({
            "> | a   | b   |",
            "> | --- | --- |",
            "",
            "- item",
            "",
            "  | a   | b   |",
            "  | --- | --- |",
        }, align())
    end)

    it("fills missing cells", function()
        set({ "| a | b | c |", "|---|---|---|", "| 1 |", "| 1 | 2 | 3 | 4 |" })
        assert.are.same({
            "| a   | b   | c   |     |",
            "| --- | --- | --- | --- |",
            "| 1   |     |     |     |",
            "| 1   | 2   | 3   | 4   |",
        }, align())
    end)

    it("keeps cells after empty cells in their columns", function()
        set({ "| a | b | c | d |", "|---|---|---|---|", "| a || c ||", "| a ||| d |" })
        assert.are.same({
            "| a   | b   | c   | d   |",
            "| --- | --- | --- | --- |",
            "| a   |     | c   |     |",
            "| a   |     |     | d   |",
        }, align())
    end)

    it("aligns multibyte and wide characters by display width", function()
        set({ "| é | 漢字 |", "|---|---|", "| abcd | x |" })
        assert.are.same({
            "| é    | 漢字 |",
            "| ---- | ---- |",
            "| abcd | x    |",
        }, align())
    end)

    it("adds pipes to a blockquote row without them", function()
        set({ "> |a|b|", "> |---|---|", "> c | d" })
        assert.are.same({
            "> | a   | b   |",
            "> | --- | --- |",
            "> | c   | d   |",
        }, align())
    end)

    it("returns 0 in a buffer without a markdown parser", function()
        vim.bo[buf].filetype = ""
        set({ "|a|b|", "|-|-|" })
        assert.are.equal(0, md.align(buf, 0, 1))
    end)

    it("skips a table with a row starting with `||` and warns", function()
        local lines = { "| a | b |", "|---|---|", "|| c |" }
        set(lines)
        local notify = vim.notify
        local level
        vim.notify = function(_, l) level = l end
        local out = align()
        vim.notify = notify
        assert.are.equal(vim.log.levels.WARN, level)
        assert.are.same(lines, out)
    end)

    it("formats only, and all of, the tables a range intersects", function()
        set({ "|a|b|", "|-|-|", "", "|c|d|", "|-|-|", "|e|f|" })
        assert.are.same({
            "|a|b|",
            "|-|-|",
            "",
            "| c   | d   |",
            "| --- | --- |",
            "| e   | f   |",
        }, align(nil, 4, 4))
    end)

    it("returns 0 outside a table and on the blank line below it", function()
        set({ "text", "|a|b|", "|-|-|", "", "more" })
        assert.are.equal(0, md.align(buf, 0, 0))
        assert.are.equal(0, md.align(buf, 3, 3))
        assert.are.equal(1, md.align(buf, 2, 2))
    end)

    it("keeps an extmark before the first changed byte", function()
        set({ "| aaa | b |", "| --- | --- |" })
        local ns = vim.api.nvim_create_namespace("markdown_table_spec")
        local id = vim.api.nvim_buf_set_extmark(buf, ns, 0, 3, {})
        align()
        assert.are.same({ 0, 3 }, vim.api.nvim_buf_get_extmark_by_id(buf, ns, id, {}))
    end)

    it("is idempotent", function()
        set({ "> | L | C |", "> |:--|:-:|", "> | longer | x |" })
        local once = align({ align_cells = true })
        assert.are.same(once, align({ align_cells = true }))
    end)
end)

describe("markdown table keymap", function()
    vim.cmd.runtime("lua/autocmds/markdown_table.lua")

    before_each(function()
        buf = vim.api.nvim_create_buf(true, false)
        vim.api.nvim_set_current_buf(buf)
        vim.bo[buf].filetype = "markdown"
    end)

    after_each(function()
        pcall(vim.api.nvim_buf_delete, buf, { force = true })
    end)

    it("uses a count as max_width in visual mode", function()
        set({ "|a|b|", "|-|-|", "|long cell|x|" })
        vim.api.nvim_feedkeys(vim.keycode("ggVG4<localleader>a"), "xt", false)
        assert.are.same({
            "| a    | b   |",
            "| ---- | --- |",
            "| long cell | x   |",
        }, get())
    end)
end)
