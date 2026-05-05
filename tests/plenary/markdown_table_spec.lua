---@diagnostic disable: undefined-global
-- Tests for lua/utils/markdown_table.lua

local md = require("utils.markdown_table")

describe("markdown_table.is_table_line", function()
    it("recognises lines starting with `|`", function()
        assert.is_true(md.is_table_line("| a | b |"))
        assert.is_true(md.is_table_line("  | a | b |"))
        assert.is_true(md.is_table_line("|---|---|"))
    end)
    it("rejects non-table lines", function()
        assert.is_false(md.is_table_line(""))
        assert.is_false(md.is_table_line("text"))
        assert.is_false(md.is_table_line("- bullet"))
    end)
end)

describe("markdown_table.format_block", function()
    it("pads columns to the widest cell", function()
        local out = md.format_block({
            "| Name | Age |",
            "| --- | --- |",
            "| Alice | 30 |",
        })
        -- All rows must have identical visual width when measured by raw byte
        -- length (no concealed delimiters in this fixture).
        assert.are.equal(#out[1], #out[2])
        assert.are.equal(#out[1], #out[3])
        -- Header padded so "Name" cell width matches "Alice".
        assert.is_truthy(out[1]:find("| Name  | Age |"))
    end)

    it("preserves alignment markers", function()
        local out = md.format_block({
            "| L | C | R |",
            "|:---|:---:|---:|",
            "| a | b | c |",
        })
        -- Separator must keep the colon markers and only widen the dashes.
        local sep = out[2]
        assert.is_truthy(sep:find(":%-+ |"), "left marker preserved: " .. sep)
        assert.is_truthy(sep:find(":%-+:"), "centre marker preserved: " .. sep)
        assert.is_truthy(sep:find("%-+:"), "right marker preserved: " .. sep)
    end)

    it("is idempotent", function()
        local input = {
            "| a | b |",
            "| --- | --- |",
            "| longer | x |",
        }
        local once = md.format_block(input)
        local twice = md.format_block(once)
        assert.are.same(once, twice)
    end)

    it("strips concealed delimiters when measuring width", function()
        -- Each col-1 cell has visual width 1 after concealing the
        -- delimiters, so all three rows pad to the same column width.
        -- Raw byte counts differ by the delimiter overhead per cell.
        local out = md.format_block({
            "| H | H |",
            "| --- | --- |",
            "| `x` | a |",
            "| **x** | a |",
            "| ~~x~~ | a |",
        })
        local function len(i) return #out[i] end
        -- `**x**` adds 2 bytes vs `` `x` `` (4 stars vs 2 backticks); same
        -- for `~~x~~` (4 tildes vs 2 backticks).
        assert.are.equal(len(3) + 2, len(4))
        assert.are.equal(len(3) + 2, len(5))
    end)

    it("caps column width at max_width without truncating cells", function()
        local out = md.format_block({
            "| ID | Description |",
            "| --- | --- |",
            "| 1 | short |",
            "| 2 | this is significantly longer than five chars |",
        }, { max_width = 5 })
        -- Column 2 capped to 5: short row is padded to 5, long row keeps full
        -- text and overflows the column (no truncation).
        local short_row = out[3]
        local long_row = out[4]
        assert.is_truthy(short_row:find("| short |"))
        assert.is_truthy(long_row:find("longer than five chars"),
            "long cell content preserved")
    end)

    it("handles escaped pipes inside cells", function()
        local out = md.format_block({
            "| Cmd | Note |",
            "| --- | --- |",
            "| `a \\| b` | piped |",
        })
        -- Escaped pipe should stay literal in the cell, not split it.
        assert.is_truthy(out[3]:find("a \\| b"))
        -- Both rows should have the same number of cells (3 outer pipes plus
        -- one escaped pipe in row 3, so 4 raw `|` total vs 3 in the header).
        local function unescaped_pipes(s)
            local n = 0
            local i = 1
            while i <= #s do
                if s:sub(i, i) == "\\" and s:sub(i + 1, i + 1) == "|" then
                    i = i + 2
                elseif s:sub(i, i) == "|" then
                    n = n + 1; i = i + 1
                else
                    i = i + 1
                end
            end
            return n
        end
        assert.are.equal(unescaped_pipes(out[1]), unescaped_pipes(out[3]))
    end)
end)

describe("markdown_table.find_block_range", function()
    local buf

    before_each(function()
        buf = vim.api.nvim_create_buf(true, false)
    end)

    after_each(function()
        pcall(vim.api.nvim_buf_delete, buf, { force = true })
    end)

    it("finds a contiguous range from a row inside the table", function()
        vim.api.nvim_buf_set_lines(buf, 0, -1, false, {
            "prose",
            "",
            "| a | b |",
            "| --- | --- |",
            "| 1 | 2 |",
            "",
            "more prose",
        })
        local s, e = md.find_block_range(buf, 3) -- separator row
        assert.are.equal(2, s)
        assert.are.equal(4, e)
    end)

    it("returns nil when row is not on a table line", function()
        vim.api.nvim_buf_set_lines(buf, 0, -1, false, { "prose only" })
        local s, e = md.find_block_range(buf, 0)
        assert.is_nil(s)
        assert.is_nil(e)
    end)

    it("handles a table at the start of the buffer", function()
        vim.api.nvim_buf_set_lines(buf, 0, -1, false, {
            "| a | b |",
            "| --- | --- |",
            "| 1 | 2 |",
            "after",
        })
        local s, e = md.find_block_range(buf, 0)
        assert.are.equal(0, s)
        assert.are.equal(2, e)
    end)
end)

describe("markdown_table.align_at", function()
    local buf

    before_each(function()
        buf = vim.api.nvim_create_buf(true, false)
        vim.api.nvim_set_current_buf(buf)
    end)

    after_each(function()
        pcall(vim.api.nvim_buf_delete, buf, { force = true })
    end)

    it("returns false when not on a table", function()
        vim.api.nvim_buf_set_lines(buf, 0, -1, false, { "prose" })
        assert.is_false(md.align_at(buf, 0))
    end)

    it("aligns the block at the cursor and returns true", function()
        vim.api.nvim_buf_set_lines(buf, 0, -1, false, {
            "preamble",
            "|a|bbbb|",
            "|---|---|",
            "|longer|x|",
        })
        assert.is_true(md.align_at(buf, 1))
        local lines = vim.api.nvim_buf_get_lines(buf, 0, -1, false)
        -- Preamble untouched.
        assert.are.equal("preamble", lines[1])
        -- All table rows now have the same raw width.
        assert.are.equal(#lines[2], #lines[3])
        assert.are.equal(#lines[2], #lines[4])
    end)
end)
