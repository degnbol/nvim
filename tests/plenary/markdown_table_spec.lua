---@diagnostic disable: undefined-global
-- Tests for lua/utils/markdown_table.lua

local md = require "utils.markdown_table"

-- ftplugin/markdown.lua requires mini.hipatterns, which the harness lacks.
vim.cmd("filetype plugin off")
-- Registers the #head! directive of the markdown_inline conceal overrides.
vim.cmd.runtime("plugin/treesitter.lua")

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

    it("puts one space next to a pipe where the input has whitespace", function()
        set({ "| a | b |", "|---|---|", "| 1 | 22 |" })
        assert.are.same({
            "| a | b  |",
            "|---|----|",
            "| 1 | 22 |",
        }, align())
    end)

    it("keeps a table without spaces next to its pipes tight", function()
        set({ "|a|bb|", "|-|-|", "|ccc|d|" })
        assert.are.same({
            "|a  |bb |",
            "|---|---|",
            "|ccc|d  |",
        }, align())
    end)

    it("uses hidden markup at a cell edge as the gap next to the pipe", function()
        set({ "| a |", "|---|", "| `x` |" })
        assert.are.same({ "| a |", "|---|", "|`x`|" }, align())
    end)

    it("starts every body in a column with a `**bold**` cell 2 columns from the pipe", function()
        set({ "| a | b |", "|---|---|", "| **x** | yy |" })
        assert.are.same({
            "|  a  | b  |",
            "|-----|----|",
            "|**x**| yy |",
        }, align())
    end)

    it("lines up a header body with a code span body", function()
        set({ "| file |", "|---|", "| `name` |" })
        assert.are.same({ "| file |", "|------|", "|`name`|" }, align())
    end)

    it("keeps the space before a backslash escape at a cell edge", function()
        set({ "| a |", "|---|", "| \\* |" })
        assert.are.same({ "| a  |", "|----|", "| \\* |" }, align())
    end)

    it("lets a cell over max_width overflow without spaces while other rows stay aligned", function()
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
            "| 2   |much longer than five|",
        }, align({ max_width = 5 }))
    end)

    it("drops the spaces of a cell over max_width so it fills its slot", function()
        set({ "| a | b | c |", "| --- | --- | --- |", "| 1234567 | 1 | x |", "| 1 | 2 | y |" })
        assert.are.same({
            "| a     | b   | c   |",
            "| ----- | --- | --- |",
            "|1234567| 1   | x   |",
            "| 1     | 2   | y   |",
        }, align({ max_width = 5 }))
    end)

    it("drops only the trailing space of a cell over max_width one short of its slot", function()
        set({ "| a | b |", "| --- | --- |", "| 123456 | 1 |", "| 1 | 2 |" })
        assert.are.same({
            "| a     | b   |",
            "| ----- | --- |",
            "| 123456| 1   |",
            "| 1     | 2   |",
        }, align({ max_width = 5 }))
    end)

    it("shrinks the cells after a cell over max_width until the row is back in its columns", function()
        set({ "| a | b | c |", "| --- | --- | --- |", "| 12345678901 | 1 | x |", "| 1 | 2 | y |" })
        assert.are.same({
            "| a     | b   | c   |",
            "| ----- | --- | --- |",
            "|12345678901|1| x   |",
            "| 1     | 2   | y   |",
        }, align({ max_width = 5 }))
    end)

    it("keeps the end of a right-aligned cell after a cell over max_width on its column end", function()
        set({ "| a | b | c |", "|---|--:|---|", "| 123456789 | 1 | x |", "| 1 | 22222 | y |" })
        assert.are.same({
            "| a     |     b | c |",
            "|-------|------:|---|",
            "|123456789|   1 | x |",
            "| 1     | 22222 | y |",
        }, align({ max_width = 5, align_cells = true }))
    end)

    it("keeps a centred cell after a cell over max_width on its column centre where it fits", function()
        set({ "| a | b | c |", "|---|:-:|---|", "| 123456789 | 1 | x |", "| 1 | 22222 | y |" })
        assert.are.same({
            "| a     |   b   | c |",
            "|-------|:-----:|---|",
            "|123456789| 1   | x |",
            "| 1     | 22222 | y |",
        }, align({ max_width = 5, align_cells = true }))
    end)

    it("keeps delimiter colons and pads all cells left by default", function()
        set({ "| L | C | R |", "|:--|:-:|--:|", "| a | b | cc |" })
        assert.are.same({
            "| L | C | R  |",
            "|:--|:-:|---:|",
            "| a | b | cc |",
        }, align())
    end)

    it("pads right and centre columns with align_cells", function()
        set({ "| L | C | R |", "|:--|:-:|--:|", "| a | b | cc |", "| a | bb | c |" })
        assert.are.same({
            "| L | C  |  R |",
            "|:--|:--:|---:|",
            "| a | b  | cc |",
            "| a | bb |  c |",
        }, align({ align_cells = true }))
    end)

    it("keeps an escaped pipe inside its cell", function()
        set({ "| a | b |", "|---|---|", "| x \\| y | z |" })
        assert.are.same({
            "| a      | b |",
            "|--------|---|",
            "| x \\| y | z |",
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
        assert.are.same({ "|a  |b  |", "|---|---|", "```" }, vim.list_slice(out, 7, 9))
    end)

    it("adds missing leading and trailing pipes", function()
        set({ "a | b", "--|--", "1 | 2" })
        assert.are.same({
            "|a  | b |",
            "|---|---|",
            "|1  | 2 |",
        }, align())
    end)

    it("keeps the prefix of a table in a blockquote and in a list item", function()
        set({ "> |a|b|", "> |---|---|", "", "- item", "", "  |a|b|", "  |---|---|" })
        assert.are.same({
            "> |a  |b  |",
            "> |---|---|",
            "",
            "- item",
            "",
            "  |a  |b  |",
            "  |---|---|",
        }, align())
    end)

    it("fills missing cells", function()
        set({ "| a | b | c |", "|---|---|---|", "| 1 |", "| 1 | 2 | 3 | 4 |" })
        assert.are.same({
            "| a | b | c |   |",
            "|---|---|---|---|",
            "| 1 |   |   |   |",
            "| 1 | 2 | 3 | 4 |",
        }, align())
    end)

    it("keeps cells after empty cells in their columns", function()
        set({ "| a | b | c | d |", "|---|---|---|---|", "| a || c ||", "| a ||| d |" })
        assert.are.same({
            "| a | b | c | d |",
            "|---|---|---|---|",
            "| a |   | c |   |",
            "| a |   |   | d |",
        }, align())
    end)

    it("aligns multibyte and wide characters by display width", function()
        set({ "| é | 漢字 |", "|---|---|", "| abcd | x |" })
        assert.are.same({
            "| é    | 漢字 |",
            "|------|------|",
            "| abcd | x    |",
        }, align())
    end)

    it("adds pipes to a blockquote row without them", function()
        set({ "> |a|b|", "> |---|---|", "> c | d" })
        assert.are.same({
            "> |a  | b |",
            "> |---|---|",
            "> |c  | d |",
        }, align())
    end)

    it("returns 0 in a buffer without a markdown parser", function()
        -- The markdown FileType has already started a parser on `buf`.
        vim.api.nvim_buf_delete(buf, { force = true })
        buf = vim.api.nvim_create_buf(true, false)
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
            "|c  |d  |",
            "|---|---|",
            "|e  |f  |",
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

    it("is idempotent with cells over max_width", function()
        set({
            "| a | b | c |",
            "| --- | --- | --- |",
            "| 1234567 | 1 | x |",
            "| 123456 | 1 | x |",
            "| 12345678901 | 1 | x |",
            "| 1 | 2 | y |",
        })
        local once = align({ max_width = 5 })
        assert.are.same(once, align({ max_width = 5 }))
    end)

    it("is idempotent with markup at cell edges", function()
        set({ "| L | C | R |", "|:--|:-:|--:|", "| **a** | `b` | *c* |", "| plain | text | here |" })
        local once = align({ align_cells = true })
        assert.are.same(once, align({ align_cells = true }))
    end)
end)

describe("markdown_table.layout and pads", function()
    local Conceal = require "agentic.utils.conceal"
    local OPTS = { max_width = 30, align_cells = true }

    before_each(function()
        buf = vim.api.nvim_create_buf(true, false)
        vim.bo[buf].filetype = "markdown"
    end)

    after_each(function()
        pcall(vim.api.nvim_buf_delete, buf, { force = true })
    end)

    --- The table on the first row of `lines`, its conceals and layout.
    --- @param lines string[]
    --- @return utils.markdown_table.Table
    --- @return table<integer, agentic.utils.Conceal.Range[]>
    --- @return utils.markdown_table.Layout
    local function measure(lines)
        set(lines)
        vim.treesitter.start(buf, "markdown")
        local tbl = md.parse(buf, md.tables(buf, 0, 0)[1])
        local conceals = md.conceals(buf, tbl)
        return tbl, conceals, md.layout(tbl, conceals, OPTS)
    end

    --- Screen vcols per row with the pads drawn: of each cell's right pipe,
    --- of its first visible character, and after its last one. The delimiter
    --- row has only pipes.
    --- @param tbl utils.markdown_table.Table
    --- @param conceals table<integer, agentic.utils.Conceal.Range[]>
    --- @param layout utils.markdown_table.Layout
    --- @param level integer
    --- @param unconcealed_row integer|nil
    --- @return { pipes: integer[], starts: integer[], ends: integer[] }[]
    local function screen(tbl, conceals, layout, level, unconcealed_row)
        local pads = md.pads(tbl, layout, conceals, level, unconcealed_row)
        local out = {}
        for i, row in ipairs(tbl.rows) do
            local inlines, pad_at = {}, {}
            for _, pad in ipairs(pads) do
                if pad.row == row.row then
                    inlines[#inlines + 1] = { col = pad.col, width = pad.width }
                    pad_at[pad.col] = (pad_at[pad.col] or 0) + pad.width
                end
            end
            local row_level = row.row == unconcealed_row and 0 or level
            local shown = Conceal.columns(row.line, conceals[i], inlines, 8, row_level).shown
            local pos = { pipes = {}, starts = {}, ends = {} }
            for c, cell in ipairs(row.cells) do
                local cl = layout.cells[i][c]
                pos.pipes[c] = shown[cell.to]
                if not row.delimiter then
                    pos.starts[c] = shown[cl.vis_col]
                    pos.ends[c] = shown[cl.vis_end] - (pad_at[cl.vis_end] or 0)
                end
            end
            out[i] = pos
        end
        return out
    end

    --- Assert every row has the pipes of the first row, and every row but
    --- the delimiter the same `starts` in `left` columns and `ends` in
    --- `right` columns.
    local function assert_aligned(rows, left, right)
        for i = 2, #rows do
            assert.are.same(rows[1].pipes, rows[i].pipes, "pipes of row " .. i)
            if next(rows[i].starts) then
                for _, c in ipairs(left or {}) do
                    assert.are.equal(rows[1].starts[c], rows[i].starts[c], ("start of row %d col %d"):format(i, c))
                end
                for _, c in ipairs(right or {}) do
                    assert.are.equal(rows[1].ends[c], rows[i].ends[c], ("end of row %d col %d"):format(i, c))
                end
            end
        end
    end

    local EXAMPLE = {
        "| correlation with K50 | trypsin    | chymotrypsin |",
        "| -------------------- | ---------- | ------------ |",
        "| K/R count            | **-0.215** | -0.064       |",
        "| F/Y/W/L count        | -0.168     | **-0.414**   |",
    }

    it("lines up a formatted table at every conceal level, at the same columns", function()
        local tbl, conceals, layout = measure(EXAMPLE)
        assert.are_not.same(md.pads(tbl, layout, conceals, 0), md.pads(tbl, layout, conceals, 3))
        assert.are.same({ 23, 38, 55 }, screen(tbl, conceals, layout, 0)[1].pipes)
        for level = 0, 3 do
            local rows = screen(tbl, conceals, layout, level)
            assert_aligned(rows, { 1, 2, 3 })
            assert.are.same({ 23, 38, 55 }, rows[1].pipes, "level " .. level)
            assert.are.same({ 2, 27, 42 }, rows[1].starts, "level " .. level)
        end
    end)

    it("keeps the columns on an unconcealed row", function()
        local tbl, conceals, layout = measure(EXAMPLE)
        local rows = screen(tbl, conceals, layout, 3, tbl.rows[3].row)
        assert_aligned(rows, { 1, 2, 3 })
        assert.are.same({ 23, 38, 55 }, rows[1].pipes)
    end)

    it("lines up an unformatted table", function()
        local tbl, conceals, layout = measure({
            "|a|**bold**|",
            "|-|-|",
            "|   long text | x |",
        })
        for level = 0, 3 do
            assert_aligned(screen(tbl, conceals, layout, level), { 1, 2 })
        end
    end)

    it("fills the hidden width of mid-cell markup", function()
        local tbl, conceals, layout = measure({
            "| a | b |",
            "| - | - |",
            "| x **y** z | 1 |",
            "| plain | 2 |",
        })
        for level = 0, 3 do
            assert_aligned(screen(tbl, conceals, layout, level), { 1, 2 })
        end
    end)

    it("lines up a long link over the cap once its markup is concealed", function()
        local tbl, conceals, layout = measure({
            "| name | b |",
            "| - | - |",
            "| [docs](https://example.com/a/very/long/path) | 1 |",
            "| plain | 2 |",
        })
        assert.is_false(layout.cells[3][1].slotted)
        for level = 2, 3 do
            assert_aligned(screen(tbl, conceals, layout, level), { 1, 2 })
        end
        local raw = screen(tbl, conceals, layout, 0)
        assert.are.same(raw[1].pipes, raw[4].pipes)
        assert.is_true(raw[3].pipes[1] > raw[1].pipes[1])
    end)

    it("puts a row back in its columns after a cell over the cap", function()
        local tbl, conceals, layout = measure({
            "| name | b | c |",
            "| - | - | - |",
            "| " .. ("a"):rep(34) .. " | 1 | 2 |",
            "| plain | " .. ("b"):rep(30) .. " | 3 |",
        })
        for _, level in ipairs({ 0, 3 }) do
            local rows = screen(tbl, conceals, layout, level)
            assert.is_true(rows[3].pipes[1] > rows[1].pipes[1])
            assert.are.same(vim.list_slice(rows[1].pipes, 2), vim.list_slice(rows[3].pipes, 2))
        end
    end)

    it("keeps a centred cell after a cell over the cap on its column centre", function()
        local tbl, conceals, layout = measure({
            "| name | b | c |",
            "| - | :-: | - |",
            "| " .. ("a"):rep(34) .. " | 1 | 2 |",
            "| " .. ("x"):rep(28) .. " | " .. ("b"):rep(30) .. " | 3 |",
        })
        for _, level in ipairs({ 0, 3 }) do
            local rows = screen(tbl, conceals, layout, level)
            assert.are.equal(rows[1].starts[2], rows[3].starts[2])
            assert.are.same(vim.list_slice(rows[1].pipes, 2), vim.list_slice(rows[3].pipes, 2))
        end
    end)

    it("needs no pads at conceallevel 0 on a table formatted with a cell over the cap", function()
        set({
            "| a | b | c |",
            "| - | --: | - |",
            "| " .. ("a"):rep(34) .. " | 1 | 2 |",
            "| x | " .. ("b"):rep(30) .. " | 3 |",
        })
        md.align(buf, 0, 3, OPTS)
        local tbl, conceals, layout = measure(get())
        assert.are.same({}, md.pads(tbl, layout, conceals, 0))
    end)

    --- Pads of `lines` after formatting them with `OPTS`, at `level`.
    --- @param lines string[]
    --- @param level integer
    --- @return { row: integer, col: integer, width: integer }[]
    local function formatted_pads(lines, level)
        set(lines)
        md.align(buf, 0, #lines - 1, OPTS)
        local tbl, conceals, layout = measure(get())
        return md.pads(tbl, layout, conceals, level)
    end

    it("needs no pads at conceallevel 0 on a formatted table with a centred column", function()
        assert.are.same({}, formatted_pads({ "| a | b |", "|---|:-:|", "| 1 | 1 |", "| 1 | 2222 |" }, 0))
    end)

    it("needs no pads at conceallevel 0 on a formatted table with empty and missing cells", function()
        assert.are.same({}, formatted_pads({ "| a | b | c |", "|---|---|---|", "| 1 |", "| 1 || 3 | 4 |" }, 0))
    end)

    it("needs no pads at conceallevel 0 on a formatted table with markup at cell edges", function()
        assert.are.same({}, formatted_pads(EXAMPLE, 0))
    end)

    it("needs no pads at conceallevel 0 and 1 on a formatted table with code spans", function()
        local lines = { "| name | b |", "|:-:|--:|", "| `x` | `yy` |", "| plain | 1 |" }
        assert.are.same({}, formatted_pads(lines, 0))
        assert.are.same({}, formatted_pads(lines, 1))
    end)

    it("adds no pads at conceallevel 0 to a table without spaces next to its pipes", function()
        local tbl, conceals, layout = measure({ "|a|b|", "|-|-|" })
        assert.are.same({}, md.pads(tbl, layout, conceals, 0))
    end)

    it("lines up the ends of a right-aligned column", function()
        local tbl, conceals, layout = measure({
            "| n |",
            "| --: |",
            "| **12** |",
            "| 345 |",
            "|    6 |",
        })
        for level = 0, 3 do
            assert_aligned(screen(tbl, conceals, layout, level), nil, { 1 })
        end
    end)

    it("lines up a struck-through cell over the strikethrough override", function()
        local tbl, conceals, layout = measure({
            "| a | b |",
            "| - | - |",
            "| ~~gone~~ | 1 |",
            "| ~one~ | 2 |",
            "| kept | 3 |",
        })
        for level = 0, 3 do
            assert_aligned(screen(tbl, conceals, layout, level), { 1, 2 })
        end
    end)

    it("keeps an escaped pipe in its cell", function()
        local tbl, conceals, layout = measure({
            "| a | b |",
            "| - | - |",
            "| x \\| y | 1 |",
            "| plain | 2 |",
        })
        for level = 1, 2 do
            assert_aligned(screen(tbl, conceals, layout, level), { 1, 2 })
        end
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
            "|a   |b  |",
            "|----|---|",
            "|long cell|x|",
        }, get())
    end)
end)
