---@diagnostic disable: undefined-global
local util = require "utils/init"

--- Replace the buffer with `lines`, then make and leave a selection with `keys`
--- typed from the start of the first line.
local function make_selection(lines, keys)
    vim.api.nvim_buf_set_lines(0, 0, -1, true, lines)
    vim.api.nvim_win_set_cursor(0, { 1, 0 })
    vim.cmd("normal! " .. vim.keycode(keys .. "<Esc>"))
end

--- Text of the last visual selection, from the first byte to the last byte inclusive.
local function range_text()
    local r1, c1, r2, c2 = util.get_visual_range()
    return vim.api.nvim_buf_get_text(0, r1 - 1, c1, r2 - 1, c2 + 1, {})
end

describe("get_visual_range", function()
    before_each(function() vim.cmd.enew { bang = true } end)

    local cases = {
        { "ASCII", { "abc" }, "vl", { 1, 0, 1, 1 }, { "ab" } },
        { "3-byte last char", { "a€" }, "vl", { 1, 0, 1, 3 }, { "a€" } },
        { "2-byte only char", { "é" }, "v", { 1, 0, 1, 1 }, { "é" } },
        { "v$", { "abc" }, "v$", { 1, 0, 1, 2 }, { "abc" } },
        { "v$ on 3-byte last char", { "a€" }, "v$", { 1, 0, 1, 3 }, { "a€" } },
        { "V", { "abc", "de" }, "lV", { 1, 0, 1, 2 }, { "abc" } },
        { "V ending on an empty line", { "abc", "" }, "Vj", { 1, 0, 2, -1 }, { "abc", "" } },
        { "V starting on an empty line", { "", "abc" }, "Vj", { 1, 0, 2, 2 }, { "", "abc" } },
        { "<C-v> ending on 3-byte char", { "abc", "x€y" }, "<C-v>jl", { 1, 0, 2, 3 }, { "abc", "x€" } },
    }
    for _, case in ipairs(cases) do
        local name, lines, keys, range, text = unpack(case)
        it("gives the last byte of the last char for " .. name, function()
            make_selection(lines, keys)
            assert.are.same(range, { util.get_visual_range() })
            assert.are.same(text, range_text())
        end)
    end
end)

describe("get_char", function()
    before_each(function()
        vim.cmd.enew { bang = true }
        vim.api.nvim_buf_set_lines(0, 0, -1, true, { "aé€" })
    end)

    it("gives nothing at column 0", function()
        assert.are.same({ "", 0 }, { util.get_char(0, 0) })
    end)

    it("gives the char ending before the column, and its start", function()
        assert.are.same({ "a", 0 }, { util.get_char(0, 1) })
        assert.are.same({ "é", 1 }, { util.get_char(0, 3) })
        assert.are.same({ "€", 3 }, { util.get_char(0, 6) })
    end)
end)

describe("readtext", function()
    it("gives the reason when the file can't be opened", function()
        local text, err = util.readtext("/nonexistent/file")
        assert.is_nil(text)
        assert.matches("No such file", err)
    end)
end)

describe("schedule_notify", function()
    local notify = vim.notify
    local notes

    before_each(function()
        notes = {}
        vim.notify = function(msg, level) notes[#notes + 1] = { msg = msg, level = level } end
    end)

    after_each(function() vim.notify = notify end)

    it("notifies stderr at ERROR level by default", function()
        util.schedule_notify { stderr = "oops\n", stdout = "" }
        vim.wait(50, function() return #notes > 0 end)
        assert.are.same({ { msg = "oops", level = vim.log.levels.ERROR } }, notes)
    end)

    it("falls back to stdout at the given level", function()
        util.schedule_notify({ stderr = "", stdout = "done\n" }, vim.log.levels.INFO)
        vim.wait(50, function() return #notes > 0 end)
        assert.are.same({ { msg = "done", level = vim.log.levels.INFO } }, notes)
    end)
end)

describe("crosses_wrap", function()
    local win, buf
    before_each(function()
        vim.cmd("20vnew")
        win, buf = vim.api.nvim_get_current_win(), vim.api.nvim_get_current_buf()
        vim.wo[win].number = false
        vim.wo[win].signcolumn = "no"
        vim.wo[win].wrap = true
        vim.api.nvim_buf_set_lines(buf, 0, -1, true, { ("x"):rep(16) .. " $abcdefgh$ tail" })
    end)
    after_each(function() vim.api.nvim_buf_delete(buf, { force = true }) end)

    it("is true for a span crossing the row end", function()
        assert.is_true(util.crosses_wrap(win, 0, 17, 10))
    end)

    it("is false for a span within the first row", function()
        assert.is_false(util.crosses_wrap(win, 0, 0, 5))
    end)

    it("is false for a span ending exactly at the row end", function()
        assert.is_false(util.crosses_wrap(win, 0, 10, 10))
    end)

    it("is false with nowrap", function()
        vim.wo[win].wrap = false
        assert.is_false(util.crosses_wrap(win, 0, 17, 10))
    end)

    it("is unchanged by an inline mark concealing the span to the same footprint", function()
        vim.wo[win].conceallevel = 1
        local ns = vim.api.nvim_create_namespace("crosses_wrap_spec")
        vim.api.nvim_buf_set_extmark(buf, ns, 0, 17, {
            end_col = 27, conceal = "", virt_text_pos = "inline", virt_text = { { ("#"):rep(9) } },
        })
        assert.is_true(util.crosses_wrap(win, 0, 17, 10))
        vim.api.nvim_buf_set_extmark(buf, ns, 0, 0, {
            end_col = 5, conceal = "", virt_text_pos = "inline", virt_text = { { ("#"):rep(4) } },
        })
        assert.is_false(util.crosses_wrap(win, 0, 0, 5))
    end)
end)

describe("line_hl", function()
    local buf, ns_a, ns_b
    before_each(function()
        buf = vim.api.nvim_create_buf(false, true)
        vim.api.nvim_buf_set_lines(buf, 0, -1, true, { "a", "b", "c" })
        ns_a = vim.api.nvim_create_namespace("line_hl_spec_a")
        ns_b = vim.api.nvim_create_namespace("line_hl_spec_b")
    end)
    after_each(function() vim.api.nvim_buf_delete(buf, { force = true }) end)

    it("is nil on a row without marks", function()
        assert.is_nil(util.line_hl(buf, 1))
    end)

    it("returns the mark's line_hl_group", function()
        vim.api.nvim_buf_set_extmark(buf, ns_a, 1, 0, { line_hl_group = "LnA" })
        assert.are.equal("LnA", util.line_hl(buf, 1))
    end)

    it("picks the highest priority across namespaces", function()
        vim.api.nvim_buf_set_extmark(buf, ns_a, 1, 0, { line_hl_group = "LnA", priority = 20 })
        vim.api.nvim_buf_set_extmark(buf, ns_b, 1, 0, { line_hl_group = "LnB", priority = 10 })
        assert.are.equal("LnA", util.line_hl(buf, 1))
    end)

    it("picks the higher extmark id at equal priority, whatever the creation order", function()
        for _ = 1, 5 do vim.api.nvim_buf_set_extmark(buf, ns_a, 2, 0, {}) end
        vim.api.nvim_buf_set_extmark(buf, ns_a, 0, 0, { line_hl_group = "LnA" }) -- id 6
        vim.api.nvim_buf_set_extmark(buf, ns_b, 0, 0, { line_hl_group = "LnB" }) -- id 1
        assert.are.equal("LnA", util.line_hl(buf, 0))
        vim.api.nvim_buf_set_extmark(buf, ns_b, 1, 0, { line_hl_group = "LnB" }) -- id 2
        vim.api.nvim_buf_set_extmark(buf, ns_a, 1, 0, { line_hl_group = "LnA" }) -- id 7
        assert.are.equal("LnA", util.line_hl(buf, 1))
    end)

    it("ignores marks on other rows", function()
        vim.api.nvim_buf_set_extmark(buf, ns_a, 0, 0, { line_hl_group = "LnA" })
        vim.api.nvim_buf_set_extmark(buf, ns_a, 2, 0, { line_hl_group = "LnA" })
        assert.is_nil(util.line_hl(buf, 1))
    end)

    it("counts a range mark ending at column 0 of the row", function()
        vim.api.nvim_buf_set_extmark(buf, ns_a, 0, 0, { end_row = 1, end_col = 0, line_hl_group = "LnA" })
        assert.are.equal("LnA", util.line_hl(buf, 1))
    end)
end)
