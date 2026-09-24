---@diagnostic disable: undefined-global
local qf = require "utils/qf"

local lines = {
    "local x = f(x)",
    "print(x)",
    "print(x)",
    "print(x)",
    "print(x)",
    "print(x)",
}

local path, other_path

--- Quickfix item on line `lnum` from column `col` to `end_col`, both 1-based,
--- end exclusive.
--- @param filename string|nil defaults to `path`
local function item(lnum, col, end_col, filename)
    return { filename = filename or path, lnum = lnum, col = col, end_lnum = lnum, end_col = end_col, text = "" }
end

--- Item for the `x` on `print(x)` at `lnum`.
local function print_x(lnum) return item(lnum, 7, 8) end

local function cursor() return vim.api.nvim_win_get_cursor(0) end

local function qf_idx() return vim.fn.getqflist({ idx = 0 }).idx end

local function qf_size() return vim.fn.getqflist({ size = 0 }).size end

describe("qf", function()
    before_each(function()
        path = vim.fn.tempname() .. ".lua"
        other_path = vim.fn.tempname() .. ".lua"
        vim.fn.writefile(lines, path)
        vim.fn.writefile(lines, other_path)
        -- Buffer names resolve symlinks, and macOS tempnames are under the /var symlink.
        path, other_path = vim.fn.resolve(path), vim.fn.resolve(other_path)
        vim.cmd.edit(path)
        vim.cmd.clearjumps()
        vim.fn.setqflist({}, 'f')
        vim.fn.setqflist({}, ' ', { title = "previous", items = { item(1, 1, 2) } })
    end)

    after_each(function()
        vim.cmd("silent! %bwipeout!")
    end)

    describe("contains_cursor", function()
        it("is true with the cursor in the range", function()
            vim.api.nvim_win_set_cursor(0, { 1, 6 })
            assert.is_true(qf.contains_cursor(item(1, 7, 8)))
        end)

        it("is false on the same line outside the range", function()
            vim.api.nvim_win_set_cursor(0, { 1, 12 })
            assert.is_false(qf.contains_cursor(item(1, 7, 8)))
        end)

        it("is false at end_col", function()
            vim.api.nvim_win_set_cursor(0, { 1, 7 })
            assert.is_false(qf.contains_cursor(item(1, 7, 8)))
        end)

        it("is false for another file", function()
            vim.api.nvim_win_set_cursor(0, { 1, 6 })
            assert.is_false(qf.contains_cursor(item(1, 7, 8, other_path)))
        end)

        it("is true on the start of an item with no end", function()
            vim.api.nvim_win_set_cursor(0, { 2, 6 })
            assert.is_true(qf.contains_cursor { filename = path, lnum = 2, col = 7 })
        end)

        it("is true on the start of a zero-width range", function()
            vim.api.nvim_win_set_cursor(0, { 2, 6 })
            assert.is_true(qf.contains_cursor(item(2, 7, 7)))
        end)
    end)

    describe("jump_or_load", function()
        local function assert_previous_list()
            assert.are.equal("previous", vim.fn.getqflist({ title = 0 }).title)
        end

        it("keeps the previous list for 0 items", function()
            qf.jump_or_load { items = {} }
            assert_previous_list()
        end)

        it("keeps the previous list for only self", function()
            vim.api.nvim_win_set_cursor(0, { 2, 6 })
            qf.jump_or_load { items = { print_x(2) } }
            assert_previous_list()
        end)

        it("jumps to the one other item after self", function()
            vim.api.nvim_win_set_cursor(0, { 2, 6 })
            qf.jump_or_load { items = { print_x(2), print_x(4) } }
            assert.are.same({ 4, 6 }, cursor())
            assert.are.equal(2, qf_idx())
            local jumps = vim.fn.getjumplist()[1]
            assert.are.equal(2, jumps[#jumps].lnum)
        end)

        it("jumps to the one other item before self", function()
            vim.api.nvim_win_set_cursor(0, { 4, 6 })
            qf.jump_or_load { items = { print_x(2), print_x(4) } }
            assert.are.same({ 2, 6 }, cursor())
            assert.are.equal(1, qf_idx())
            local jumps = vim.fn.getjumplist()[1]
            assert.are.equal(4, jumps[#jumps].lnum)
        end)

        it("loads self + 3 others without moving or opening the window", function()
            vim.api.nvim_win_set_cursor(0, { 3, 6 })
            local win = vim.api.nvim_get_current_win()
            qf.jump_or_load { items = { print_x(2), print_x(3), print_x(4), print_x(5) } }
            assert.are.same({ 3, 6 }, cursor())
            assert.are.equal(win, vim.api.nvim_get_current_win())
            assert.are.equal(2, qf_idx())
            assert.are.equal(0, vim.fn.getqflist({ winid = 0 }).winid)
            vim.cmd.cnext()
            assert.are.same({ 4, 6 }, cursor())
        end)

        it("adds self and jumps to a single item", function()
            vim.api.nvim_win_set_cursor(0, { 2, 0 })
            qf.jump_or_load { items = { print_x(5) } }
            assert.are.same({ 5, 6 }, cursor())
            assert.are.equal(2, qf_size())
            vim.cmd.cprevious()
            assert.are.same({ 2, 0 }, cursor())
        end)

        it("inserts self in sort order when no item contains the cursor", function()
            vim.api.nvim_win_set_cursor(0, { 3, 0 })
            qf.jump_or_load { items = { print_x(2), print_x(4), print_x(5) } }
            assert.are.same({ 3, 0 }, cursor())
            assert.are.equal(4, qf_size())
            assert.are.equal(2, qf_idx())
            assert.are.equal(3, vim.fn.getqflist()[2].lnum)
            vim.cmd.cnext()
            assert.are.same({ 4, 6 }, cursor())
        end)

        it("appends self when the cursor sorts after all items", function()
            vim.api.nvim_win_set_cursor(0, { 6, 0 })
            qf.jump_or_load { items = { print_x(2), print_x(3), print_x(4) } }
            assert.are.equal(4, qf_idx())
            assert.are.equal(4, qf_size())
        end)

        it("does not change what.items", function()
            vim.api.nvim_win_set_cursor(0, { 3, 0 })
            local items = { print_x(2), print_x(4), print_x(5) }
            local copy = vim.deepcopy(items)
            qf.jump_or_load { items = items }
            assert.are.same(copy, items)
        end)

        it("counts only the ref under the cursor as self", function()
            vim.api.nvim_win_set_cursor(0, { 1, 12 })
            qf.jump_or_load { items = { item(1, 7, 8), item(1, 13, 14) } }
            assert.are.same({ 1, 6 }, cursor())
            assert.are.equal(1, qf_idx())
        end)
    end)
end)
