---@diagnostic disable: undefined-global
local ws = require "utils/scroll"

--- Line i is 100 `a`s, then `<Tab>Bi<Tab>C`. Line 3 is `x<Tab>B2<Tab>C2`.
--- With `vartabstop=102,4,4`, the first Tab of line 3 spans screen columns
--- 2..102.
--- @param n integer|nil number of lines, default 40
--- @return string[]
local function tsv_lines(n)
    local lines = {}
    for i = 1, n or 40 do
        lines[i] = string.rep("a", 100) .. "\tB" .. i .. "\tC"
    end
    lines[3] = "x\tB2\tC2"
    return lines
end

--- Show `lines` in a fresh scratch buffer in the current window, `nowrap`,
--- `vartabstop=102,4,4`, no scrolloff/sidescrolloff.
--- @param lines string[]
--- @return integer winid
local function scratch(lines)
    local buf = vim.api.nvim_create_buf(false, true)
    vim.api.nvim_buf_set_lines(buf, 0, -1, false, lines)
    vim.bo[buf].vartabstop = "102,4,4"
    local win = vim.api.nvim_get_current_win()
    vim.api.nvim_win_set_buf(win, buf)
    vim.wo[win].wrap = false
    vim.wo[win].number = false
    vim.wo[win].scrolloff = 0
    vim.wo[win].sidescrolloff = 0
    vim.wo[win].foldmethod = "manual"
    vim.wo[win].smoothscroll = false
    return win
end

local function view(win)
    return vim.api.nvim_win_call(win, vim.fn.winsaveview)
end

local getmousepos = vim.fn.getmousepos

--- Make `getmousepos()` report `win` as the window under the mouse.
local function mouse_over(win)
    vim.fn.getmousepos = function() return { winid = win } end
end

describe("scroll", function()
    before_each(function()
        vim.o.columns = 80
        vim.o.lines = 24
        vim.o.mousescroll = "ver:1,hor:1"
        vim.cmd "silent! only!"
    end)

    after_each(function()
        vim.fn.getmousepos = getmousepos
    end)

    describe("nearest_line", function()
        local nearest_line = ws.nearest_line

        --- @param hits table<integer, boolean>
        local function finder(hits)
            return function(lnum) return hits[lnum] and lnum * 10 or nil end
        end

        it("returns the cursor line when it qualifies", function()
            local lnum, value = nearest_line(finder { [1] = true, [2] = true }, 2, 1, 3, 2)
            assert.are.equal(2, lnum)
            assert.are.equal(20, value)
        end)

        it("picks the nearest qualifying line", function()
            assert.are.equal(4, nearest_line(finder { [1] = true, [4] = true }, 3, 1, 4, 2.5))
        end)

        it("breaks distance ties toward the centre", function()
            local hits = { [1] = true, [3] = true }
            assert.are.equal(3, nearest_line(finder(hits), 2, 1, 5, 3))
            assert.are.equal(1, nearest_line(finder(hits), 2, 1, 5, 1))
        end)

        it("returns nil when nothing in the band qualifies", function()
            assert.is_nil(nearest_line(finder {}, 2, 1, 3, 2))
        end)

        it("ignores qualifying lines outside the band", function()
            assert.is_nil(nearest_line(finder { [1] = true, [5] = true }, 3, 2, 4, 3))
        end)

        it("searches the band from a cursor line outside it", function()
            assert.are.equal(5, nearest_line(finder { [5] = true, [8] = true }, 1, 5, 8, 6))
        end)

        it("does not call find beyond the first hit", function()
            local called = {}
            local function find(lnum)
                called[#called + 1] = lnum
                return lnum == 4 or nil
            end
            assert.are.equal(4, nearest_line(find, 3, 1, 9, 5))
            assert.are.same({ 3, 4 }, called)
        end)
    end)

    describe("col_in_band", function()
        local col_in_band = ws.col_in_band

        it("returns the character under want when it fits", function()
            local win = scratch(tsv_lines())
            assert.are.equal(9, col_in_band(win, 1, 10, 1, 80))
        end)

        it("returns the first fitting character right of the band start", function()
            local win = scratch(tsv_lines())
            assert.are.equal(3, col_in_band(win, 1, 1, 4, 83))
        end)

        it("returns the last fitting character left of the band end", function()
            local win = scratch(tsv_lines())
            assert.are.equal(99, col_in_band(win, 2, 103, 22, 101))
        end)

        it("skips a Tab partly left of the band", function()
            -- Line 3: x(1) Tab(2..102) B(103) 2(104) …
            local win = scratch(tsv_lines())
            assert.are.equal(2, col_in_band(win, 3, 1, 4, 110))
        end)

        it("steps off a Tab partly right of the band", function()
            local win = scratch(tsv_lines())
            assert.are.equal(0, col_in_band(win, 3, 61, 1, 80))
        end)

        it("returns nil when a Tab is wider than the band", function()
            local win = scratch(tsv_lines())
            assert.is_nil(col_in_band(win, 3, 1, 4, 83))
            assert.is_nil(col_in_band(win, 3, 103, 22, 101))
        end)

        it("fits an empty line only at the first column", function()
            local win = scratch { "" }
            assert.are.equal(0, col_in_band(win, 1, 5, 1, 80))
            assert.is_nil(col_in_band(win, 1, 5, 2, 81))
        end)

        it("skips a double-width character cut by either band edge", function()
            -- 日(1..2) 本(3..4) 語(5..6)
            local win = scratch { "日本語" }
            assert.are.equal(3, col_in_band(win, 1, 1, 2, 10))
            assert.are.equal(3, col_in_band(win, 1, 6, 1, 5))
        end)

        it("clamps want to the line end", function()
            local win = scratch(tsv_lines())
            local last = #tsv_lines()[1] - 1
            assert.are.equal(last, col_in_band(win, 1, vim.v.maxcol, 30, 200))
        end)
    end)

    describe("fit_leftcol", function()
        local fit_leftcol = ws.fit_leftcol

        it("returns leftcol when a character fits there", function()
            local win = scratch(tsv_lines())
            assert.are.equal(3, fit_leftcol(win, 1, 3, "l"))
            assert.are.equal(3, fit_leftcol(win, 1, 3, "h"))
        end)

        it("skips a Tab wider than the window", function()
            -- Line 3: x(1) Tab(2..102) B(103); B fits from leftcol 23 = 103 − 80.
            local win = scratch(tsv_lines())
            assert.are.equal(23, fit_leftcol(win, 3, 3, "l"))
            assert.are.equal(0, fit_leftcol(win, 3, 21, "h"))
        end)

        it("returns nil when no character fits in that direction", function()
            local win = scratch { "x\t" }
            assert.is_nil(fit_leftcol(win, 1, 3, "l"))
        end)

        it("stops at inline virtual text after the last character", function()
            local win = scratch { "ab" }
            local ns = vim.api.nvim_create_namespace "scroll_spec"
            vim.api.nvim_buf_set_extmark(0, ns, 0, 2,
                { virt_text = { { ("v"):rep(10) } }, virt_text_pos = "inline" })
            assert.is_nil(fit_leftcol(win, 1, 2, "l"))
        end)
    end)

    describe("place", function()
        local place = ws.place

        it("keeps the cursor line and applies leftcol", function()
            local win = scratch(tsv_lines())
            assert.is_true(place(win, 10, 1, 1))
            local v = view(win)
            assert.are.same({ 1, 10, 10, 10 }, { v.lnum, v.col, v.leftcol, v.curswant })
        end)

        it("moves to the nearest line with a fitting character", function()
            local win = scratch(tsv_lines())
            assert.is_true(place(win, 3, 3, 1, 0))
            vim.cmd "redraw"
            local v = view(win)
            assert.are.same({ 4, 3, 3, 0 }, { v.lnum, v.col, v.leftcol, v.curswant })
        end)

        it("honours sidescrolloff", function()
            local win = scratch(tsv_lines())
            vim.wo[win].sidescrolloff = 5
            assert.is_true(place(win, 10, 1, 1))
            assert.are.equal(15, view(win).col)
            assert.is_true(place(win, 0, 1, 1))
            assert.are.equal(0, view(win).col)
            assert.is_true(place(win, 0, 1, 80))
            assert.are.equal(74, view(win).col)
        end)

        it("narrows the band by the number column", function()
            local win = scratch(tsv_lines())
            assert.is_true(place(win, 0, 1, 80))
            assert.are.equal(79, view(win).col)
            vim.wo[win].number = true
            assert.is_true(place(win, 0, 1, 80))
            assert.are.equal(80 - vim.fn.getwininfo(win)[1].textoff - 1, view(win).col)
        end)

        it("skips lines in a closed fold", function()
            local lines = tsv_lines()
            lines[1] = "x\tB\tC"
            local win = scratch(lines)
            vim.cmd "2,4fold"
            assert.is_true(place(win, 3, 1, 1))
            assert.are.equal(5, view(win).lnum)
        end)

        it("keeps the cursor within scrolloff", function()
            local win = scratch(tsv_lines())
            vim.wo[win].scrolloff = 999
            vim.fn.winrestview { topline = 10, lnum = 20 }
            vim.cmd "redraw"
            assert.is_true(place(win, 3, 1, 1))
            vim.cmd "redraw"
            local v = view(win)
            assert.are.equal(10, v.topline)
            assert.are.equal(3, v.leftcol)
        end)

        it("returns false when sidescrolloff leaves no column", function()
            local win = scratch(tsv_lines())
            vim.wo[win].sidescrolloff = 40
            local before = view(win)
            assert.is_false(place(win, 3, 1, 1))
            assert.are.same(before, view(win))
        end)
    end)

    describe("hscroll", function()
        it("scrolls right by moving to a neighbour line", function()
            local win = scratch(tsv_lines())
            mouse_over(win)
            vim.api.nvim_win_set_cursor(win, { 3, 0 })
            ws.hscroll("l")()
            vim.cmd "redraw"
            local v = view(win)
            assert.are.equal(3, v.leftcol)
            assert.are_not.equal(3, v.lnum)
        end)

        it("falls back to a line in the scrolloff margin", function()
            local lines = {}
            for i = 1, 40 do lines[i] = "x\tB\tC" end
            lines[12] = string.rep("a", 100)
            local win = scratch(lines)
            vim.wo[win].scrolloff = 5
            mouse_over(win)
            vim.fn.winrestview { topline = 10, lnum = 20, col = 0 }
            vim.cmd "redraw"
            ws.hscroll("l")()
            vim.cmd "redraw"
            local v = view(win)
            assert.are.same({ 12, 3, 3 }, { v.lnum, v.col, v.leftcol })
            assert.are.equal(7, v.topline)
        end)

        it("jumps right to the next fitting character, at the right edge", function()
            local lines = {}
            for i = 1, 40 do lines[i] = "x\tB2\tC2" end
            local win = scratch(lines)
            mouse_over(win)
            ws.hscroll("l")()
            vim.cmd "redraw"
            local v = view(win)
            -- B at column 103 = leftcol + 80
            assert.are.same({ 1, 2, 23 }, { v.lnum, v.col, v.leftcol })

            vim.wo[win].sidescrolloff = 5
            vim.fn.winrestview { col = 0, leftcol = 0 }
            ws.hscroll("l")()
            vim.cmd "redraw"
            assert.are.same({ 2, 28 }, { view(win).col, view(win).leftcol })
        end)

        it("keeps a cursor that still fits, with its curswant", function()
            local win = scratch(tsv_lines())
            mouse_over(win)
            vim.fn.winrestview { lnum = 1, col = 50, curswant = vim.v.maxcol }
            ws.hscroll("l")()
            local v = view(win)
            assert.are.same({ 1, 50, 3, vim.v.maxcol }, { v.lnum, v.col, v.leftcol, v.curswant })
        end)

        it("keeps a cursor after the end of the line", function()
            local win = scratch { string.rep("a", 50) }
            vim.wo[win].virtualedit = "onemore"
            mouse_over(win)
            vim.fn.winrestview { lnum = 1, col = 50 }
            ws.hscroll("l")()
            local v = view(win)
            assert.are.same({ 50, 3 }, { v.col, v.leftcol })
        end)

        it("scrolls left by the count", function()
            local win = scratch(tsv_lines())
            mouse_over(win)
            vim.fn.winrestview { lnum = 3, col = 2, leftcol = 24, curswant = 102 }
            vim.cmd "redraw"
            ws.hscroll("h")()
            vim.cmd "redraw"
            assert.are.equal(21, view(win).leftcol)
        end)

        it("does nothing rightward when nothing fits", function()
            local win = scratch { "x", "y", "z" }
            mouse_over(win)
            local before = view(win)
            ws.hscroll("l")()
            assert.are.same(before, view(win))
        end)

        it("jumps left as little as a visible line allows", function()
            -- Line 1: b 103..107, Tab 108..162, C 163. Line 2: b 103..117, Tab 118..162, C 163.
            vim.o.columns = 40
            local win = scratch { "x\t" .. ("b"):rep(5) .. "\tC", "x\t" .. ("b"):rep(15) .. "\tC" }
            vim.bo.vartabstop = "102,60,4"
            mouse_over(win)
            vim.fn.winrestview { lnum = 1, col = 8, leftcol = 123 }
            vim.cmd "redraw"
            ws.hscroll("h")()
            vim.cmd "redraw"
            local v = view(win)
            assert.are.same({ 2, 16, 116 }, { v.lnum, v.col, v.leftcol })
        end)

        it("jumps left over a Tab wider than the window", function()
            local win = scratch { "x\tB2\tC2" }
            mouse_over(win)
            vim.fn.winrestview { lnum = 1, col = 2, leftcol = 24, curswant = 102 }
            vim.cmd "redraw"
            ws.hscroll("h")()
            vim.cmd "redraw"
            assert.are.equal(0, view(win).leftcol)
        end)
    end)

    describe("vscroll", function()
        it("keeps leftcol when the cursor is pushed onto a Tab", function()
            local win = scratch(tsv_lines())
            mouse_over(win)
            vim.fn.winrestview { topline = 2, lnum = 2, col = 60, curswant = 60, leftcol = 0 }
            vim.cmd "redraw"
            ws.vscroll("j")()
            vim.cmd "redraw"
            local v = view(win)
            assert.are.same({ 3, 0, 3 }, { v.topline, v.leftcol, v.lnum })
            vim.cmd.normal { "j", bang = true }
            assert.are.equal(60, view(win).col)
        end)

        it("keeps leftcol when the cursor is pushed onto a short line", function()
            local lines = tsv_lines()
            lines[3] = ""
            local win = scratch(lines)
            mouse_over(win)
            vim.fn.winrestview { topline = 2, lnum = 2, col = 60, curswant = 60, leftcol = 30 }
            vim.cmd "redraw"
            ws.vscroll("j")()
            vim.cmd "redraw"
            local v = view(win)
            assert.are.same({ 3, 30, 4, 60 }, { v.topline, v.leftcol, v.lnum, v.col })
        end)

        it("changes leftcol as little as needed when nothing fits it", function()
            -- Short lines: b 103..107, Tab 108..162 (wider than the window), C 163.
            vim.o.columns = 40
            local lines = {}
            for i = 1, 29 do lines[i] = "x\t" .. ("b"):rep(5) .. "\tC" end
            lines[30] = "x\t" .. ("b"):rep(50) .. "\tC"
            local win = scratch(lines)
            vim.bo.vartabstop = "102,60,4"
            mouse_over(win)
            vim.fn.winrestview { topline = 9, lnum = 30, col = 39, curswant = 139, leftcol = 120 }
            vim.cmd "redraw"
            ws.vscroll("k")()
            vim.cmd "redraw"
            local v = view(win)
            -- C at the right edge: 163 = 123 + 40
            assert.are.same({ 8, 29, 8, 123, 139 }, { v.topline, v.lnum, v.col, v.leftcol, v.curswant })
        end)

        it("goes to leftcol 0 when only empty lines are in view", function()
            local lines = { string.rep("a", 100) }
            for i = 2, 40 do lines[i] = "" end
            local win = scratch(lines)
            mouse_over(win)
            vim.fn.winrestview { topline = 1, lnum = 1, col = 60, leftcol = 20 }
            vim.cmd "redraw"
            ws.vscroll("j")()
            vim.cmd "redraw"
            local v = view(win)
            assert.are.same({ 2, 0, 2 }, { v.topline, v.leftcol, v.lnum })
        end)

        it("scrolls a wrap window plainly", function()
            local win = scratch(tsv_lines())
            vim.wo[win].wrap = true
            mouse_over(win)
            ws.vscroll("j")()
            assert.are.equal(2, view(win).topline)
        end)

        it("scrolls the window under the mouse, not the current one", function()
            local other = scratch(tsv_lines())
            vim.cmd "vsplit"
            local current = vim.api.nvim_get_current_win()
            vim.cmd "redraw"
            local before = view(current)
            mouse_over(other)
            ws.vscroll("j")()
            assert.are.equal(current, vim.api.nvim_get_current_win())
            assert.are.same(before, view(current))
            assert.are.equal(2, view(other).topline)
        end)

        it("does nothing with ver:0", function()
            local win = scratch(tsv_lines())
            mouse_over(win)
            vim.o.mousescroll = "ver:0,hor:1"
            ws.vscroll("j")()
            assert.are.equal(1, view(win).topline)
        end)

        it("scrolls by the ver count", function()
            local win = scratch(tsv_lines())
            mouse_over(win)
            vim.o.mousescroll = "ver:5,hor:1"
            ws.vscroll("j")()
            assert.are.equal(6, view(win).topline)
            ws.vscroll("k")()
            assert.are.equal(1, view(win).topline)
        end)

        it("scrolls past the end of the buffer keeping leftcol", function()
            local win = scratch(tsv_lines())
            mouse_over(win)
            vim.fn.winrestview { topline = 39, lnum = 39, col = 60, curswant = 60, leftcol = 20 }
            vim.cmd "redraw"
            ws.vscroll("j")()
            vim.cmd "redraw"
            local v = view(win)
            assert.are.same({ 40, 20, 40 }, { v.topline, v.leftcol, v.lnum })
        end)

        -- Input is consumed only by the main loop, which a running spec never
        -- reaches, so the event goes to a child nvim.
        it("handles a real wheel event through a mapping", function()
            local addr = vim.fn.tempname()
            local child = vim.system { "nvim", "--headless", "--clean", "-i", "NONE", "--listen", addr }
            assert.is_true(vim.wait(5000, function() return vim.uv.fs_stat(addr) ~= nil end))
            local chan = vim.fn.sockconnect("pipe", addr, { rpc = true })
            local function lua(code, ...)
                return vim.rpcrequest(chan, "nvim_exec_lua", code, { ... })
            end
            lua([[
                local rtp, lines = ...
                vim.opt.rtp:prepend(rtp)
                vim.o.columns, vim.o.lines = 80, 24
                vim.api.nvim_buf_set_lines(0, 0, -1, false, lines)
                vim.bo.vartabstop = "102,4,4"
                vim.wo.wrap, vim.wo.scrolloff, vim.o.mousescroll = false, 0, "ver:1,hor:1"
                vim.keymap.set("n", "<ScrollWheelDown>", require("utils/scroll").vscroll("j"))
                vim.fn.winrestview { topline = 2, lnum = 2, col = 60, curswant = 60, leftcol = 0 }
                vim.cmd "redraw"
            ]], vim.fn.getcwd(), tsv_lines())
            vim.rpcrequest(chan, "nvim_input_mouse", "wheel", "down", "", 0, 5, 10)
            local v = lua "vim.cmd 'redraw' return vim.fn.winsaveview()"
            child:kill "sigterm"
            assert.are.same({ 3, 0, 3 }, { v.topline, v.leftcol, v.lnum })
        end)
    end)
end)
