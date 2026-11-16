--- Wheel and trackpad scroll (`<ScrollWheel…>` keys) that moves the view by
--- the requested amount, keeps `leftcol`, and moves the cursor only as far as
--- necessary.
---
--- With 'nowrap', Vim keeps the whole screen span of the cursor character
--- inside the window, less 'sidescrolloff' at each side. If the span is not
--- inside, Vim moves the cursor or changes `leftcol`. Thus a vertical scroll
--- onto a short line or a wide Tab moves the view sideways, and `zl`/`zh` are
--- blocked or jump too far.
--- Rationale: neovim skill, `references/rendering.md`, section "Horizontal
--- scroll (`leftcol`) is clamped by the cursor line".
local set_view = require("utils/init").set_view

local M = {}

--- Columns per horizontal wheel step.
local count = 3

--- Line in [lo, hi] nearest to `curline` for which `find` is non-nil.
--- If two lines are at equal distance, the line nearer `centre` wins.
--- Calls `find` in order of distance and stops at the first hit.
--- @generic T
--- @param find fun(lnum: integer): T|nil
--- @param curline integer can be outside [lo, hi]
--- @param lo integer inclusive
--- @param hi integer inclusive
--- @param centre number line that breaks distance ties
--- @return integer|nil lnum
--- @return T|nil value result of `find(lnum)`
function M.nearest_line(find, curline, lo, hi, centre)
    for d = 0, math.max(curline - lo, hi - curline) do
        local near, far = curline - d, curline + d
        if math.abs(far - centre) < math.abs(near - centre) then
            near, far = far, near
        end
        for _, lnum in ipairs(d == 0 and { near } or { near, far }) do
            if lnum >= lo and lnum <= hi then
                local value = find(lnum)
                if value ~= nil then return lnum, value end
            end
        end
    end
end

--- Character of `lnum` in `winid` that covers screen column `vcol`; past the
--- line end, the last character. An empty line has one character slot at
--- column 1.
--- @return integer col 1-indexed byte column
--- @return integer start 1-indexed screen column
--- @return integer end 1-indexed screen column
local function char_span(winid, lnum, vcol)
    local col = math.max(vim.fn.virtcol2col(winid, lnum, vcol), 1)
    local span = vim.fn.virtcol({ lnum, col }, true, winid)
    return col, span[1], span[2]
end

local function folded(winid, lnum)
    return vim.api.nvim_win_call(winid, function()
        return vim.fn.foldclosed(lnum) ~= -1
    end)
end

--- Character on `lnum` whose whole screen span is in [vlo, vhi], nearest to
--- screen column `want`.
--- @param winid integer
--- @param lnum integer 1-indexed
--- @param want integer 1-indexed screen column
--- @param vlo integer 1-indexed screen column, inclusive
--- @param vhi integer 1-indexed screen column, inclusive
--- @return integer|nil col 0-indexed byte column, nil if none fits
function M.col_in_band(winid, lnum, want, vlo, vhi)
    local width = vim.fn.virtcol({ lnum, "$" }, false, winid) - 1
    local col, s, e = char_span(winid, lnum, math.max(math.min(want, vhi, width), vlo, 1))
    if s < vlo then
        col, s, e = char_span(winid, lnum, e + 1)
    elseif e > vhi then
        col, s, e = char_span(winid, lnum, s - 1)
    end
    if s >= vlo and e <= vhi then return col - 1 end
end

--- Screen columns where Vim lets the cursor character be with this
--- `leftcol`: the text area of `winid`, less 'sidescrolloff' at each side.
--- @param winid integer
--- @param leftcol integer
--- @return integer vlo 1-indexed screen column, inclusive
--- @return integer vhi 1-indexed screen column, inclusive
function M.cursor_band(winid, leftcol)
    local info = vim.fn.getwininfo(winid)[1]
    local siso = vim.api.nvim_get_option_value("sidescrolloff", { win = winid })
    -- 'sidescrolloff' is suspended at the start of the line.
    return leftcol + 1 + (leftcol > 0 and siso or 0),
        leftcol + info.width - info.textoff - siso
end

--- `leftcol` nearest to `leftcol` in direction `dir`, `leftcol` included, at
--- which a character of `lnum` is fully visible in `winid`.
--- @param winid integer
--- @param lnum integer 1-indexed
--- @param leftcol integer
--- @param dir "h"|"l" toward smaller or larger leftcol
--- @return integer|nil leftcol nil if no character fits in that direction
function M.fit_leftcol(winid, lnum, leftcol, dir)
    local info = vim.fn.getwininfo(winid)[1]
    local width = info.width - info.textoff
    local siso = vim.api.nvim_get_option_value("sidescrolloff", { win = winid })
    local vlo, vhi = M.cursor_band(winid, leftcol)
    -- Walk characters away from the view in `dir`; the needed leftcol only
    -- grows (rightward) or shrinks (leftward), so the first fit is nearest.
    local col, s, e = char_span(winid, lnum, dir == "l" and vlo or vhi)
    while true do
        local fit = dir == "l" and math.max(leftcol, e - width + siso)
            or math.max(0, math.min(leftcol, s - 1 - siso))
        local fit_lo, fit_hi = M.cursor_band(winid, fit)
        if s >= fit_lo and e <= fit_hi then return fit end
        local next_col, next_s, next_e = char_span(winid, lnum, dir == "l" and e + 1 or s - 1)
        if next_col == col then return end
        col, s, e = next_col, next_s, next_e
    end
end

--- Visible lines of `winid` at least `so` lines from the window edges. At the
--- start and end of the buffer the edge lines count too, as for 'scrolloff'.
--- @param winid integer
--- @param so integer|nil lines; nil: 'scrolloff'
--- @return integer lo 1-indexed, inclusive
--- @return integer hi 1-indexed, inclusive
function M.line_band(winid, so)
    local info = vim.fn.getwininfo(winid)[1]
    local last = vim.api.nvim_buf_line_count(vim.api.nvim_win_get_buf(winid))
    so = math.min(
        so or vim.api.nvim_get_option_value("scrolloff", { win = winid }),
        math.floor((info.height - 1) / 2))
    local lo = info.topline == 1 and 1 or math.min(info.topline + so, last)
    return lo, info.botline >= last and last or math.max(info.botline - so, lo)
end

--- `leftcol` nearest to `leftcol`, `leftcol` included, at which a line of
--- `line_band(winid, so)` outside closed folds has a fully visible character.
--- @param winid integer
--- @param leftcol integer
--- @param dir "h"|"l"|nil search toward smaller or larger leftcol; nil: both
--- @param so integer|nil see `line_band`
--- @return integer|nil leftcol nil if no such line
function M.view_fit_leftcol(winid, leftcol, dir, so)
    local lo, hi = M.line_band(winid, so)
    local nearest
    for lnum = lo, hi do
        for _, d in ipairs(dir and { dir } or { "h", "l" }) do
            local fit = not folded(winid, lnum) and M.fit_leftcol(winid, lnum, leftcol, d) or nil
            if fit and (nearest == nil or math.abs(fit - leftcol) < math.abs(nearest - leftcol)) then
                nearest = fit
                if nearest == leftcol then return nearest end
            end
        end
    end
    return nearest
end

--- Set `leftcol` of `winid` and move the cursor to a character that is fully
--- visible with that `leftcol`. The cursor goes to the line of
--- `line_band(winid, so)` nearest `curline` that has such a character, and
--- to the character nearest `want` on that line. Lines in closed folds are
--- skipped. With `so` below 'scrolloff', Vim can then scroll the window
--- vertically to keep the cursor 'scrolloff' lines from the edge.
--- @param winid integer
--- @param leftcol integer
--- @param curline integer
--- @param want integer 1-indexed screen column
--- @param curswant integer|nil nil: start of the new cursor span, 0-indexed
--- @param so integer|nil see `line_band`
--- @return boolean placed false if no line qualifies (nothing changed)
function M.place(winid, leftcol, curline, want, curswant, so)
    local info = vim.fn.getwininfo(winid)[1]
    local vlo, vhi = M.cursor_band(winid, leftcol)
    local lo, hi = M.line_band(winid, so)

    local function find(lnum)
        if not folded(winid, lnum) then return M.col_in_band(winid, lnum, want, vlo, vhi) end
    end

    local lnum, col = M.nearest_line(find, curline, lo, hi, (info.topline + info.botline) / 2)
    if lnum == nil then return false end
    set_view {
        winid = winid,
        lnum = lnum,
        col = col,
        leftcol = leftcol,
        curswant = curswant or vim.fn.virtcol({ lnum, col + 1 }, true, winid)[1] - 1,
    }
    return true
end

--- Set `leftcol` of `winid` to `view_fit_leftcol(winid, leftcol, dir, so)` and
--- move the cursor there with `place`, preferring the 'scrolloff' band.
--- @param winid integer
--- @param leftcol integer
--- @param dir "h"|"l"|nil see `view_fit_leftcol`
--- @param so integer|nil see `line_band`
--- @param cursor { lnum: integer, want: integer, curswant: integer|nil } see `place`
--- @return boolean placed false if no line qualifies (nothing changed)
function M.settle(winid, leftcol, dir, so, cursor)
    local fit = M.view_fit_leftcol(winid, leftcol, dir, so)
    if fit == nil then return false end
    return M.place(winid, fit, cursor.lnum, cursor.want, cursor.curswant)
        or M.place(winid, fit, cursor.lnum, cursor.want, cursor.curswant, so)
end

--- Build a keymap rhs that scrolls the window under the mouse horizontally
--- by `count` columns. Does nothing with 'wrap'. The cursor moves to a fully
--- visible character, preferably on a line in the 'scrolloff' band, else on
--- any visible line (Vim then scrolls vertically to restore 'scrolloff').
--- If no visible line has one, the scroll goes further in the same direction,
--- as little as needed for one to have it (`view_fit_leftcol`). If none has
--- one at all, a leftward scroll falls back to `zh`.
--- @param dir "h"|"l" left or right
--- @return fun() rhs
function M.hscroll(dir)
    return function()
        local winid = vim.fn.getmousepos().winid
        if winid == 0 or vim.wo[winid].wrap then return end
        local view = vim.api.nvim_win_call(winid, vim.fn.winsaveview)
        local leftcol = math.max(0, view.leftcol + (dir == "h" and -count or count))
        if leftcol == view.leftcol then return end
        local span = vim.fn.virtcol({ view.lnum, view.col + 1 }, true, winid)
        local vlo, vhi = M.cursor_band(winid, leftcol)
        if span[1] >= vlo and span[2] <= vhi then
            set_view { winid = winid, leftcol = leftcol }
            return
        end
        if not M.settle(winid, leftcol, dir, 0, { lnum = view.lnum, want = span[1] })
            and dir == "h"
        then
            -- Lets the view always reach leftcol 0.
            vim.api.nvim_win_call(winid, function()
                vim.cmd.normal { count .. "zh", bang = true }
            end)
        end
    end
end

--- Build a keymap rhs that scrolls the window under the mouse vertically by
--- the `ver:` count of 'mousescroll'. With 'nowrap', keeps `leftcol` if a
--- line in the 'scrolloff' band has a fully visible character, else changes
--- it as little as needed for one to have it (`settle`).
--- @param dir "k"|"j" up or down
--- @return fun() rhs
function M.vscroll(dir)
    local key = vim.keycode(dir == "j" and "<C-e>" or "<C-y>")
    return function()
        local winid = vim.fn.getmousepos().winid
        if winid == 0 then return end
        local n = tonumber(vim.o.mousescroll:match "ver:(%d+)") or 3
        if n == 0 then return end
        local old = vim.api.nvim_win_call(winid, vim.fn.winsaveview)
        local new = vim.api.nvim_win_call(winid, function()
            vim.cmd.normal { n .. key, bang = true }
            return vim.fn.winsaveview()
        end)
        if vim.wo[winid].wrap or new.leftcol == old.leftcol then return end
        M.settle(winid, old.leftcol, nil, nil,
            { lnum = new.lnum, want = old.curswant + 1, curswant = old.curswant })
    end
end

return M
