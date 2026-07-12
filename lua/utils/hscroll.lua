--- Horizontal scroll that walks the cursor onto a scrollable line first.
---
--- `zh`/`zl` scroll the view but clamp `leftcol` to keep the cursor on a
--- real character *of the cursor line*, so on a short/empty cursor line a
--- rightward scroll silently no-ops even with longer lines in view. The
--- mechanism here moves the cursor onto a line long enough to permit the
--- scroll, then scrolls — guarded on scrollability so a stray horizontal
--- drift during a vertical trackpad swipe cannot move the cursor.
--- Rationale: neovim skill, `references/rendering.md § Horizontal scroll`.
---
--- `leftcol` stays coupled to the cursor line, so the offset resets to 0
--- whenever the cursor lands on a shorter line by any other means (vertical
--- scroll, j/k, search). Not preserved here; that needs a change to native
--- horizontal-scroll behaviour, out of scope for this wrapper.
local M = {}

--- Shared between the qualification test and the `normal!` scroll distance
--- so they cannot drift apart.
local count = 3

--- Choose a buffer line to move onto so a rightward scroll of `need`
--- columns is not clamped by too short a cursor line.
---
--- Scrolling right, the cursor rides the left edge until it sits on the
--- line's last character, so a line permits the scroll iff its display
--- width exceeds `need`. Staying on the cursor line cannot cause vertical
--- drift, so it is preferred whenever it qualifies; otherwise the nearest
--- qualifying line in the scrolloff-safe band is chosen, ties broken toward
--- the screen centre.
---
--- @param widths integer[] display width of each visible buffer line;
---   index i is buffer line (top + i - 1), covering the viewport w0..w$
--- @param top integer 1-indexed buffer line at the viewport top (w0)
--- @param curline integer 1-indexed current cursor line
--- @param need integer columns that must be scrollable (leftcol + count);
---   a line qualifies iff its width > need
--- @param lo integer lowest candidate buffer line (scrolloff-safe, inclusive)
--- @param hi integer highest candidate buffer line (scrolloff-safe, inclusive)
--- @return integer|nil target 1-indexed line to move onto, nil if none qualifies
function M.scroll_target(widths, top, curline, need, lo, hi)
    local function qualifies(ln)
        local w = widths[ln - top + 1]
        return w ~= nil and w > need
    end

    if qualifies(curline) then
        return curline
    end

    local centre = top + (#widths - 1) / 2
    local target, best_dist, best_centre
    for ln = lo, hi do
        if qualifies(ln) then
            local dist = math.abs(ln - curline)
            local cdist = math.abs(ln - centre)
            if target == nil
                or dist < best_dist
                or (dist == best_dist and cdist < best_centre)
            then
                target, best_dist, best_centre = ln, dist, cdist
            end
        end
    end
    return target
end

--- Build a keymap rhs that scrolls the view horizontally by `count`.
--- @param dir "h"|"l" scroll direction ('h' left, 'l' right)
--- @return fun() rhs
function M.hscroll(dir)
    return function()
        if vim.wo.wrap then
            return -- leftcol is pinned at 0 when wrapping
        end

        local leftcol = vim.fn.winsaveview().leftcol

        if dir == 'h' then
            -- Leftward is never blocked: content always exists to the left
            -- once leftcol > 0, so no scan or cursor move is needed.
            if leftcol == 0 then return end
            vim.cmd('normal! ' .. count .. 'zh')
            return
        end

        local w0 = vim.fn.line('w0')
        local wbot = vim.fn.line('w$')
        local curline = vim.fn.line('.')
        local last = vim.fn.line('$')

        -- v1: virtcol ignores conceal, so widths are over-reported on a
        -- nowrap buffer with conceal active and the target may be off.
        local widths = {}
        for lnum = w0, wbot do
            widths[lnum - w0 + 1] = vim.fn.virtcol({ lnum, '$' }) - 1
        end

        -- Keep the move inside the vertical viewport (mirrors set_view's
        -- clamp) so relocating the cursor can't trigger a vertical scroll
        -- to satisfy scrolloff.
        local scrolloff = vim.api.nvim_get_option_value('scrolloff', { win = 0 })
        local lo = math.min(w0 + scrolloff, last)
        local hi = math.max(math.min(wbot - scrolloff, last), lo)

        local target = M.scroll_target(widths, w0, curline, leftcol + count, lo, hi)
        if target == nil then return end

        local delta = target - curline
        local motion = delta > 0 and (delta .. 'j')
            or delta < 0 and (-delta .. 'k')
            or ''
        vim.cmd('normal! ' .. motion .. count .. 'z' .. dir)
    end
end

return M
