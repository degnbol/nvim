--- A focused float for editing the term under the cursor.
---
--- `require("wordinput").open(opts, on_confirm)` prefills `opts.default` and
--- opens in normal mode with the cursor at the same char offset within the term,
--- anchoring the float to the word start rather than the cursor. Assumes
--- `opts.default` is the term at the parent cursor (true for LSP rename); derives
--- offset/anchor from the live cursor. Uses a plain scratch buffer, so mode and
--- cursor are owned outright: no prompt buffer, no `startinsert!`, no auto-expand
--- reconfigure.
local M = {}

--- Byte offset of the cursor within the word (per 'iskeyword') it sits on, in
--- window `win`; 0 when the cursor is not on a word character. Mirrors <cword>,
--- which is what LSP rename prefills its input with.
--- @param win number
--- @return number offset
local function cursor_word_offset(win)
    local buf = vim.api.nvim_win_get_buf(win)
    local row, col = unpack(vim.api.nvim_win_get_cursor(win))
    local before = vim.api.nvim_buf_get_lines(buf, row - 1, row, false)[1]:sub(1, col)
    -- \k respects 'iskeyword', which is buffer-local, so match in `buf`'s context.
    return #vim.api.nvim_buf_call(buf, function()
        return vim.fn.matchstr(before, [[\k*$]])
    end)
end

--- @param opts { default?: string }|nil
--- @param on_confirm fun(value: string|nil)
function M.open(opts, on_confirm)
    assert(type(on_confirm) == "function", "`on_confirm` must be a function")
    opts = opts or {}

    -- Capture origin before any window work: parent is unambiguously current.
    local offset = cursor_word_offset(vim.api.nvim_get_current_win())
    local default = opts.default or ""

    local buf = vim.api.nvim_create_buf(false, true)
    vim.api.nvim_buf_set_lines(buf, 0, -1, false, { default })

    local win = vim.api.nvim_open_win(buf, true, {
        relative = "cursor",
        -- Anchor the buffer's first column to the word start: offset cells left
        -- of the cursor, minus 1 more for the border column, so the prefilled
        -- term overlays the original word.
        col = -offset - 1,
        row = -1,
        width = math.max(30, #default + 3),
        height = 1,
        style = "minimal",
        border = "rounded",
    })
    -- Reuse snacks' input highlight groups so the look matches the generic
    -- vim.ui.input float (snacks defines these; it's always loaded).
    vim.wo[win].winhighlight =
        "NormalFloat:SnacksInputNormal,FloatBorder:SnacksInputBorder,FloatTitle:SnacksInputTitle"

    -- Normal mode, cursor at the char offset within the term. Plain buffer, so
    -- this is one synchronous call — nothing resets it.
    vim.api.nvim_win_set_cursor(win, { 1, offset })

    local done = false
    --- Close the float, then hand the result to on_confirm. Closing first means
    --- nvim_get_current_win() inside the ensuing rename resolves to the parent.
    --- @param value string|nil
    local function finish(value)
        if done then return end
        done = true
        if vim.api.nvim_win_is_valid(win) then
            pcall(vim.api.nvim_win_close, win, true)
        end
        on_confirm(value)
    end

    -- Defer buffer deletion past the closing keystroke (see neovim skill,
    -- references/ui.md: wiping a float buffer mid-keymap lets the key fall
    -- through to the parent window).
    vim.api.nvim_create_autocmd("WinClosed", {
        pattern = tostring(win),
        once = true,
        callback = function()
            vim.schedule(function()
                if vim.api.nvim_buf_is_valid(buf) then
                    vim.api.nvim_buf_delete(buf, { force = true })
                end
            end)
        end,
    })
    -- Focus loss is an implicit cancel, so a stray click doesn't strand the float.
    vim.api.nvim_create_autocmd("WinLeave", {
        buffer = buf,
        once = true,
        callback = function() finish(nil) end,
    })

    local function keymap(mode, lhs, fn)
        vim.keymap.set(mode, lhs, fn, { buffer = buf, nowait = true, silent = true })
    end
    -- <CR> in both normal and insert (an insert-mode <CR> would otherwise split
    -- the one-line buffer) → confirm the current line.
    keymap({ "n", "i" }, "<CR>", function()
        finish(vim.api.nvim_buf_get_lines(buf, 0, 1, false)[1])
    end)
    keymap({ "n", "i" }, "<C-c>", function() finish(nil) end)
    -- <Esc> cancels from normal mode; from insert it falls through to leave insert.
    keymap("n", "<Esc>", function() finish(nil) end)
end

return M
