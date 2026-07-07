--- An alternative to `vim.ui.input` for the term under the cursor.
---
--- Differences to stock `vim.ui.input`:
--- - Doesn't change the mode. In practise this mean opening in normal mode instead of insert mode.
--- - Placed such that the word appears in the float exactly where it was 
---   visible in the buffer, instead of placing the float relative to the cursor.
---
--- Similar differences to plugin alternatives such as snacks' input.

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

--- Float width that fits `text`: its display width plus one cell for the insert
--- caret's room past the last char at EOL. That trailing cell is left
--- transparent (see M.open) so it shows the buffer beneath rather than a chrome
--- strip — only the term itself carries the float's bg_dim.
--- @param text string
--- @return number width
local function width_for(text)
    return vim.api.nvim_strwidth(text) + 1
end

--- `on_confirm` receives the confirmed text (or nil if cancelled) and the parent
--- buffer column matching the float cursor's offset within the word: after an
--- in-place replacement of the word (e.g. LSP rename), placing the cursor there
--- lands it at the same offset within the new term as it had in the float.
--- @param opts { default?: string }|nil
--- @param on_confirm fun(value: string|nil, col: number)
function M.open(opts, on_confirm)
    assert(type(on_confirm) == "function", "`on_confirm` must be a function")
    opts = opts or {}

    -- Capture origin before any window work: parent is unambiguously current.
    local parent_win = vim.api.nvim_get_current_win()
    local offset = cursor_word_offset(parent_win)
    local word_start = vim.api.nvim_win_get_cursor(parent_win)[2] - offset
    local default = opts.default or ""

    local buf = vim.api.nvim_create_buf(false, true)
    vim.api.nvim_buf_set_lines(buf, 0, -1, false, { default })

    local win = vim.api.nvim_open_win(buf, true, {
        relative = "cursor",
        -- Borderless: anchor the first column to the word start (offset cells
        -- left of the cursor) and sit on the cursor's own line, so the prefilled
        -- term overlays the original word in place.
        col = -offset,
        row = 0,
        width = width_for(default),
        height = 1,
        style = "minimal",
        border = "none",
    })
    -- `minimal` turns off the float's own cursorline, so it renders as a slim
    -- strip. winblend=100 makes the whole float transparent; the trailing caret
    -- cell (past the term) is left that way so the buffer shows through it
    -- instead of a chrome strip. The term itself is repainted opaque below.
    vim.wo[win].winblend = 100

    -- Restore opaque chrome over the term only: a highlight whose blend=0
    -- overrides winblend for its cells (see :h highlight-blend), copying an
    -- existing group's colours. Read now so it tracks the current colorscheme.
    local ns = vim.api.nvim_create_namespace("wordinput")
    local chrome = vim.api.nvim_get_hl(0, { name = "PmenuSel" })
    vim.api.nvim_set_hl(0, "WordInputTerm", { fg = chrome.fg, bg = chrome.bg, blend = 0 })

    -- Normal mode, cursor at the char offset within the term. Plain buffer, so
    -- this is one synchronous call — nothing resets it.
    vim.api.nvim_win_set_cursor(win, { 1, offset })

    -- Fit the float to the term. A width-only reconfigure keeps the word-start
    -- `col` and the buffer cursor (unlike a full reconfigure), so the float
    -- widens rightward from the anchor.
    local function set_width(text)
        vim.api.nvim_win_set_config(win, { width = width_for(text) })
    end

    -- Paint the opaque chrome over the term only: the extmark spans exactly the
    -- term, so the trailing caret cell stays transparent and the buffer shows
    -- through it. The insert caret at EOL therefore takes its colour from that
    -- transparent cell (kitty's `cursor none` reverse-video); leaving it as-is is
    -- accepted — the caret reading the buffer beneath is consistent with the
    -- float overlaying the buffer in place.
    local function paint()
        if not vim.api.nvim_win_is_valid(win) then return end
        local line = vim.api.nvim_buf_get_lines(buf, 0, 1, false)[1] or ""
        vim.api.nvim_buf_set_extmark(buf, ns, 0, 0, {
            id = 1,
            end_col = #line,
            hl_group = "WordInputTerm",
        })
    end
    set_width(default)
    paint()

    -- Widen *ahead* of the typed char: InsertCharPre fires before v:char lands,
    -- so the float already fits it when neovim recomputes scroll at insertion.
    -- Repaint waits for the ensuing TextChangedI, once v:char is in the buffer.
    vim.api.nvim_create_autocmd("InsertCharPre", {
        buffer = buf,
        callback = function()
            local line = vim.api.nvim_buf_get_lines(buf, 0, 1, false)[1] or ""
            set_width(line .. vim.v.char)
        end,
    })

    -- Authoritative resize+repaint for everything InsertCharPre misses:
    -- insert-mode deletion (shrink) and any normal-mode edit (paste, x, ciw…).
    vim.api.nvim_create_autocmd({ "TextChanged", "TextChangedI" }, {
        buffer = buf,
        callback = function()
            local line = vim.api.nvim_buf_get_lines(buf, 0, 1, false)[1] or ""
            set_width(line)
            paint()
        end,
    })

    local done = false
    --- Close the float, then hand the result to on_confirm. Closing first means
    --- nvim_get_current_win() inside the ensuing rename resolves to the parent.
    --- @param value string|nil
    local function finish(value)
        if done then return end
        done = true
        -- Read the float cursor offset before closing; map it to the parent
        -- column so the word start it anchors to is preserved after an edit.
        local col = word_start
        if vim.api.nvim_win_is_valid(win) then
            col = word_start + vim.api.nvim_win_get_cursor(win)[2]
            pcall(vim.api.nvim_win_close, win, true)
        end
        on_confirm(value, col)
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
