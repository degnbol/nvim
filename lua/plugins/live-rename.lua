-- Place the buffer cursor at the caret's in-word offset within the NEW name
-- after a rename, matching where it sat in live-rename's float — the behaviour
-- the old wordinput.nvim had. live-rename's own submit leaves the cursor at its
-- original absolute column (apply_text_edits preserves it via extmark gravity),
-- which loses the offset when the name length changes.
--
-- The float caret's byte offset is the target-within-word. live-rename applies
-- the edit synchronously (firing on_bytes on the doc buffer) *before* closing
-- the float, so the first on_bytes while the float is still valid is the rename
-- itself, and the float caret is still readable then. A later edit fires with
-- the float already gone — that's the disarm signal, so no stale attach hijacks
-- an unrelated edit.
local function rename_keeping_cursor()
    local win = vim.api.nvim_get_current_win()
    local buf = vim.api.nvim_win_get_buf(win)
    local row, col = unpack(vim.api.nvim_win_get_cursor(win))
    -- Recover the word start: \k respects buffer-local 'iskeyword', so match in
    -- the buffer's context (mirrors <cword>, which is what rename prefills).
    local before = vim.api.nvim_buf_get_lines(buf, row - 1, row, false)[1]:sub(1, col)
    local offset = #vim.api.nvim_buf_call(buf, function()
        return vim.fn.matchstr(before, [[\k*$]])
    end)
    local word_start = col - offset

    require("live-rename").rename()

    -- rename() focuses its float; if it bailed (no client / invalid position)
    -- the window is unchanged and there's nothing to track.
    local float_win = vim.api.nvim_get_current_win()
    if float_win == win then return end

    vim.api.nvim_buf_attach(buf, false, {
        on_bytes = function()
            if not vim.api.nvim_win_is_valid(float_win) then return true end
            local caret = vim.api.nvim_win_get_cursor(float_win)[2]
            vim.schedule(function()
                pcall(vim.api.nvim_win_set_cursor, win, { row, word_start + caret })
            end)
            return true
        end,
    })
end

return {
    -- Live LSP rename in a float overlaid on the word in place, with the line
    -- reflowing as you type. Replaced our own modules/wordinput.nvim, which did
    -- the same overlay but not the reflow.
    {
        "live-rename.nvim",
        keys = {
            { "grn", rename_keeping_cursor, desc = "Rename" },
        },
    },
}
