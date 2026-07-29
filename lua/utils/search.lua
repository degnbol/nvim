local util = require "utils/init"

local M = {}

---Horizontally scroll the current window so the whole match under the cursor
---is visible. Native search only guarantees the match *start* (the cursor) is
---on screen, so a match reaching past the right edge stays clipped. No-op when
---'wrap' is on, when there is no single-line match at the cursor, or when the
---match already fits.
function M.reveal()
    if vim.wo.wrap then return end
    local pattern = vim.fn.getreg('/')
    if pattern == '' then return end
    local match_end = vim.fn.searchpos(pattern, 'cenz') -- {lnum, col}, no cursor move
    if match_end[1] ~= vim.fn.line('.') then return end -- no match, or spans lines

    local win = vim.fn.getwininfo(vim.api.nvim_get_current_win())[1]
    local text_width = win.width - win.textoff
    local leftcol = vim.fn.winsaveview().leftcol
    local end_col = vim.fn.virtcol({ match_end[1], match_end[2] })

    if end_col > leftcol + text_width then
        util.set_view({ leftcol = end_col - text_width })
    end
end

return M
