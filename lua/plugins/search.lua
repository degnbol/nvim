local util = require "utils/init"
local map = require "utils/keymap"

---Horizontally scroll the current window so the whole match under the cursor
---is visible, its end flush against the right edge. Native search already keeps
---the match *start* (the cursor) on screen, so only a right-clipped end needs
---revealing; an already-visible match is left where native put it. No-op when
---'wrap' is on or there is no single-line match at the cursor.
local function reveal()
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

-- Shared reveal trigger. Defined at module top level (not in a plugin's `after`)
-- so both the / commit chain and vim-asterisk can reference it regardless of
-- load order.
map.n('<Plug>(RevealMatch)', reveal)

---Repeat search with `key` (n/N): jump natively, keep vim-searchindex's match
---counter, then scroll to reveal the whole match. Factory so n and N share it.
local function search_repeat(key)
    return function()
        local ok, err = pcall(vim.cmd.normal, { vim.v.count1 .. key, bang = true })
        if not ok then
            vim.api.nvim_echo({ { (err:gsub('^Vim:', '')), 'ErrorMsg' } }, true, {})
            return
        end
        -- <Plug>SearchIndex restores a view snapshot taken with the match at
        -- the left edge, so it must run (synchronously, hence 'x') before the
        -- reveal or it clobbers the anchored leftcol.
        vim.api.nvim_feedkeys(vim.keycode('<Plug>SearchIndex'), 'mx', false)
        reveal()
    end
end

return {
    -- Let's search result box show number of matches when there's >99 matches.
    {
        "vim-searchindex",
        -- keys = { "/", "g/", "*", "#", "g*", "g#" }
        after = function()
            map.n('n', search_repeat('n'), "Search next")
            map.n('N', search_repeat('N'), "Search prev")
            -- Chain the reveal onto searchindex's own <CR> cmdline map (which
            -- appends <Plug>SearchIndex after a / or ? search) so the jump,
            -- count and reveal run in one synchronous pass, same as n/N.
            map.c('<CR>', function()
                if vim.fn.getcmdtype():match('[/?]') then
                    return '<CR><Plug>SearchIndex<Plug>(RevealMatch)'
                end
                return '<CR>'
            end, "Search: jump, count, reveal", { expr = true, remap = true })
        end,
    },
    -- Show a scrollbar, mostly in order to show search results far away in file.
    -- Requires hlslens to show search results in scrollbar.
    {
        "nvim-scrollbar",
        cmd = { "ScrollbarToggle", "ScrollbarShow" },
        before = function()
            require("lz.n").trigger_load("nvim-hlslens")
            map.n("yoS", function()
                require "scrollbar.utils".toggle()
                -- The solution from scrollbar to not auto start hlslens virtual texts doesn't seem to work.
                require "hlslens".disable()
            end, "Scrollbar")
        end,
        after = function()
            require("scrollbar").setup {
                show = false, -- enable with :ScrollbarToggle etc.
                handlers = {
                    diagnostics = true,
                    search = true
                }
            }
        end,
    },
    -- Improvements to nmap *, etc.
    {
        "vim-asterisk",
        before = function()
            -- Option to keep cursor at the same location within searched cword.
            -- Useful for refactoring, e.g. FOO_|BAR.
            vim.g["asterisk#keeppos"] = 1
            -- Version of normal map * that doesn't move cursor is <Plug>(asterisk-z*).
            -- Instead of having mappings with z prefix, I think a nice
            -- solution is * searched as default vim on all repeated presses.
            local function is_searching_for_cword()
                local cword = vim.fn.expand("<cword>")
                -- Check if cword is the most recent search term.
                -- Not currently checking if we are currently searching or if this was just from the last search.
                local search_term = vim.fn.getreg("/")
                return cword == search_term or "\\<" .. cword .. "\\>" == search_term
            end
            ---Get the keymap expression from vim-asterisk given the stock vim expression.
            ---Returns function since the string expression will depend on if
            ---we have started a cword search or not.
            ---@param stock_expr string *, #, g*, or g#
            ---@return function
            local function asterisk_expr(stock_expr)
                local prefix, suffix
                if #stock_expr > 1 then
                    prefix = stock_expr:sub(1, -2)
                    suffix = stock_expr:sub(-1)
                else
                    prefix = ""
                    suffix = stock_expr
                end
                return function()
                    local z = is_searching_for_cword() and "" or "z"
                    local keys = "<Plug>(asterisk-" .. prefix .. z .. suffix .. ")"
                    -- Reveal only in normal mode; in o-pending/visual * is a
                    -- motion and a trailing reveal would corrupt the operator.
                    if vim.fn.mode() == "n" then
                        keys = keys .. "<Plug>(RevealMatch)"
                    end
                    return keys
                end
            end
            map.nox("*", asterisk_expr("*"), "Search next \\<cword\\>", { expr = true })
            map.nox("#", asterisk_expr("#"), "Search previous \\<cword\\>", { expr = true })
            map.nox("g*", asterisk_expr("g*"), "Search next cword", { expr = true })
            map.nox("g#", asterisk_expr("g#"), "Search previous cword", { expr = true })
        end,
    },
    -- Show counter for how many n or N's a search result is away from cursor.
    {
        "nvim-hlslens",
        lazy = true,
        cmd = { "HlSearchLensEnable", "HlSearchLensToggle" },
        before = function()
            map.n("yoH", function() require 'hlslens'.toggle() end, "HlSearchLens", { silent = true })
        end,
        after = function()
            require("hlslens").setup({
                auto_enable = false,
                -- :nohlsearch when moving out of search term.
                -- cons: disables when scrolling far since cursor moves. Doesn't disable right away on insert, only after a change.
                -- better solution made manually in ModeChanged and specific keypresses, such as Esc.
                -- calm_down = true,
            })
        end,
    },
}
