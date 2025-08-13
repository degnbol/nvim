local map = require "utils/keymap"

return {
    -- Let's search result box show number of matches when there's >99 matches.
    {
        "google/vim-searchindex",
        -- keys = { "/", "g/", "*", "#", "g*", "g#" }
    },
    -- Show a scrollbar, mostly in order to show search results far away in file.
    -- Requires hlslens to show search results in scrollbar.
    {
        "petertriho/nvim-scrollbar",
        cmd = { "ScrollbarToggle", "ScrollbarShow" },
        dependencies = { "kevinhwang91/nvim-hlslens" },
        init = function()
            map.n("yoS", function()
                require "scrollbar.utils".toggle()
                -- The solution from scrollbar to not auto start hlslens virtual texts doesn't seem to work.
                require "hlslens".disable()
            end, "Scrollbar")
        end,
        config = function()
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
        "haya14busa/vim-asterisk",
        init = function()
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
                    if is_searching_for_cword() then
                        return "<Plug>(asterisk-" .. prefix .. suffix .. ")"
                    else
                        return "<Plug>(asterisk-" .. prefix .. "z" .. suffix .. ")"
                    end
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
        "kevinhwang91/nvim-hlslens",
        lazy = true,
        cmd = { "HlSearchLensEnable", "HlSearchLensToggle" },
        init = function()
            map.n("yoH", function() require 'hlslens'.toggle() end, "HlSearchLens", { silent = true })
        end,
        opts = {
            auto_enable = false,
            -- :nohlsearch when moving out of search term.
            -- cons: disables when scrolling far since cursor moves. Doesn't disable right away on insert, only after a change.
            -- better solution made manually in ModeChanged and specific keypresses, such as Esc.
            -- calm_down = true,
        },
    },
}
