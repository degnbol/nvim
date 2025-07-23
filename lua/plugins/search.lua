

return {
    -- let's search result box show number of matches when there's >99 matches
    {
        "google/vim-searchindex",
        -- keys={"/", "g/"} -- keys doesn't seem to work
    },
    -- show a scrollbar, mostly in order to show search results far away in file
    -- requires hlslens to show search results in scrollbar
    {
        "petertriho/nvim-scrollbar",
        cmd = { "ScrollbarToggle", "ScrollbarShow" },
        dependencies = { "kevinhwang91/nvim-hlslens" },
        init = function()
            vim.keymap.set("n", "yoS", "<Cmd>ScrollbarToggle<CR>", { desc = "Scrollbar" })
        end,
        opts = {
            show = false, -- enable with :ScrollbarToggle etc.
            marks = { Search = { color = "yellow" } },
            handlers = {
                diagnostics = false,
                search = true
            }
        }
    },
    -- show counter for how many n or N's a search result is away from cursor
    -- also lots of other search highlight customizations possible
    {
        "kevinhwang91/nvim-hlslens",
        cmd = { "HlSearchLensEnable", "HlSearchLensToggle" },
        init = function()
            function nmap(lhs, rhs, desc, opts)
                opts = opts or {}
                opts.desc = desc
                vim.keymap.set('n', lhs, rhs, opts)
            end

            -- integrate with haya14busa/vim-asterisk
            nmap("yoH", function() require 'hlslens'.toggle() end, "HlSearchLens", { silent = true })
            nmap('*', [[<Plug>(asterisk-z*) <Cmd>lua require('hlslens').start()<CR>]], "Search word under cursor")
            nmap('#', [[<Plug>(asterisk-z#) <Cmd>lua require('hlslens').start()<CR>]], "Search word under cursor")
            nmap('g*', [[<Plug>(asterisk-gz*)<Cmd>lua require('hlslens').start()<CR>]], "Search word under cursor")
            nmap('g#', [[<Plug>(asterisk-gz#)<Cmd>lua require('hlslens').start()<CR>]], "Search word under cursor")
            nmap('*', [[<Plug>(asterisk-z*) <Cmd>lua require('hlslens').start()<CR>]], "Search word under cursor")
            nmap('#', [[<Plug>(asterisk-z#) <Cmd>lua require('hlslens').start()<CR>]], "Search word under cursor")
            nmap('g*', [[<Plug>(asterisk-gz*)<Cmd>lua require('hlslens').start()<CR>]], "Search word under cursor")
            nmap('g#', [[<Plug>(asterisk-gz#)<Cmd>lua require('hlslens').start()<CR>]], "Search word under cursor")
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
