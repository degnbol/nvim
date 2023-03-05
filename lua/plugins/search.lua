#!/usr/bin/env lua
return {
    -- let's search result box show number of matches when there's >99 matches
    "google/vim-searchindex",
    -- show a scrollbar, mostly in order to show search results far away in file
    -- requires hlslens to show search results in scrollbar
    {"petertriho/nvim-scrollbar", dependencies={"kevinhwang91/nvim-hlslens"}, opts = {
        show = false, -- enable with :ScrollbarToggle etc.
        marks = { Search = { color = "yellow" } },
        handlers = {
            diagnostics = false,
            search = true
        }
    }},
    -- show counter for how many n or N's a search result is away from cursor
    -- also lots of other search highlight customizations possible
    {"kevinhwang91/nvim-hlslens", config=function()
        local api = vim.api

        local config = require('hlslens.config')
        local render = require('hlslens.render')
        local extmark = require('hlslens.render.extmark')

        local bufs = {}
        local ns = api.nvim_create_namespace('hlslens')

        -- switch if you like or don't like the virtual text when searching
        if true then
            config.override_lens = function(render, posList, nearest, idx, relIdx) end
        else
            config.override_lens = function(render, posList, nearest, idx, relIdx)
                local sfw = vim.v.searchforward == 1
                local indicator, text, chunks
                local absRelIdx = math.abs(relIdx)
                if absRelIdx > 1 then
                    indicator = ('%d%s'):format(absRelIdx, sfw ~= (relIdx > 1) and '▲' or '▼')
                elseif absRelIdx == 1 then
                    indicator = sfw ~= (relIdx == 1) and '▲' or '▼'
                else
                    indicator = ''
                end

                local lnum, col = unpack(posList[idx])
                if nearest and indicator == '' then
                    text = ''
                else
                    text = indicator
                end

                bufnr = api.nvim_get_current_buf()
                render.setVirt(bufnr, lnum, col, text, nearest)
            end

            render.setVirt = function(bufnr, lnum, col, text, nearest)
                if nearest then
                    hlgroup = 'HlSearchLensNear'
                else
                    hlgroup = 'HlSearchLens'
                end

                bufs[bufnr] = true
                api.nvim_buf_set_extmark(bufnr, ns, lnum-1, col-1, {
                    id = nil,
                    virt_text = {{text, hlgroup}},
                    virt_text_pos = 'overlay',
                    priority = 100,
                })
            end

            -- just copy-pasted to access bufs
            extmark.clearBuf = function(bufnr)
                if not bufnr then
                    return
                end
                bufnr = bufnr == 0 and api.nvim_get_current_buf() or bufnr
                if bufs[bufnr] then
                    if api.nvim_buf_is_valid(bufnr) then
                        api.nvim_buf_clear_namespace(bufnr, ns, 0, -1)
                        -- seemly a bug for nvim_create_namespace, can't clear extmarks totally
                        for _, extm in pairs(api.nvim_buf_get_extmarks(bufnr, ns, 0, -1, {})) do
                            api.nvim_buf_del_extmark(bufnr, ns, extm[1])
                        end
                    end
                    bufs[bufnr] = nil
                end
            end

            extmark.clearAllBuf = function()
                for bufnr in pairs(bufs) do
                    extmark.clearBuf(bufnr)
                end
                bufs = {}
            end

        end

        -- g/ should highlight last search (using the google search plugin maybe)
        -- but it gets broken by this plugin so we mimmic its behavior:
        local utils = require('utils')
        utils.map("n", "g/", "//<CR>`'", {noremap=true, silent=true})


        -- integrate with haya14busa/vim-asterisk
        vim.api.nvim_set_keymap('n', '*', [[<Plug>(asterisk-z*)<Cmd>lua require('hlslens').start()<CR>]], {})
        vim.api.nvim_set_keymap('n', '#', [[<Plug>(asterisk-z#)<Cmd>lua require('hlslens').start()<CR>]], {})
        vim.api.nvim_set_keymap('n', 'g*', [[<Plug>(asterisk-gz*)<Cmd>lua require('hlslens').start()<CR>]], {})
        vim.api.nvim_set_keymap('n', 'g#', [[<Plug>(asterisk-gz#)<Cmd>lua require('hlslens').start()<CR>]], {})
        vim.api.nvim_set_keymap('x', '*', [[<Plug>(asterisk-z*)<Cmd>lua require('hlslens').start()<CR>]], {})
        vim.api.nvim_set_keymap('x', '#', [[<Plug>(asterisk-z#)<Cmd>lua require('hlslens').start()<CR>]], {})
        vim.api.nvim_set_keymap('x', 'g*', [[<Plug>(asterisk-gz*)<Cmd>lua require('hlslens').start()<CR>]], {})
        vim.api.nvim_set_keymap('x', 'g#', [[<Plug>(asterisk-gz#)<Cmd>lua require('hlslens').start()<CR>]], {})

        -- integreate with mg979/vim-visual-multi which is probably my intention to add later
        vim.cmd([[
        aug VMlens
        au!
        au User visual_multi_start lua require('vmlens').start()
        au User visual_multi_exit lua require('vmlens').exit()
        aug END
        ]])
    end},
}
