#!/usr/bin/env lua
return {
    -- popular plugin to take advantage of some exposed UI hooks
    -- doesn't seem to apply changes
    { "stevearc/dressing.nvim", enabled = false, event = "VeryLazy", opts = {} },
    -- :Capture hi to call :hi where you can search etc.
    {"tyru/capture.vim", cmd="Capture"},
    {
        dir = "$XDG_CONFIG_HOME/nvim/kittyREPL.nvim", dev=true, opts={
            keymap={
                focus="<C-CR>",
                setlast="<leader>rr",
                new="<leader>rs",
                run="<CR>",
                paste="<S-CR>",
                help="<leader>K",
                runLine="<CR><CR>",
                pasteLine="<S-CR><S-CR>",
                runVisual="<CR>",
                pasteVisual="<S-CR>",
                q="<leader>rq",
                cr="<leader>r<S-CR>",
                ctrld="<leader>rd",
                ctrlc="<leader>rc",
                interrupt="<leader>rk",
                scrollStart="[r",
                scrollUp="[r",
                scrollDown="]r",
                progress="<leader>rp",
                editPaste="<leader>re",
            },
            exclude = {tex=true, text=true, tsv=true, markdown=true},
            progress = true,
            editpaste = true,
            closepager = true,
        },
    },
    -- ultra fold
    {'kevinhwang91/nvim-ufo', dependencies={'kevinhwang91/promise-async'}, config=function()
        -- https://github.com/kevinhwang91/nvim-ufo

        -- 1 char width column showing icons to indicate folding levels
        -- vim.o.foldcolumn = '1'
        -- pretty icons maybe although default plus and minus are pretty clear
        -- vim.o.fillchars = [[eob: ,fold: ,foldopen:,foldsep: ,foldclose:]]
        -- vim.o.fillchars = [[eob: ,fold: ,foldopen:▾,foldsep: ,foldclose:▸]]
        -- vim.o.fillchars = [[eob: ,fold: ,foldopen:,foldsep: ,foldclose:]]
        -- vim.o.fillchars = [[eob: ,fold: ,foldopen:-,foldsep: ,foldclose:+]]

        vim.keymap.set('n', 'zR', require('ufo').openAllFolds, {desc="Open all folds"})
        vim.keymap.set('n', 'zM', require('ufo').closeAllFolds, {desc="Close all folds"})

        -- vim.api.nvim_set_hl(0, "UfoFoldedEllipsis", {link="Statement", default=true})
        -- vim.cmd "hi! default link UfoFoldedEllipsis Statement"

        require('ufo').setup {
            open_fold_hl_timeout = 0, --disable hl
            close_fold_kinds = {},
            fold_virt_text_handler = function(virtText, lnum, endLnum, width, truncate)
                local newVirtText = {}
                local suffix = (' … %d'):format(endLnum - lnum)
                local sufWidth = vim.fn.strdisplaywidth(suffix)
                local targetWidth = width - sufWidth
                local curWidth = 0
                for _, chunk in ipairs(virtText) do
                    local chunkText = chunk[1]
                    local chunkWidth = vim.fn.strdisplaywidth(chunkText)
                    if targetWidth > curWidth + chunkWidth then
                        table.insert(newVirtText, chunk)
                    else
                        chunkText = truncate(chunkText, targetWidth - curWidth)
                        local hlGroup = chunk[2]
                        table.insert(newVirtText, {chunkText, hlGroup})
                        chunkWidth = vim.fn.strdisplaywidth(chunkText)
                        -- str width returned from truncate() may less than 2nd argument, need padding
                        if curWidth + chunkWidth < targetWidth then
                            suffix = suffix .. (' '):rep(targetWidth - curWidth - chunkWidth)
                        end
                        break
                    end
                    curWidth = curWidth + chunkWidth
                end
                table.insert(newVirtText, {suffix, 'NonText'})
                return newVirtText
            end
        }
    end},
    -- file explorer as a buffer
    {"stevearc/oil.nvim",
    dependencies = { "nvim-tree/nvim-web-devicons" },
    config=function ()
        local oil = require "oil"
        oil.setup {
            -- no prompt on rename, including folder change.
            -- The prompt is a nice list of the performed changes.
            -- skip_confirm_for_simple_edits = true,
        }
        vim.keymap.set("n", "-", oil.open, { desc = "Open parent directory" })
        -- Move the builtin - to _ since - is used here above and it also is natural to use shift for both + and -
        -- -+ are different from jk since they go at start of line
        vim.keymap.set('n', "_", "-")
    end}
}
