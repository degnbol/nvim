#!/usr/bin/env lua
return {
    "ojroques/nvim-bufdel", -- :BufDel that deletes a buffer better than built-in :bdelete and :bwipeout, by preserving layout and closing terminal buffers better.
    -- {"moll/vim-bbye", config=function()
    --     vim.api.nvim_set_keymap("n", "<leader>x", [[<Cmd>Bdelete<CR>]], {silent = true})
    -- end},
    "tyru/capture.vim", -- :Capture hi to call :hi where you can search etc.
    {
        "degnbol/kittyREPL.nvim", opts={
            keymap={
                line="<CR><CR>",
                visual="<CR>",
                operator="<CR>",
                win="<leader><CR>",
            },
            exclude = {tex=true, text=true, tsv=true, markdown=true},
        }
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

        vim.keymap.set('n', 'zR', require('ufo').openAllFolds)
        vim.keymap.set('n', 'zM', require('ufo').closeAllFolds)

        require('ufo').setup()
    end},
    -- file explorer as a buffer
    {"stevearc/oil.nvim", config=function ()
        local oil = require "oil"
        oil.setup {
            -- no prompt on rename, including folder change.
            -- The prompt is a nice list of the performed changes.
            -- skip_confirm_for_simple_edits = true,
        }
        vim.keymap.set("n", "-", oil.open, { desc = "Open parent directory" })
    end}
}
