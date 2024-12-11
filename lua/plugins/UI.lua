#!/usr/bin/env lua
return {
    -- popular plugin to take advantage of some exposed UI hooks
    {
        "stevearc/dressing.nvim",
        -- enabled = false, -- was disabled since rename didn't seem to take effect.
        event = "VeryLazy",
        opts = {},
    },
    -- :Capture hi to call :hi where you can search etc.
    {"tyru/capture.vim", cmd="Capture"},
    {
        dir = vim.opt.runtimepath:get()[1] .. "/kittyREPL.nvim", dev=true, opts={
            keymap={
                focus="<C-CR>",
                setlast="<leader>rr",
                new="<leader>rs",
                run="<CR>",
                paste="<S-CR>",
                help="<leader>K",
                runLine="<CR><CR>",
                runLineFor="<leader>ro",
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
