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
    -- File explorer as a buffer with manipulation abilities.
    -- Differes from mini.files by only having a single view taking up the whole screen and has different default keymaps and preview behaviour.
    {
        "stevearc/oil.nvim",
        lazy = true,
        dependencies = { "nvim-tree/nvim-web-devicons" },
        init = function ()
            -- use <leader>e for explore followed by E since e is for mini.files
            -- vim.keymap.set("n", "<leader>eE", function () require"oil".open() end, { desc = "Oil" })
            -- similar to the go up one level keymap within oil
            vim.keymap.set("n", "<leader>-", function () require"oil".open() end, { desc = "Oil" })
        end,
        opts = {
            -- no prompt on rename, including folder change.
            -- The prompt is a nice list of the performed changes.
            -- skip_confirm_for_simple_edits = true,
            keymaps = {
                -- I regularly /-search for filename and don't need to hl that inside the buffer after opening
                ["<CR>"] = function ()
                    vim.cmd "nohl"
                    require"oil".select()
                end
            }
        },
    },
    {
        "sphamba/smear-cursor.nvim",
        enabled = false,
        opts = {                         -- Default  Range
            stiffness = 0.8,               -- 0.6      [0, 1]
            trailing_stiffness = 0.5,      -- 0.3      [0, 1]
            distance_stop_animating = 0.5, -- 0.1      > 0
            hide_target_hack = false,      -- true     boolean
            legacy_computing_symbols_support = true,
        },
    }
}
