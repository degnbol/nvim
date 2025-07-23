
return {
    -- popular plugin to take advantage of some exposed UI hooks
    {
        "stevearc/dressing.nvim",
        -- enabled = false, -- was disabled since rename didn't seem to take effect.
        event = "VeryLazy",
        opts = {},
    },
    -- :Capture hi to call :hi where you can search etc.
    { "tyru/capture.vim", cmd = "Capture" },
    {
        dir = vim.opt.runtimepath:get()[1] .. "/kittyREPL.nvim",
        dev = true,
        opts = {
            keymap = {
                focus = "<C-CR>",
                setlast = "<leader>rr",
                new = "<leader>rs",
                run = "<CR>",
                paste = "<S-CR>",
                help = "<leader>K",
                runLine = "<CR><CR>",
                runLineFor = "<leader>ro",
                runLineForI = "<leader>ri",
                pasteLine = "<S-CR><S-CR>",
                runVisual = "<CR>",
                pasteVisual = "<S-CR>",
                q = "<leader>rq",
                cr = "<leader>r<S-CR>",
                ctrld = "<leader>rd",
                ctrlc = "<leader>rc",
                interrupt = "<leader>rk",
                scrollStart = "[r",
                scrollUp = "[r",
                scrollDown = "]r",
                progress = "<leader>rp",
                editPaste = "<leader>re",
            },
            exclude = { tex = true, text = true, tsv = true, markdown = true },
            progress = true,
            editpaste = true,
            closepager = true,
        },
    },
    -- File explorer as a buffer with manipulation abilities.
    -- Differes from mini.files by only having a single view taking up the whole screen and has different default keymaps and preview behaviour.
    {
        "stevearc/quicker.nvim",
        keys = { "<leader>qq", "<leader>qo" },
        lazy = true,
        opts = {
            keys = {
                {
                    ">",
                    function()
                        require("quicker").expand({ before = 2, after = 2, add_to_existing = true })
                    end,
                    desc = "Expand quickfix context",
                },
                {
                    "<",
                    function()
                        require("quicker").collapse()
                    end,
                    desc = "Collapse quickfix context",
                },
            },
        },
    },
    {
        "stevearc/oil.nvim",
        lazy = true,
        dependencies = { "nvim-tree/nvim-web-devicons" },
        init = function()
            -- use <leader>e for explore followed by E since e is for mini.files
            -- vim.keymap.set("n", "<leader>eE", function () require"oil".open() end, { desc = "Oil" })
            -- similar to the go up one level keymap within oil
            vim.keymap.set("n", "<leader>-", function() require "oil".open() end, { desc = "Oil" })
        end,
        opts = {
            -- no prompt on rename, including folder change.
            -- The prompt is a nice list of the performed changes.
            -- skip_confirm_for_simple_edits = true,
            keymaps = {
                -- I regularly /-search for filename and don't need to hl that inside the buffer after opening
                ["<CR>"] = function()
                    vim.cmd "nohl"
                    require "oil".select()
                end
            }
        },
    },
    {
        "A7Lavinraj/fyler.nvim",
        dependencies = { "nvim-tree/nvim-web-devicons" },
        opts = {
            icon_provider = "nvim-web-devicons",
        },
    },
    {
        "jake-stewart/auto-cmdheight.nvim",
        lazy = false,
        opts = {
            -- max cmdheight before displaying hit enter prompt.
            max_lines = 5,

            -- number of seconds until the cmdheight can restore.
            duration = 2,

            -- whether key press is required to restore cmdheight.
            remove_on_key = true,

            -- always clear the cmdline after duration and key press.
            -- by default it will only happen when cmdheight changed.
            clear_always = false,
        }
    },
}
