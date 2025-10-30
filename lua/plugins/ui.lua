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
	lazy = true,
	cmd = "Fyler",
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
    {
        'b0o/incline.nvim',
        init = function()
            -- This plugin only makes sense with global statusline,
            -- or no statusline at all.
            -- The latter is possible with a trick.
            -- By setting cmdheight=0 and (optionally) clearing the
            -- 'statusline' setting we can use a single line for cmdline and
            -- (secret) global statusline.
            -- The problem with this approach is how buggy cmdheight=0 still is.
            -- A different solution is laststatus=0 where the statusline is redesigned as a simple border.
            vim.opt.laststatus = 0
            -- TODO:
            -- decide if there should be more of standard statusline info (such as "HELP" for help files) that should be added to render.
            vim.opt.statusline = ("─"):rep(vim.api.nvim_win_get_width(0))
            vim.api.nvim_create_autocmd("WinResized", {
                group = vim.api.nvim_create_augroup("statusline-update", { clear = true }),
                callback = function()
                    vim.opt.statusline = ("─"):rep(vim.api.nvim_win_get_width(0))
                end
            })
        end,
        config = function()
            local helpers = require 'incline.helpers'
            local devicons = require 'nvim-web-devicons'
            require('incline').setup {
                ignore = {
                    -- Also display for help etc.
                    buftypes = {},
                    wintypes = {},
                    unlisted_buffers = false,
                },
                window = {
                    padding = 0,
                    margin = { horizontal = 0, vertical = 0, },
                    placement = {
                        -- Keep at top since it fits well with TS-context that also provides breadcrumps.
                        -- vertical = "bottom",
                    },
                    overlap = {
                        borders = false,
                        -- statusline = true,
                        -- tabline = true,
                        -- winbar = true,
                    },
                    -- Above fidget if margin.vertical =1 and placement.vertical = "bottom" but not if margin.vertical = 0
                    -- zindex = 100,
                },
                render = function(props)
                    local filename = vim.fn.fnamemodify(vim.api.nvim_buf_get_name(props.buf), ':t')
                    if filename == '' then
                        filename = '[No Name]'
                    end
                    local ft_icon, ft_color = devicons.get_icon_color(filename)
                    local modified = vim.bo[props.buf].modified
                    return {
                        ft_icon and { ft_icon, guifg = ft_color } or '',
                        ' ',
                        { filename, gui = modified and 'bold' or nil, guifg = props.focused and "white" or "gray" },
                    }
                end,
            }
        end,
        event = 'VeryLazy',
    },
}
