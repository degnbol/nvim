#!/usr/bin/env lua
return {
    -- recommended compiled fuzzy finder for telescope. Cannot be opt=true when needed by tzachar/cmp-fuzzy-path
    {
        "nvim-telescope/telescope-fzf-native.nvim", build='make',
        lazy = true, -- load as (fake) dependency of telescope
        config=function ()
            -- load the native fzf as recommended
            -- require'telescope'.load_extension('fzf')
            -- pcall is protected call, i.e. doesn't make a big deal out of errors
            -- taken from https://github.com/nvim-lua/kickstart.nvim/blob/master/init.lua
            pcall(require('telescope').load_extension, 'fzf')
        end
    },
    {
        "nvim-telescope/telescope.nvim",
        cmd = "Telescope",
        dependencies={
            "nvim-lua/plenary.nvim",
            "nvim-telescope/telescope-fzf-native.nvim",
        },
        init = function ()
            vim.keymap.set("n", "<leader>fb", "<Cmd>Telescope buffers<CR>", { desc="Buffer" })
            vim.keymap.set("n", "<leader>ff", "<Cmd>Telescope find_files<CR>", { desc="File" })
            vim.keymap.set("n", "<leader>fK", "<Cmd>Telescope help_tags<CR>", { desc="Help tags" })
            vim.keymap.set("n", "<leader>fh", "<Cmd>Telescope highlights<CR>", { desc="Highlights" })
            vim.keymap.set("n", "<leader>fj", "<Cmd>Telescope jumplist<CR>", { desc="Highlights" })
            vim.keymap.set("n", "<leader>fo", "<Cmd>Telescope oldfiles<CR>", { desc="Recent" })
            vim.keymap.set("n", "<leader>fr", "<Cmd>Telescope resume<CR>", { desc="Resume" })
            vim.keymap.set("n", "<leader>fw", "<Cmd>Telescope live_grep<CR>", { desc="Grep" })
            vim.keymap.set("n", "<leader>fW", "<Cmd>Telescope grep_string<CR>", { desc="Grep under cursor" })
            vim.keymap.set("n", "<C-/>", function ()
                require'telescope.builtin'.keymaps(require'telescope.themes'.get_ivy {
                    show_plug = false, -- simple solution to reduce width of the lhs column
                })
            end, { desc="Keymaps" })
        end,
        opts={
            defaults = {
                layout_strategy = "vertical",
                prompt_prefix = '',
                selection_caret = '',
                entry_prefix = '',
                multi_icon = '',
            },
            pickers = {
                colorscheme = {
                    -- window has to be big enough to show the preview window
                    enable_preview = true,
                }
            },
            extensions = {
                dash = {
                    file_type_keywords = {
                        python = {"python", "numpy", "scipy", "pandas"}
                    }
                }
            }
        }
    },
    {
        -- enables command :Telescope symbols
        "nvim-telescope/telescope-symbols.nvim",
        enabled = false,
    }, 
    {'sudormrfbin/cheatsheet.nvim', init = function ()
        vim.keymap.set("n", "<leader>?", "<Cmd>Cheatsheet<CR>", {desc="Cheatsheet"})
    end, cmd="Cheatsheet", dependencies={'nvim-telescope/telescope.nvim'}},
    -- TODO we don't actually use this yet
    {"lalitmee/browse.nvim", enabled = false, dependencies={"nvim-telescope/telescope.nvim"}}, -- search stackoverflow quicker
    {"nvim-telescope/telescope-bibtex.nvim",
    dependencies={'nvim-telescope/telescope.nvim'}, ft="tex",
    config=function()
        -- https://github.com/nvim-telescope/telescope-bibtex.nvim
        telescope = require "telescope"
        telescope.load_extension("bibtex")
        -- the following is used to detect *.bib file used
        -- by looking for \bibliography and \addbibresource.
        telescope.setup { context = true, }
    end},
    -- consider papis as well.
    -- bibliography references, mostly relevant for citations in .tex documents.
    -- { "jghauser/papis.nvim", 
    --     dependencies = { "kkharji/sqlite.lua", "nvim-lua/plenary.nvim", "MunifTanjim/nui.nvim", "nvim-treesitter/nvim-treesitter", "nvim-telescope/telescope.nvim", "hrsh7th/nvim-cmp" },
    --     rocks="lyaml", config = function() require("papis").setup() end,
    -- },
    {
        "debugloop/telescope-undo.nvim",
        dependencies = { -- note how they're inverted to above example
            {
                "nvim-telescope/telescope.nvim",
                dependencies = { "nvim-lua/plenary.nvim" },
            },
        },
        keys = {
            { -- lazy style key map
                "<leader>fu",
                "<cmd>Telescope undo<CR>",
                desc = "undo history",
            },
        },
        opts = {
            -- don't use `defaults = { }` here, do this in the main telescope spec
            extensions = {
                undo = {
                    -- telescope-undo.nvim config, see below
                },
                -- no other extensions here, they can have their own spec too
            },
        },
        config = function(_, opts)
            -- Calling telescope's setup from multiple specs does not hurt, it will happily merge the
            -- configs for us. We won't use data, as everything is in it's own namespace (telescope
            -- defaults, as well as each extension).
            require("telescope").setup(opts)
            require("telescope").load_extension("undo")

            vim.keymap.set("n", "<leader>fu", "<cmd>Telescope undo<CR>")
        end,
    },
}
