#!/usr/bin/env lua
return {
    -- recommended compiled fuzzy finder for telescope. Cannot be opt=true when needed by tzachar/cmp-fuzzy-path
    {
        "nvim-telescope/telescope-fzf-native.nvim", build='make', config=function ()
            -- load the native fzf as recommended
            -- require'telescope'.load_extension('fzf')
            -- pcall is protected call, i.e. doesn't make a big deal out of errors
            -- taken from https://github.com/nvim-lua/kickstart.nvim/blob/master/init.lua
            pcall(require('telescope').load_extension, 'fzf')
        end
    },
    {
        "nvim-telescope/telescope.nvim", dependencies={"nvim-lua/plenary.nvim", "nvim-telescope/telescope-fzf-native.nvim"},
        opts={
            defaults = {
                layout_strategy = "vertical",
                prompt_prefix = '',
                selection_caret = '',
                entry_prefix = '',
                multi_icon = '',
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
    {'sudormrfbin/cheatsheet.nvim', dependencies={'nvim-telescope/telescope.nvim'}}, -- <leader>? to give cheatsheet popup. 
    {"lalitmee/browse.nvim", dependencies={"nvim-telescope/telescope.nvim"}}, -- search stackoverflow quicker
    {"nvim-telescope/telescope-bibtex.nvim", dependencies={'nvim-telescope/telescope.nvim'}, config=function()
        -- https://github.com/nvim-telescope/telescope-bibtex.nvim
        telescope = require "telescope"
        telescope.load_extension("bibtex")
        -- the following is used to detect *.bib file used
        -- by looking for \bibliography and \addbibresource.
        telescope.setup {
            context = true,
        }
    end},
    -- consider papis as well.
    -- bibliography references, mostly relevant for citations in .tex documents.
    -- { "jghauser/papis.nvim", 
    --     dependencies = { "kkharji/sqlite.lua", "nvim-lua/plenary.nvim", "MunifTanjim/nui.nvim", "nvim-treesitter/nvim-treesitter", "nvim-telescope/telescope.nvim", "hrsh7th/nvim-cmp" },
    --     rocks="lyaml", config = function() require("papis").setup() end,
    -- },
}
