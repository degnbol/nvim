#!/usr/bin/env lua
return {
    {
        "norcalli/nvim-base16.lua",
        -- TODO: find way to use each of the builtin ones from this with colorscheme cmd
        lazy = true,
        -- event = "User ColorSchemeLoad", -- see whichkey
        dependencies={"norcalli/nvim.lua"},
    },

    {
        "NTBBloodbath/sweetie.nvim",
        event = "User ColorSchemeLoad", -- see whichkey
        opts = {
            overrides = {
                Comment = { italic = false },
                CommentBold = { italic = false, },
            },
        },
    },

    {
        'maxmx03/fluoromachine.nvim',
        -- simply loading this will set it. We define code in colors/ instead of the config here.
        event = "User ColorSchemeLoad", -- see whichkey
    },

    {
        "rebelot/kanagawa.nvim",
        event = "User ColorSchemeLoad", -- see whichkey
        config = function ()
            require'kanagawa'.setup {
                commentStyle = { italic = false },
            }
        end
    },

    {
        'rose-pine/neovim',
        event = "User ColorSchemeLoad", -- see whichkey
        name = 'rose-pine',
        opts = { }
    },

    {
        "xero/miasma.nvim",
        event = "User ColorSchemeLoad", -- see whichkey
    },

    {
        "ribru17/bamboo.nvim",
        event = "User ColorSchemeLoad", -- see whichkey
        opts = {
            code_style = {
                comments = 'none',
                keywords = 'italic',
            },
        },
    },

    {
        "EdenEast/nightfox.nvim",
        event = "User ColorSchemeLoad", -- see whichkey
        opts = {
            options = { styles = { keywords = "italic" } }
        }
    },

    {
        "projekt0n/github-nvim-theme",
        event = "User ColorSchemeLoad", -- see whichkey
        config = function()
            require'github-theme'.setup { }
        end
    },

    {
        "catppuccin/nvim",
        event = "User ColorSchemeLoad", -- see whichkey
        name = "catppuccin",
    }
}
