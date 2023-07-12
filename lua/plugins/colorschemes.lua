#!/usr/bin/env lua
return {
    {
        "norcalli/nvim-base16.lua",
        -- TODO: find a way to set these schemes with a cmd, ideally the colorscheme cmd
        lazy = true,
        dependencies={"norcalli/nvim.lua"},
        config = function ()
            base16 = require'base16'
            base16.themes["gigavolt"] = require("base16/gigavolt")
            base16.themes["decaf"] = require("base16/decaf")
            -- local theme_names = base16.theme_names()
            -- base16(base16.themes.decaf)
        end,
    },

    {
        "NTBBloodbath/sweetie.nvim",
        cmd = "Colorscheme",
        opts = {
            overrides = {
                Comment = { italic = false },
                CommentBold = { italic = false, },
            },
        },
    },

    {
        'maxmx03/fluoromachine.nvim',
        lazy=true,
        opts = {
            glow = true,
            -- theme = 'fluoromachine',
            -- theme = 'retrowave',
            theme = 'delta',
            overrides = {
                 ['Comment'] = { italic = false },
             }
        }
    },

    {
        "rebelot/kanagawa.nvim",
        lazy=false,
        priority = 1000,
        config = function ()
            require'kanagawa'.setup {
                commentStyle = { italic = false },
            }
            vim.cmd "colorscheme kanagawa"
        end
    },

    {
        'rose-pine/neovim',
        lazy=true,
        name = 'rose-pine',
        opts = { }
    },

    {"xero/miasma.nvim", lazy=true},

    {
        "ribru17/bamboo.nvim",
        lazy = true,
        opts = {
            code_style = {
                comments = 'none',
                keywords = 'italic',
            },
        },
    },

    {
        "EdenEast/nightfox.nvim",
        lazy = true,
        opts = {
            options = { styles = { keywords = "italic" } }
        }
    },

    {
        "projekt0n/github-nvim-theme",
        lazy = true,
        config = function()
            require'github-theme'.setup { }
        end
    },

    {
        "catppuccin/nvim",
        lazy = true,
        name = "catppuccin",
    }
}
