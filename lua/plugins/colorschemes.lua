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
        event = "User ColorSchemeLoadDark", -- see whichkey
    },

    {
        "rebelot/kanagawa.nvim",
        event = "User ColorSchemeLoad", -- see whichkey
        opts = {
            commentStyle = { italic = false },
            colors = {
                theme = {
                    lotus = {
                        -- too similar to comments
                        ui = { nontext = "#acc295", },
                    }
                }
            }
        }

    },

    {
        -- minimal colors
        -- zenwritten, neobones, vimbones, rosebones, forestbones, nordbones, 
        -- tokyobones, seoulbones, duckbones, zenburned, kanagawabones, 
        -- randombones
        -- kanagawa-bones is from here and currently is the default dark theme for kitty
        "mcchrish/zenbones.nvim",
        dependencies = { "rktjmp/lush.nvim" },
        event = "User ColorSchemeLoad", -- see whichkey
    },

    {
        -- rose-pine, rose-pine-main, rose-pine-moon, rose-pine-dawn
        'rose-pine/neovim',
        event = "User ColorSchemeLoad", -- see whichkey
        name = 'rose-pine',
        opts = { }
    },

    {
        "xero/miasma.nvim",
        event = "User ColorSchemeLoadDark", -- see whichkey
    },

    {
        "ribru17/bamboo.nvim",
        event = "User ColorSchemeLoadDark", -- see whichkey
        opts = {
            code_style = {
                comments = 'none',
                keywords = 'italic',
            },
        },
    },

    {
        -- nightfox, dayfox, dawnfox, duskfox, nordfox, terafox, carbonfox
        "EdenEast/nightfox.nvim",
        event = "User ColorSchemeLoad", -- see whichkey
        opts = {
            options = { styles = { keywords = "italic" } }
        }
    },

    {
        "projekt0n/github-nvim-theme",
        name = "github-theme",
        event = "User ColorSchemeLoad", -- see whichkey
    },

    {
        "catppuccin/nvim",
        event = "User ColorSchemeLoad", -- see whichkey
        name = "catppuccin",
    },
}
