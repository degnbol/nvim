#!/usr/bin/env lua

vim.keymap.set("n", "<leader>cs", function()
    -- Load colorschemes. Using lazy keys didn't work
    -- https://www.reddit.com/r/neovim/comments/12tcx0b/attempt_at_adding_color_schemes_to_list_of/
    vim.api.nvim_exec_autocmds("User", { pattern = "ColorSchemeLoad" })
    -- simple attempt at only showing dark themes if we are in dark-mode
    if vim.o.background == "dark" then
        vim.api.nvim_exec_autocmds("User", { pattern = "ColorSchemeLoadDark" })
    end
    require("telescope.builtin").colorscheme()
end, { desc="Colorscheme" })

return {
    {
        "norcalli/nvim-base16.lua",
        -- TODO: find way to use each of the builtin ones from this with colorscheme cmd
        lazy = true,
        -- event = "User ColorSchemeLoad",
        dependencies={"norcalli/nvim.lua"},
    },

    {
        "NTBBloodbath/sweetie.nvim",
        event = "User ColorSchemeLoad",
        config = function ()
            vim.g.sweetie = { overrides = {
                Comment = { italic = false },
                CommentBold = { italic = false, },
            }}
        end,
    },

    {
        'maxmx03/fluoromachine.nvim',
        event = "User ColorSchemeLoadDark",
    },

    {
        "rebelot/kanagawa.nvim",
        event = "User ColorSchemeLoad",
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
        event = "User ColorSchemeLoad",
    },

    {
        -- rose-pine, rose-pine-main, rose-pine-moon, rose-pine-dawn
        'rose-pine/neovim',
        event = "User ColorSchemeLoad",
        name = 'rose-pine',
        opts = { }
    },

    {
        "xero/miasma.nvim",
        event = "User ColorSchemeLoadDark",
    },

    {
        "ribru17/bamboo.nvim",
        event = "User ColorSchemeLoadDark",
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
        event = "User ColorSchemeLoad",
        opts = {
            options = { styles = { keywords = "italic" } }
        }
    },

    {
        "projekt0n/github-nvim-theme",
        name = "github-theme",
        event = "User ColorSchemeLoad",
    },

    {
        "catppuccin/nvim",
        event = "User ColorSchemeLoad",
        name = "catppuccin",
    },
}
