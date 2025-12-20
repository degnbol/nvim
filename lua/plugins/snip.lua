local map = require "utils/keymap"

return {
    {
        'L3MON4D3/LuaSnip',
        lazy = true,                                      -- load as completion dependency
        dependencies = "nvim-treesitter/nvim-treesitter", -- depend on treesitter for the ft_func
        -- the make command is optional: https://github.com/L3MON4D3/LuaSnip
        build = "make install_jsregexp",
    },
    -- these default snippets can be replaced with my custom snippets when I have enough
    {
        "honza/vim-snippets",
        enabled = true,
        lazy = true, -- load as cmp dependency
        -- dependencies = { 'saadparwaiz1/cmp_luasnip' },
        config = function()
            require("luasnip.loaders.from_snipmate").lazy_load { exclude = { "tex", "julia", "all", "_", "python" } }
        end
    },
    {
        "rafamadriz/friendly-snippets",
        enabled = true,
        lazy = true, -- load as cmp dependency
        -- dependencies = { 'saadparwaiz1/cmp_luasnip' },
        config = function()
            require("luasnip.loaders.from_vscode").lazy_load { exclude = { "tex", "julia", "license", "global", "all", "python" } }
        end
    },
    {
        "chrisgrieser/nvim-scissors",
        lazy = true,
        dependencies = {"stevearc/dressing.nvim", "ibhagwan/fzf-lua"},
        opts = {},
        init = function ()
            map.n(
                "<leader>xe",
                function() require("scissors").editSnippet() end,
                "Snippet: Edit"
            )
            -- when used in visual mode, prefills the selection as snippet body
            map.nx(
                "<leader>xa",
                function() require("scissors").addNewSnippet() end,
                "Snippet: Add"
            )
        end,
    }
}
