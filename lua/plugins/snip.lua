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
        enabled = false,
        lazy = true, -- load as cmp dependency
        -- dependencies = { 'saadparwaiz1/cmp_luasnip' },
        config = function()
            require("luasnip.loaders.from_snipmate").lazy_load { exclude = { "tex", "julia", "all", "_", "python" } }
        end
    },
    {
        "rafamadriz/friendly-snippets",
        enabled = false,
        lazy = true, -- load as cmp dependency
        -- dependencies = { 'saadparwaiz1/cmp_luasnip' },
        config = function()
            require("luasnip.loaders.from_vscode").lazy_load { exclude = { "tex", "julia", "license", "global", "all", "python" } }
        end
    },
}
