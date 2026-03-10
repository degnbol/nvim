return {
    {
        'L3MON4D3/LuaSnip',
        lazy = true,                                      -- load as completion dependency
        dependencies = "nvim-treesitter/nvim-treesitter", -- depend on treesitter for the ft_func
        -- the make command is optional: https://github.com/L3MON4D3/LuaSnip
        build = "make install_jsregexp",
    },
}
