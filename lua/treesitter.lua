require"nvim-treesitter.configs".setup {
    ensure_installed = {
        "bash",
        "lua",
        "json",
        "python",
        "julia",
        "latex",
        "java",
        "r"
    },
    highlight = {
        enable = true,
        use_languagetree = true
    },
}


-- vim.wo.foldmethod = 'indent'
-- vim.wo.foldlevel = 99 -- so we don't fold from the start
-- I would like to use treesitter folding but treesitter doesn't work if I write the two lines as instructed on their readme
-- vim.wo.foldmethod = 'expr'
-- vim.wo.foldexpr = 'nvim_treesitter#foldexpr()'
