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
        use_languagetree = true,
        -- sometimes highlight fails with julia, so having both highlighters active helps.
        -- again, this affects indent somehow.
        additional_vim_regex_highlighting = "julia",
    },
    incremental_selection = {
        enable = true,
        keymaps = {
            init_selection = "<leader><up>",
            node_incremental = "<leader><up>",
            scope_incremental = "<leader><S-up>",
            node_decremental = "<leader><down>",
        },
    },
}

vim.wo.foldlevel = 99 -- so we don't fold from the start
-- Fallback if treesitter folding doesn't work:
-- vim.wo.foldmethod = 'indent'
vim.wo.foldmethod = 'expr'
vim.wo.foldexpr = 'nvim_treesitter#foldexpr()'
