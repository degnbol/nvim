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
    indent = {
        -- when treesitter runs for julia the indent is wrong after inline for loop,
        -- but setting the indentexpr to the one from treesitter it works!
        -- The problem with indent happening when writing a bracket is actually due to
        -- indentkeys which are a set of keys that when written will trigger autoindent
        -- of a line. Setting the indent of a line is ALWAYS done by calling
        -- indentexpr which is GetJuliaIndent() by default. 
        -- Setting this enable flag sets indentexpr=nvim_treesitter#indent()
        enable = true
    }
}

vim.wo.foldlevel = 99 -- so we don't fold from the start
-- Fallback if treesitter folding doesn't work:
-- vim.wo.foldmethod = 'indent'
vim.wo.foldmethod = 'expr'
vim.wo.foldexpr = 'nvim_treesitter#foldexpr()'
