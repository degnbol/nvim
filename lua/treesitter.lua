require"nvim-treesitter.configs".setup {
    ensure_installed = {
        "bash",
        "lua",
        "json",
        "python",
        "julia",
        "latex",
        -- "java",
        -- "kotlin",
        "help", -- vim help files https://github.com/neovim/tree-sitter-vimdoc
        "r",
        "markdown",
        "toml",
        "vim",
        "regex",
        "make",
        "cmake",
        "cpp",
        "bibtex",
        "gitignore",
        "gitattributes",
        "diff", -- for diff output https://github.com/the-mikedavis/tree-sitter-diff
        "scheme", -- what treesitter queries (*.scm) are written in
    },
    highlight = {
        enable = true,
        disable = {
            "vim", -- not perfect
            -- "help", -- removes useful colors from :h group-name
            "markdown", -- not as good as regular syntax highlight
        }, 
        additional_vim_regex_highlighting = {
            "julia", -- basic things like true and false are not recognized as bool and I couldn't fix it with a custom highlights.scm
            "help", -- treesitter version removes useful colors from :h group-name
            "bash", -- spending too much time writing treesitter query. Also covers zsh.
        },
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
    -- opt-in to using treesitter for https://github.com/andymass/vim-matchup
    matchup = {
        enable = true,
        disable = {}, -- optional, list of language that will be disabled
    }
}

vim.wo.foldlevel = 99 -- so we don't fold from the start
-- Fallback if treesitter folding doesn't work:
-- vim.wo.foldmethod = 'indent'
vim.wo.foldmethod = 'expr'
vim.wo.foldexpr = 'nvim_treesitter#foldexpr()'

-- use bash treesitter for zsh since zsh is very basic
local ft_to_parser = require"nvim-treesitter.parsers".filetype_to_parsername
ft_to_parser.zsh = "bash"




