require"nvim-treesitter.configs".setup {
    ensure_installed = {
        "bash",
        "lua",
        "json",
        "python",
        "julia",
        "latex",
        -- "java",
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
local ft_to_lang = require('nvim-treesitter.parsers').ft_to_lang
require('nvim-treesitter.parsers').ft_to_lang = function(ft)
    if ft == 'zsh' then return 'bash' end
    return ft_to_lang(ft)
end

