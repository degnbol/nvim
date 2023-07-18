#!/usr/bin/env lua
return {
    -- treesitter
    -- language coloring and ensuring of installation
    {'nvim-treesitter/nvim-treesitter', build=':TSUpdate', config=function()
        require"nvim-treesitter.configs".setup {
            ensure_installed = {
                "bash",
                "c_sharp",
                "lua",
                "json",
                "python",
                "julia",
                "latex",
                -- "java",
                -- "kotlin",
                "help", -- vim help files https://github.com/neovim/tree-sitter-vimdoc
                "r",
                "markdown", -- for block code
                "markdown_inline", -- for inline code
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
                "awk",
                "rust",
                "javascript",
                "scala",
                "sql",
                "graphql", --ext .gql, e.g. schema for graph databases
            },
            highlight = {
                enable = true,
                disable = {
                    "vim", -- not perfect
                    "latex", -- messes with vimtex in lots of ways, e.g. conceal, detection of mathzone, cycling with ts$
                    -- "help", -- removes useful colors from :h group-name
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
                    init_selection = "<c-space>",
                    node_incremental = "<c-space>",
                    scope_incremental = "<c-s>",
                    node_decremental = "<c-backspace>",
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
        vim.treesitter.language.register("bash", "zsh")

        -- TODO: add
        -- https://github.com/nvim-treesitter/nvim-treesitter/#adding-parsers
        -- https://github.com/Beaglefoot/tree-sitter-awk
        -- then make bash/injections.scm that takes command awk raw_string and captures the raw_string with @awk
        -- maybe mlr but would probs have to write it or something


    end},
    -- refactor
    {"nvim-treesitter/nvim-treesitter-refactor", dependencies='nvim-treesitter/nvim-treesitter', config=function()
        require"nvim-treesitter.configs".setup {
            -- for https://github.com/nvim-treesitter/nvim-treesitter-refactor
            refactor = {
                highlight_definitions = { enable = true },
                smart_rename = {
                    enable = true,
                    keymaps = {
                        smart_rename = "<leader>rn",
                    },
                },
                navigation = {
                    enable = true,
                    keymaps = {
                        goto_definition = "gnd",
                        list_definitions = "gnD",
                        list_definitions_toc = "gO",
                        goto_next_usage = "<a-*>",
                        goto_previous_usage = "<a-#>",
                    },
                },
            },
        }
    end},
    -- selecting, moving functions etc.
    {"nvim-treesitter/nvim-treesitter-textobjects", dependencies='nvim-treesitter/nvim-treesitter', config=function()
        require"nvim-treesitter.configs".setup {
            -- for https://github.com/nvim-treesitter/nvim-treesitter-textobjects
            textobjects = {
                select = {
                    enable = true,
                    -- Automatically jump forward to textobj, similar to targets.vim 
                    lookahead = true,

                    keymaps = {
                        -- You can use the capture groups defined in textobjects.scm
                        ["af"] = "@function.outer",
                        ["if"] = "@function.inner",
                        ["ac"] = "@comment.outer",
                        -- ["ic"] = "@comment.inner", -- doesn't exist
                        -- these are so cool:
                        ["aa"] = "@parameter.outer",
                        ["ia"] = "@parameter.inner",
                    }
                },
                swap = {
                    enable = true,
                    swap_next = {
                        ["<leader>a"] = "@parameter.inner",
                    },
                    swap_previous = {
                        ["<leader>A"] = "@parameter.inner",
                    }
                },
                move = {
                    enable = true,
                    set_jumps = true, -- whether to set jumps in the jumplist
                    goto_next_start = {
                        ["]m"] = "@function.outer",
                        ["]o"] = "@loop.outer",
                    },
                    goto_next_end = {
                        ["]M"] = "@function.outer",
                        ["]O"] = "@loop.outer",
                    },
                    goto_previous_start = {
                        ["[m"] = "@function.outer",
                        ["[o"] = "@loop.outer",
                    },
                    goto_previous_end = {
                        ["[M"] = "@function.outer",
                        ["[O"] = "@loop.outer",
                    }
                },
                lsp_interop = {
                    enable = true,
                    border = 'none',
                    peek_definition_code = {
                        ["<leader>df"] = "@function.outer",
                        ["<leader>dF"] = "@class.outer",
                    }
                },
            },
        }
    end},
    -- in vis mode use . , ; i; to select based on treesitter 
    {"RRethy/nvim-treesitter-textsubjects", dependencies='nvim-treesitter/nvim-treesitter', config=function()
        -- while in visual mode these keybindings will change what is selected
        require('nvim-treesitter.configs').setup {
            textsubjects = {
                enable = true,
                -- optional keymap to select the previous selection, which effectively decreases the incremental smart selection
                prev_selection = ',',
                keymaps = {
                    -- this will do incremental expand select
                    ['.'] = 'textsubjects-smart',
                    -- treesitter based container will be selected with va;
                    ['a;'] = 'textsubjects-container-outer',
                    -- inside of treesitter based container will be selected with vi;
                    ['i;'] = 'textsubjects-container-inner',
                },
            },
        }
    end},
    -- show the "context" at the top line, i.e. function name when in a function
    -- "romgrk/nvim-treesitter-context",
    -- error for julia tree-sitter:
    -- % jumps between matching coding blocks, not just single chars.
    -- {"andymass/vim-matchup", dependencies='nvim-treesitter/nvim-treesitter', config=function()
        -- require'nvim-treesitter.configs'.setup {
        --     matchup = {enable = true}
        -- }
    -- end},
    -- tree sitter based rainbow color parenthesis to easily see the matching
    {"mrjones2014/nvim-ts-rainbow", dependencies='nvim-treesitter/nvim-treesitter', config=function()
        require("nvim-treesitter.configs").setup {
            -- for the p00f/nvim-ts-rainbow plugin
            rainbow = {
                enable = true, -- update broke this plugin
                -- disable = {"julia"},
                extended_mode = true, -- Highlight also non-parentheses delimiters, boolean or table: lang -> boolean
                max_file_lines = 1000, -- Do not enable for files with more than 1000 lines, int
            },
        }
    end},
    -- alt that might be better maintained.
    -- I prefer the other for now, since it rainbows "\begin" and "\end" in latex, where this version also colors the following "{...}"
    {"HiPhish/nvim-ts-rainbow2", enabled=false, dependencies='nvim-treesitter/nvim-treesitter', config=function()
        require("nvim-treesitter.configs").setup {
            rainbow = {
                enable = true,
                -- disable = {"julia"},
                -- Which query to use for finding delimiters
                query = {'rainbow-parens', html='rainbow-tags', latex='rainbow-blocks',},
                -- Highlight the entire buffer all at once
                strategy = require('ts-rainbow').strategy.global,
            },
        }
    end},
}

