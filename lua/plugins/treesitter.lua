#!/usr/bin/env lua

return {
    -- treesitter
    -- language coloring and ensuring of installation
    {
        'nvim-treesitter/nvim-treesitter',
        build=':TSUpdate',
        config=function()
            require"nvim-treesitter.configs".setup {
                ensure_installed = {
                    "bash",
                    "c_sharp",
                    "lua",
                    "json",
                    "python",
                    "julia",
                    "matlab",
                    "latex",
                    -- "java",
                    -- "kotlin",
                    "help", -- vim help files https://github.com/neovim/tree-sitter-vimdoc
                    "r",
                    -- some error
                    -- "markdown", -- for block code
                    -- "markdown_inline", -- for inline code
                    "toml",
                    "vim",
                    "regex",
                    "make",
                    -- "norg",
                    "cmake",
                    "cpp",
                    "bibtex",
                    "gitignore",
                    "gitattributes",
                    "diff", -- for diff output https://github.com/the-mikedavis/tree-sitter-diff
                    "query", -- what treesitter queries (*.scm) are written in
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
            -- it doesn't work well enough so just disable TS for zsh.
            -- vim.treesitter.language.register("bash", "zsh")

            -- TODO: add
            -- https://github.com/nvim-treesitter/nvim-treesitter/#adding-parsers
            -- https://github.com/Beaglefoot/tree-sitter-awk
            -- then make bash/injections.scm that takes command awk raw_string and captures the raw_string with @awk
            -- maybe mlr but would probs have to write it or something


        end
    },
    -- refactor
    {
        "nvim-treesitter/nvim-treesitter-refactor",
        dependencies='nvim-treesitter/nvim-treesitter',
        config=function()
            require"nvim-treesitter.configs".setup {
                -- for https://github.com/nvim-treesitter/nvim-treesitter-refactor
                refactor = {
                    highlight_definitions = { enable = true },
                    smart_rename = {
                        enable = true,
                        keymaps = {
                            -- same as LSP which replaces the mapping if attached.
                            smart_rename = "<leader>rn",
                        },
                    },
                    navigation = {
                        enable = true,
                        keymaps = {
                            -- gets replaced if LSP is attached to buffer
                            goto_definition = "gd",
                            list_definitions = "gD",
                            list_definitions_toc = "gO",
                            goto_next_usage = "<a-*>",
                            goto_previous_usage = "<a-#>",
                        },
                    },
                },
            }
        end
    },
    -- selecting, moving functions etc.
    {
        "nvim-treesitter/nvim-treesitter-textobjects",
        dependencies='nvim-treesitter/nvim-treesitter',
        config=function()
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
                            ["ic"] = "@comment.inner",
                            ["aa"] = "@parameter.outer",
                            ["ia"] = "@parameter.inner",
                            ["ik"] = "@assignment.lhs",
                            ["ik"] = "@assignment.lhs",
                            ["iv"] = "@assignment.rhs",
                        }
                    },
                    swap = {
                        enable = true,
                        swap_next = {
                            ["<leader>a]"]       = {query="@parameter.inner", desc="Swap next arg"},
                            ["<leader>a<Right>"] = {query="@parameter.inner", desc="Swap next arg"},
                            ["<leader>al"]       = {query="@parameter.inner", desc="Swap next arg"},
                        },
                        swap_previous = {
                            ["<leader>a["]       = {query="@parameter.inner", desc="Swap prev arg"},
                            ["<leader>a<Left>"]  = {query="@parameter.inner", desc="Swap prev arg"},
                            ["<leader>ah"]       = {query="@parameter.inner", desc="Swap prev arg"},
                        }
                    },
                    move = {
                        enable = true,
                        set_jumps = true, -- whether to set jumps in the jumplist
                        goto_next_start = {
                            -- We want to go to next arg at the same level, 
                            -- e.g. skip over args that are in a nested 
                            -- function. This means we use <C-a> for the 
                            -- textobj and define a more messy keybind further 
                            -- down.
                            ["]<C-a>"] = "@parameter.inner",
                            ["]f"] = "@function.outer",
                            ["]o"] = "@loop.outer",
                            ["]z"] = { query = "@fold", query_group = "folds", desc = "Next fold" },
                        },
                        goto_next_end = {
                            ["]<C-S-a>"] = "@parameter.inner",
                            ["]F"] = "@function.outer",
                            ["]O"] = "@loop.outer",
                        },
                        goto_previous_start = {
                            ["[<C-a>"] = "@parameter.inner",
                            ["[f"] = "@function.outer",
                            ["[o"] = "@loop.outer",
                            ["[z"] = { query = "@fold", query_group = "folds", desc = "Prev fold" },
                        },
                        goto_previous_end = {
                            ["[<C-S-a>"] = "@parameter.inner",
                            ["[F"] = "@function.outer",
                            ["[O"] = "@loop.outer",
                        }
                    },
                    lsp_interop = {
                        enable = true,
                        border = 'none',
                        peek_definition_code = {
                            ["<leader>df"] = {query="@function.outer", desc="Peek function"},
                            ["<leader>dF"] = {query="@class.outer", desc="Peek class"},
                        }
                    },
                },
            }

            -- TODO: if we are not yet inside an arg, it should behave like ]<C-a>
            -- It would also be nice to write this in a cleaner way using lua 
            -- by taking inspo from the functions being called here maybe.
            -- Low priority since it does the main thing we want.
            vim.keymap.set('n', ']a', "via]<C-a><Esc>", { remap=true, desc="Next arg start" })
            vim.keymap.set('n', ']A', "via]<C-a>ia<Esc>", { remap=true, desc="Next arg end" })
            vim.keymap.set('n', '[a', "viao[<C-S-a>iao<Esc>", { remap=true, desc="Prev arg start" })
            vim.keymap.set('n', '[A', "viao[<C-S-a><Esc>", { remap=true, desc="Prev arg end" })
        end
    },
    -- in vis mode use . , ; i; to select based on treesitter 
    {
        "RRethy/nvim-treesitter-textsubjects",
        dependencies='nvim-treesitter/nvim-treesitter',
        config=function()
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
        end
    },
    -- show the "context" at the top line, i.e. function name when in a function
    -- "romgrk/nvim-treesitter-context",
    -- error for julia tree-sitter:
    -- % jumps between matching coding blocks, not just single chars.
    -- {
    -- "andymass/vim-matchup", dependencies='nvim-treesitter/nvim-treesitter', config=function()
    -- require'nvim-treesitter.configs'.setup {
    --     matchup = {enable = true}
    -- }
    -- end},
    -- tree sitter based rainbow color parenthesis to easily see the matching
    {
        "mrjones2014/nvim-ts-rainbow",
        dependencies='nvim-treesitter/nvim-treesitter',
        config=function()
            require("nvim-treesitter.configs").setup {
                -- for the p00f/nvim-ts-rainbow plugin
                rainbow = {
                    enable = true, -- update broke this plugin
                    -- disable = {"julia"},
                    extended_mode = true, -- Highlight also non-parentheses delimiters, boolean or table: lang -> boolean
                    max_file_lines = 1000, -- Do not enable for files with more than 1000 lines, int
                },
            }
        end
    },
    -- alt that might be better maintained.
    -- I prefer the other for now, since it rainbows "\begin" and "\end" in latex, where this version also colors the following "{...}"
    {
        "HiPhish/nvim-ts-rainbow2",
        enabled=false,
        dependencies='nvim-treesitter/nvim-treesitter',
        config=function()
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
        end
    },
}

