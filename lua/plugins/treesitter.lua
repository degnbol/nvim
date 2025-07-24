return {
    -- treesitter
    -- language coloring and ensuring of installation
    {
        'nvim-treesitter/nvim-treesitter',
        enabled = true,
        branch = "master", -- master is frozen for backwards compatability, main will become default in future.
        build = ':TSUpdate',
        config = function()
            -- TODO: replicate these master branch settings with what is equivalent for the main branch version.
            -- require "nvim-treesitter.configs".setup {
            --     ensure_installed = {
            --         "awk",
            --         "bash",
            --         "c_sharp",
            --         "lua",
            --         "json",
            --         "python",
            --         "julia",
            --         "matlab",
            --         "latex",
            --         -- "java",
            --         -- "kotlin",
            --         "vimdoc",
            --         "r",
            --         "markdown",        -- for block code
            --         "markdown_inline", -- for inline code
            --         "toml",
            --         "vim",
            --         "regex",
            --         "make",
            --         -- "norg",
            --         "cmake",
            --         "cpp",
            --         "bibtex",
            --         "gitcommit",
            --         "gitignore",
            --         "gitattributes",
            --         "diff",  -- for diff output https://github.com/the-mikedavis/tree-sitter-diff
            --         "query", -- what treesitter queries (*.scm) are written in
            --         "awk",
            --         "rust",
            --         "javascript",
            --         "scala",
            --         "sql",
            --         "graphql", --ext .gql, e.g. schema for graph databases
            --         "qf",      -- see below. Run :TSInstall qf
            --     },
            --     incremental_selection = {
            --         enable = true,
            --         keymaps = {
            --             init_selection = "<leader><up>",
            --             node_incremental = "<leader><up>",
            --             node_decremental = "<leader><down>",
            --             scope_incremental = "<leader><left>",
            --         },
            --     },
            --     -- opt-in to using treesitter for https://github.com/andymass/vim-matchup
            --     matchup = {
            --         enable = true,
            --         disable = {}, -- optional, list of language that will be disabled
            --     }
            -- }

            -- use bash treesitter for zsh since zsh is very basic
            -- it doesn't work well enough so just disable TS for zsh.
            -- vim.treesitter.language.register("bash", "zsh")

            -- TODO: add
            -- https://github.com/nvim-treesitter/nvim-treesitter/#adding-parsers
            -- https://github.com/Beaglefoot/tree-sitter-awk
            -- then make bash/injections.scm that takes command awk raw_string and captures the raw_string with @awk
            -- maybe mlr but would probs have to write it or something

            -- Add custom parsers.
            -- Typst parser and queries not up-to-date on https://github.com/nvim-treesitter/nvim-treesitter
            vim.api.nvim_create_autocmd('User', {
                pattern = 'TSUpdate',
                group = vim.api.nvim_create_augroup("treesitter", { clear = false }),
                callback = function()
                    local parsers = require('nvim-treesitter.parsers')
                    parsers.typst = {
                        install_info = {
                            url = "https://github.com/uben0/tree-sitter-typst",
                            files = { "src/parser.c", "src/scanner.c" },
                        }
                    }
                    parsers.qf = {
                        install_info = {
                            url = "https://github.com/OXY2DEV/tree-sitter-qf",
                            files = { "src/parser.c" },
                            branch = "main",
                        },
                    }
                end
            })
        end
    },
    -- The queries aren't up to date from TS above so we also need the repo here so we overwrite queries/typst/.
    "uben0/tree-sitter-typst",
    -- refactor
    {
        "nvim-treesitter/nvim-treesitter-refactor",
        enabled = false,
        dependencies = 'nvim-treesitter/nvim-treesitter',
        config = function()
            require "nvim-treesitter.configs".setup {
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
        enabled = true,
        dependencies = 'nvim-treesitter/nvim-treesitter',
        config = function()
            require 'nvim-treesitter.configs'.setup {
                select = {
                    enable = true,
                    -- Automatically jump forward to textobj, similar to targets.vim
                    lookahead = true,

                    keymaps = {
                        -- You can use the capture groups defined in textobjects.scm
                        ["af"] = "@function.outer",
                        ["if"] = "@function.inner",
                        ["ac"] = "@comment.outer",
                        -- doesn't seem to be supported much. We use a kana derived textobj plugin for ic.
                        ["ic"] = "@comment.inner",
                        -- iC and aC are used for multiline comment elsewhere
                        ["a?"] = "@conditional.outer",
                        ["i?"] = "@conditional.inner",
                        ["ao"] = "@loop.outer",
                        ["io"] = "@loop.inner",
                        ["aa"] = "@parameter.outer",
                        ["ia"] = "@parameter.inner",
                        ["ik"] = "@assignment.lhs",
                        ["iv"] = "@assignment.rhs",
                        ["ab"] = "@block.outer",
                        ["ib"] = "@block.inner",
                        -- TODO: textobj for whatever the current enclosing container is, similar to nmap <leader><Up>
                    }
                },
                swap = {
                    enable = true,
                    swap_next = {
                        ["<leader>a]"]       = { query = "@parameter.inner", desc = "Swap next arg" },
                        ["<leader>a<Right>"] = { query = "@parameter.inner", desc = "Swap next arg" },
                        ["<leader>al"]       = { query = "@parameter.inner", desc = "Swap next arg" },
                    },
                    swap_previous = {
                        ["<leader>a["]      = { query = "@parameter.inner", desc = "Swap prev arg" },
                        ["<leader>a<Left>"] = { query = "@parameter.inner", desc = "Swap prev arg" },
                        ["<leader>ah"]      = { query = "@parameter.inner", desc = "Swap prev arg" },
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
                        -- similar to hover help so we use similar keymap as hover.
                        -- Hover currently uses gh for "go hover", similar to gd, gf, etc.
                        -- gh and gH are used for starting select mode by default which we never use.
                        ["gH"] = { query = "@function.outer", desc = "Peek function" },
                        ["g<C-H>"] = { query = "@class.outer", desc = "Peek class" },
                    }
                },
            }

            -- TODO: if we are not yet inside an arg, it should behave like ]<C-a>
            -- It would also be nice to write this in a cleaner way using lua
            -- by taking inspo from the functions being called here maybe.
            -- Low priority since it does the main thing we want.
            vim.keymap.set('n', ']a', "via]<C-a><Esc>", { remap = true, desc = "Next arg start" })
            vim.keymap.set('n', ']A', "via]<C-a>ia<Esc>", { remap = true, desc = "Next arg end" })
            vim.keymap.set('n', '[a', "viao[<C-S-a>iao<Esc>", { remap = true, desc = "Prev arg start" })
            vim.keymap.set('n', '[A', "viao[<C-S-a><Esc>", { remap = true, desc = "Prev arg end" })
        end
    },
    -- in vis mode use . , ; i; to select based on treesitter
    {
        "RRethy/nvim-treesitter-textsubjects",
        dependencies = 'nvim-treesitter/nvim-treesitter',
        enabled = false,
        config = function()
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
    {
        "romgrk/nvim-treesitter-context",
        opts = {
            max_lines = 1,
            min_window_height = 15, -- hide on small windows
        }
    },
    -- error for julia tree-sitter:
    -- % jumps between matching coding blocks, not just single chars.
    -- {
    -- "andymass/vim-matchup", dependencies='nvim-treesitter/nvim-treesitter', config=function()
    -- require'nvim-treesitter.configs'.setup {
    --     matchup = {enable = true}
    -- }
    -- end},
    {
        "https://gitlab.com/HiPhish/rainbow-delimiters.nvim",
        name = "rainbow_delimiters",
        enabled = false,
        opts = {
            strategy = {
                [''] = 'rainbow-delimiters.strategy.global',
            },
            query = {
                [''] = 'rainbow-delimiters',
                latex = 'rainbow-blocks', -- doesn't work otherwise
            },
            highlight = {
                'RainbowDelimiterViolet',
                'RainbowDelimiterOrange',
                'RainbowDelimiterYellow',
                'RainbowDelimiterBlue',
                'RainbowDelimiterGreen',
                'RainbowDelimiterCyan',
                'RainbowDelimiterRed',
            },
        }
    }
}
