local hi = require "utils/highlights"
local map = require "utils/keymap"

return {
    -- Convenience of installing parsers for given languages.
    {
        'nvim-treesitter/nvim-treesitter',
        branch = "main", -- master is frozen for backwards compatibility, main will become default in future.
        build = ':TSUpdate',
        config = function()
            local nvim_treesitter = require 'nvim-treesitter'
            nvim_treesitter.setup {}

            -- TS installations will keep retrying unless the TS CLI is installed (with cargo).
            if vim.fn.executable("cargo") == 1 then
                nvim_treesitter.install {
                    -- "bash", -- so broken
                    "c_sharp",
                    "lua",
                    "json",
                    "python",
                    "julia",
                    "matlab",
                    "latex",
                    -- "java",
                    -- "kotlin",
                    "vimdoc",
                    "r",
                    "markdown",        -- for block code
                    "markdown_inline", -- for inline code
                    "toml",
                    "vim",
                    "regex",
                    "make",
                    -- "norg",
                    "cmake",
                    "cpp",
                    "bibtex",
                    "gitcommit",
                    "gitignore",
                    "git_config",
                    "gitattributes",
                    "diff",  -- for diff output https://github.com/the-mikedavis/tree-sitter-diff
                    "query", -- what treesitter queries (*.scm) are written in
                    "awk",
                    "rust",
                    "javascript",
                    "scala",
                    "sql",
                    "graphql", --ext .gql, e.g. schema for graph databases
                    -- "qf",      -- see below. Run :TSInstall qf
                    "yaml",
                    "ini",
                    "css",
                }
            end


            -- Add custom openscad parser (not in default list)
            vim.api.nvim_create_autocmd('User', {
                pattern = 'TSUpdate',
                callback = function()
                    require('nvim-treesitter.parsers').openscad = {
                        install_info = {
                            url = 'https://github.com/bollian/tree-sitter-openscad',
                            branch = 'master',
                            queries = 'queries', -- directory with highlights.scm etc.
                        },
                    }
                end
            })

            -- TODO: add
            -- https://github.com/nvim-treesitter/nvim-treesitter/#adding-parsers
            -- https://github.com/Beaglefoot/tree-sitter-awk
            -- then make bash/injections.scm that takes command awk raw_string and captures the raw_string with @awk
            -- maybe mlr but would probs have to write it or something
        end
    },
    -- Selecting, moving functions etc.
    {
        "nvim-treesitter/nvim-treesitter-textobjects",
        branch = "main",
        opts = {
            select = {
                -- Jump forward to textobj if not already inside it.
                lookahead = true,
                -- Other opts.
                -- selection_modes allows for using a different selection
                -- mode than v, e.g. V or <C-v>, for specific textobjects.
            },
            move = {
                -- whether to set jumps in the jumplist
                set_jumps = true,
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
        },
        init = function()
            local suffixes = {
                ['function'] = 'f',
                -- Also an option to use comma as in the non-TS version in
                -- chrisgrieser/nvim-various-textobjs
                -- The benefit is having a available for the args passed to neovim which is default.
                -- The downside is then goto end of arg would be '<', which we
                -- are currently using for goto less indented.
                parameter = 'a',
                loop = 'o',
                comment = 'c',
                conditional = '?',
                ["assignment.lhs"] = 'k',
                ["assignment.rhs"] = 'v',
                block = 'B',
            }

            local function map_select(lhs, obj)
                map.ox(lhs, function()
                    require "nvim-treesitter-textobjects.select".select_textobject(obj, "textobjects")
                end)
            end
            local function map_select_i(obj)
                map_select("i" .. suffixes[obj], "@" .. obj .. ".inner")
            end
            local function map_select_a(obj)
                map_select("a" .. suffixes[obj], "@" .. obj .. ".outer")
            end
            local function map_selects(obj)
                map_select_i(obj)
                map_select_a(obj)
            end
            local function map_swaps(lhs_next, lhs_previous, obj)
                obj = '@' .. obj
                map.n(lhs_next, function()
                    require "nvim-treesitter-textobjects.swap".swap_next(obj)
                end, "Swap " .. obj .. " next")
                map.n(lhs_previous, function()
                    require "nvim-treesitter-textobjects.swap".swap_previous(obj)
                end, "Swap " .. obj .. " previous")
            end
            local function map_move(obj, inner_or_outer)
                local _obj = "@" .. obj .. "." .. inner_or_outer
                local suffix = suffixes[obj]
                map.nox("]" .. suffix, function()
                    require("nvim-treesitter-textobjects.move").goto_next_start(_obj, "textobjects")
                end)
                map.nox("[" .. suffix, function()
                    require("nvim-treesitter-textobjects.move").goto_previous_start(_obj, "textobjects")
                end)
                suffix = suffix:upper()
                if suffix ~= suffixes[obj] then
                    map.nox("]" .. suffix, function()
                        require("nvim-treesitter-textobjects.move").goto_next_end(_obj, "textobjects")
                    end)
                    map.nox("[" .. suffix, function()
                        require("nvim-treesitter-textobjects.move").goto_previous_end(_obj, "textobjects")
                    end)
                end
            end

            map_selects("function")
            -- iC and aC are used for multiline comment elsewhere
            -- map_obj('C', "class")
            map_select_a("comment") -- comment.inner doesn't seem well supported.
            map_selects("conditional")
            map_selects("loop")
            map_selects("parameter")
            map_selects("assignment.lhs")
            map_selects("assignment.rhs")
            -- b for bracket.
            -- B already works for block textobjs builtin. see :h iB.
            -- map_selects("block")

            map_swaps("]<C-a>", "[<C-a>", "parameter.inner")

            map_move("parameter", "inner")
            map_move("function", "outer")
            map_move("loop", "outer")

            -- Repeat movement with ; and ,
            local ts_repeat_move = require "nvim-treesitter-textobjects.repeatable_move"
            -- vim way: ; goes to the direction you were moving.
            map.nox(";", ts_repeat_move.repeat_last_move)
            map.nox(",", ts_repeat_move.repeat_last_move_opposite)
            -- Also, make builtin f, F, t, T also repeatable with ; and ,
            map.nox("f", ts_repeat_move.builtin_f_expr, nil, { expr = true })
            map.nox("F", ts_repeat_move.builtin_F_expr, nil, { expr = true })
            map.nox("t", ts_repeat_move.builtin_t_expr, nil, { expr = true })
            map.nox("T", ts_repeat_move.builtin_T_expr, nil, { expr = true })
        end
    },
    -- show the "context" at the top line, i.e. function name when in a function
    {
        "romgrk/nvim-treesitter-context",
        opts = {
            max_lines = 1,
            min_window_height = 15, -- Hide on small windows.
            multiwindow = true,     -- Show context in inactive windows.
        },
        init = function()
            map.n("g<up>", function()
                require("treesitter-context").go_to_context(vim.v.count1)
            end, "Goto TS context", { silent = true })

            hi.afterColorscheme(function()
                hi.clear("TreesitterContext")
                local border = { underdashed = true, special = "gray" }
                hi.set("TreesitterContextBottom", border)
                hi.set("TreesitterContextLineNumberBottom", border)
            end)
        end,
    },
}
