local hi = require "utils/highlights"
local map = require "utils/keymap"

return {
    -- Convenience of installing parsers for given languages.
    {
        'nvim-treesitter',
        after = function()
            local nvim_treesitter = require 'nvim-treesitter'
            nvim_treesitter.setup {}

            -- Register custom parsers (not in the default registry) BEFORE
            -- nvim_treesitter.install() so install() can resolve their names.
            -- nvim-treesitter fires `User TSUpdate` from reload_parsers() at
            -- the start of install/update.
            local wesl_path = vim.pack.get({ 'tree-sitter-wesl' }, { info = false })[1].path
            vim.api.nvim_create_autocmd('User', {
                pattern = 'TSUpdate',
                callback = function()
                    local parsers = require('nvim-treesitter.parsers')
                    parsers.openscad = {
                        install_info = {
                            url = 'https://github.com/bollian/tree-sitter-openscad',
                            branch = 'master',
                            queries = 'queries',
                        },
                    }
                    -- vim.pack checks out the latest release (pack_specs.lua).
                    -- queries/ uses Helix capture names; queries/nvim/ uses nvim's.
                    parsers.wesl = {
                        install_info = {
                            path = wesl_path,
                            queries = 'queries/nvim',
                        },
                    }
                    -- The local grammars' queries resolve off the rtp: modules/ is in
                    -- runtimepath, in standard queries/<lang>/ layout.
                    parsers.pymol_select = {
                        install_info = {
                            path = vim.fn.stdpath('config') .. '/modules/tree-sitter-pymol-select',
                        },
                    }
                    parsers.miller = {
                        install_info = {
                            path = vim.fn.stdpath('config') .. '/modules/tree-sitter-miller',
                        },
                    }
                    parsers.smarts = {
                        install_info = {
                            path = vim.fn.stdpath('config') .. '/modules/tree-sitter-smarts',
                        },
                    }
                end
            })

            -- Parsers require the tree-sitter CLI (>= 0.25.0) for generate + build.
            if vim.fn.executable("tree-sitter") == 1 then
                local task = nvim_treesitter.install {
                    "jq",
                    "awk",
                    "zsh", -- also used for bash
                    "c_sharp",
                    "lua",
                    "json",
                    "python",
                    "julia",
                    "perl", -- one-liners injected from zsh (`perl -pe '…'`)
                    "matlab",
                    "latex",
                    "typst",
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
                    "just",
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
                    "tsv", -- also .tab and .bed; the csv parser requires it too
                    "graphql", --ext .gql, e.g. schema for graph databases
                    -- "qf",      -- see below. Run :TSInstall qf
                    "yaml",
                    "ini",
                    "css",
                    "wesl",         -- registered above
                    "pymol_select", -- local grammar registered above
                    "miller",       -- local grammar registered above
                    "smarts",       -- local grammar registered above
                }
                -- Expose task so headless scripts can wait on it:
                -- nvim --headless +"lua require('nvim-treesitter')._install_task:wait()" +qa
                nvim_treesitter._install_task = task
            end
        end
    },
    -- Moving to and swapping functions etc. Loads at startup: mini.ai
    -- (lua/textobjects.lua) reads its `textobjects` queries.
    {
        "nvim-treesitter-textobjects",
        after = function()
            require("nvim-treesitter-textobjects").setup {
                move = {
                    -- whether to set jumps in the jumplist
                    set_jumps = true,
                },
            }

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
        end,
        before = function()
            local function map_swaps(lhs_next, lhs_previous, obj)
                obj = '@' .. obj
                map.n(lhs_next, function()
                    require "nvim-treesitter-textobjects.swap".swap_next(obj)
                end, "Swap " .. obj .. " next")
                map.n(lhs_previous, function()
                    require "nvim-treesitter-textobjects.swap".swap_previous(obj)
                end, "Swap " .. obj .. " previous")
            end
            local function map_move(obj, inner_or_outer, suffix)
                local _obj = "@" .. obj .. "." .. inner_or_outer
                map.nox("]" .. suffix, function()
                    require("nvim-treesitter-textobjects.move").goto_next_start(_obj, "textobjects")
                end)
                map.nox("[" .. suffix, function()
                    require("nvim-treesitter-textobjects.move").goto_previous_start(_obj, "textobjects")
                end)
                local end_suffix = suffix:upper()
                if end_suffix ~= suffix then
                    map.nox("]" .. end_suffix, function()
                        require("nvim-treesitter-textobjects.move").goto_next_end(_obj, "textobjects")
                    end)
                    map.nox("[" .. end_suffix, function()
                        require("nvim-treesitter-textobjects.move").goto_previous_end(_obj, "textobjects")
                    end)
                end
            end

            map_swaps("]<C-a>", "[<C-a>", "parameter.inner")

            map_move("parameter", "inner", "a")
            map_move("function", "outer", "f")
            map_move("loop", "outer", "o")
        end,
    },
    -- show the "context" at the top line, i.e. function name when in a function
    {
        "nvim-treesitter-context",
        --- @type TSContext.UserConfig
        after = function()
            require("treesitter-context").setup {
                max_lines = 1,
                min_window_height = 15, -- Hide on small windows.
                multiwindow = true,     -- Show context in inactive windows.
                on_attach = function(buf)
                    -- Keep context off the prompt box (AgenticInput), but allow
                    -- it on the chat buffer (AgenticChat): injected code blocks
                    -- pin their own function/class context, and a scoped
                    -- queries/agentic/context.scm pins the tool-call ### header.
                    local ft = vim.bo[buf].filetype
                    return not ft:find("^AgenticInput")
                end,
            }
        end,
        before = function()
            map.n("g<up>", function()
                require("treesitter-context").go_to_context(vim.v.count1)
            end, "Goto TS context", { silent = true })

            hi.onColorScheme(function()
                hi.clear("TreesitterContext")
                local border = { underdashed = true, special = "gray" }
                hi.set("TreesitterContextBottom", border)
                hi.set("TreesitterContextLineNumberBottom", border)
            end)
        end,
    },
}
