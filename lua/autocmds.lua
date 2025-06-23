#!/usr/bin/env lua
-- default contents in new files of specific types.
-- Either lines, or function producing lines
local defaultlines = {
    -- as a function instead of separate patterns as otherwise they would all get triggered.
    sh = function(filepath)
        if filepath:match("%.activate.sh$") then
            return { "conda activate ENVIRONMENT 2> /dev/null" }
        elseif filepath:match("%.deactivate.sh$") then
            return { "[ \"$CONDA_DEFAULT_ENV\" = base ] || conda deactivate" }
        end
        return { "cd $0:h" }
    end,
    scm = {
        ";extends",
    },
    R = {
        [[if (!require("pacman")) install.packages("pacman")]],
        [[pacman::p_load(data.table, ggplot2, cowplot, ggh4x, svglite)]],
    },
    py = {
        "import numpy as np",
    },
    jl = {
        "using DataFrames, CSV",
    },
    lua = function(filepath)
        -- if writing luasnippets
        if filepath:match("/luasnippets/") then
            return {
                'local ls = require "luasnip"',
                'local s = ls.snippet',
                'local sn = ls.snippet_node',
                'local isn = ls.indent_snippet_node',
                'local t = ls.text_node',
                'local i = ls.insert_node',
                'local f = ls.function_node',
                'local c = ls.choice_node',
                'local d = ls.dynamic_node',
                'local r = ls.restore_node',
                'local events = require("luasnip.util.events")',
                'local ai = require("luasnip.nodes.absolute_indexer")',
                'local extras = require("luasnip.extras")',
                'local l = extras.lambda',
                'local rep = extras.rep',
                'local p = extras.partial',
                'local m = extras.match',
                'local n = extras.nonempty',
                'local dl = extras.dynamic_lambda',
                'local fmt = require("luasnip.extras.fmt").fmt',
                'local fmta = require("luasnip.extras.fmt").fmta',
                'local conds = require("luasnip.extras.expand_conditions")',
                'local postfix = require("luasnip.extras.postfix").postfix',
                'local types = require("luasnip.util.types")',
                'local parse = require("luasnip.util.parser").parse_snippet',
                'local ms = ls.multi_snippet',
                'local k = require("luasnip.nodes.key_indexer").new_key',
                '',
                'local lsu = require "utils/luasnip"',
                'local re = lsu.re', -- my convenience function for regex match
                '',
                'return {',
                '}',
            }
        elseif filepath:match("lua/plugins") then
            return {
                '',
                'return {',
                '}',
            }
        end
        return {}
    end,
    vim = function(filepath)
        if filepath:match("ftdetect/") then
            return { [[au BufNewFile,BufRead *.EXT set filetype=FILETYPE]] }
        end
        return {}
    end,
    -- complgen in config/complgen/*.usage
    usage = function()
        -- assume filename is command name
        local cmdname = vim.fn.expand('%:r')
        return {
            "# || syntax allows <PATH> to take completion precidence over options (otherwise they would be suggested together)",
            cmdname .. " <PATH>... || [<OPTION>]... <PATH>...;",
            "",
            "<OPTION> ::=",
            "      (--help) \"display this help and exit\"",
            "    | (--version) \"display PyMOL version and exit\"",
            "    ;",
            "",
            "# Fix PATH not adding space after completion.",
            [[<PATH@zsh> ::= {{{ _path_files | sed '/\/$/!s/$/ /' }}};]],
        }
    end,
}

local grp = vim.api.nvim_create_augroup("defaultfile", { clear = true })

for ext, lines in pairs(defaultlines) do
    local callback
    if type(lines) == "function" then
        callback = function()
            local filepath = vim.api.nvim_buf_get_name(0)
            vim.api.nvim_put(lines(filepath), "l", false, true)
        end
    else
        callback = function()
            vim.api.nvim_put(lines, "l", false, true)
        end
    end
    vim.api.nvim_create_autocmd("BufNewFile", {
        pattern = "*." .. ext,
        group = grp,
        callback = callback
    })
end


vim.api.nvim_create_autocmd("Filetype", {
    pattern = "*",
    group = vim.api.nvim_create_augroup("cmp-syntax", { clear = true }),
    callback = function()
        if vim.opt_local.omnifunc:get() == "" then
            vim.opt_local.omnifunc = "syntaxcomplete#Complete"
        end
    end
})

-- start treesitter for each new filetype
vim.api.nvim_create_autocmd("Filetype", {
    pattern = "*",
    group = vim.api.nvim_create_augroup("treesitter", { clear = false }),
    callback = function(args)
        local disabled = {
            "vim",               -- not perfect
            "latex",             -- messes with vimtex in lots of ways, e.g. conceal, detection of mathzone, cycling with ts$
            "sh", "bash", "zsh", -- broken
        }
        local additional_vim_regex_highlighting = {
            "vimdoc",   -- treesitter version doesn't contain useful colors from :h group-name
            "bash",     -- spending too much time writing treesitter query. Also covers zsh.
            "markdown", -- my custom comment syntax matches in after/syntax/markdown.vim
            -- Semicolon isn't currently highlighted in all cases by TS so we want to incl vim regex hl for jl.
            -- However, jl can get slowed down a lot in certain files from the syntax hl. The solution:
            -- We enable it, but avoid any default syntax hl and only set custom syntax hl in syntax/julia.vim.
            "julia",
            "sql",  -- custom postgres highlight in syntax/sql.vim
            "wgsl", -- custom in syntax/wgsl.vim
        }
        if disabled[vim.bo.filetype] == nil then
            -- Ignore errors for buffers where treesitter parser is not installed.
            pcall(vim.treesitter.start, args.buf)
            if additional_vim_regex_highlighting[vim.bo.filetype] then
                vim.bo[args.buf].syntax = 'on'
            end
        end
    end
})
