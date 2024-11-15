#!/usr/bin/env lua
-- default contents in new files of specific types.
-- Either lines, or function producing lines
local defaultlines = {
    sh = {
        "cd $0:h",
    },
    scm = {
        ";extends",
    },
    R = {
        "suppressPackageStartupMessages(library(data.table))",
        "suppressPackageStartupMessages(library(ggplot2))",
    },
    py = {
        "import numpy as np",
    },
    jl = {
        "using DataFrames, CSV",
    },
    ["activate.sh"] = {
        "conda activate ENVIRONMENT 2> /dev/null"
    },
    ["deactivate.sh"] = {
        "[ \"$CONDA_DEFAULT_ENV\" = base ] || conda deactivate"
    },
    lua = function ()
        local filepath = vim.api.nvim_buf_get_name(0)
        -- if writing luasnippets
        if filepath:match("/luasnippets/") then
            return {
                'local ls = require "luasnip"',
                'local extras = require "luasnip/extras"',
                'local s = ls.s',
                'local t = ls.t',
                'local f = ls.f',
                'local c = ls.c',
                'local fmta = extras.fmta',
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
    -- complgen in config/complgen/*.usage
    usage = function ()
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

local grp = vim.api.nvim_create_augroup("defaultfile", {clear=true})

for ext, lines in pairs(defaultlines) do
    local callback
    if type(lines) == "function" then
        callback = function ()
            vim.api.nvim_put(lines(), "l", false, true)
        end
    else
        callback = function ()
            vim.api.nvim_put(lines, "l", false, true)
        end
    end
    vim.api.nvim_create_autocmd("BufNewFile", {
        pattern = "*." .. ext,
        group = grp,
        callback = callback
    })
end


local grp = vim.api.nvim_create_augroup("cmp-syntax", {clear=true})
vim.api.nvim_create_autocmd("Filetype", {
    pattern = "*",
    group = grp,
    callback = function ()

        if vim.opt_local.omnifunc:get() == "" then
            vim.opt_local.omnifunc = "syntaxcomplete#Complete"
        end
    end
})

