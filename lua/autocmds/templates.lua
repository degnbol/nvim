local shebang = "#!/usr/bin/env "

--- Escape literal `$` for snippet syntax.
local function esc(s) return (s:gsub("%$", "\\$")) end

--- Join lines into a snippet body string.
--- Items are either plain strings (escaped) or raw snippet strings
--- wrapped with `raw()`.
local RAW = {} -- sentinel
local function raw(s) return setmetatable({ s }, RAW) end
local function snippet(lines)
    local parts = {}
    for _, line in ipairs(lines) do
        if getmetatable(line) == RAW then
            parts[#parts + 1] = line[1]
        else
            parts[#parts + 1] = esc(line)
        end
    end
    return table.concat(parts, "\n")
end

-- default contents in new files of specific types.
-- Either a snippet body string, or a function returning one.
-- Use raw() for lines containing snippet syntax ($0, $1, $TM_FILENAME, etc.).
local templates = {
    -- as a function instead of separate patterns as otherwise they would all get triggered.
    sh = function()
        -- Get zsh syntax
        vim.api.nvim_set_option_value("filetype", "zsh", { buf = 0 })
        return snippet {
            shebang .. "zsh",
            "set -euo pipefail",
            "cd ${0:A:h}",
            "",
        }
    end,
    zsh = snippet {
        shebang .. "zsh",
        "set -euo pipefail",
        "cd ${0:A:h}",
        "",
    },
    bash = snippet {
        shebang .. "bash",
        "cd $(dirname $0)",
        "",
    },
    scm = snippet {
        ";extends",
        raw "$0",
    },
    R = snippet {
        shebang .. "Rscript",
        [[if (!require("pacman", quiet=TRUE)) install.packages("pacman")]],
        [[pacman::p_load(data.table, ggplot2, cowplot, ggh4x, svglite)]],
        "",
    },
    py = function(filepath)
        if filepath:match("%.pml%.py$") then return nil end
        return snippet {
            shebang .. "python3",
            "import numpy as np",
            "",
        }
    end,
    ["pml.py"] = snippet {
        "#!/bin/sh",
        [[''''exec pymol -cxpqk "$0" -- "$@" # ''']],
        '"""',
        raw "Usage: ./$TM_FILENAME FILE.pdb [ARGS] > OUTFILE.tsv",
        '"""',
        "import sys, os",
        '_script_dir = os.path.dirname(os.path.abspath(sys.argv[0]))',
        "if _script_dir not in sys.path:",
        "    sys.path.insert(0, _script_dir)",
        "from pymol import cmd",
        "import pml",
        "",
        raw "$0",
        "",
        "cmd.quit()",
    },
    jl = snippet {
        shebang .. "julia",
        "using DataFrames, CSV",
        "",
    },
    lua = function(filepath)
        -- if writing luasnippets
        if filepath:match("/luasnippets/") then
            return snippet {
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
                'local re = lsu.re',
                '',
                'return {',
                raw "\t$0",
                '}',
            }
        elseif filepath:match("lua/plugins") then
            return snippet {
                '',
                'return {',
                raw "\t$0",
                '}',
            }
        end
        return nil
    end,
    vim = function(filepath)
        if filepath:match("ftdetect/") then
            return snippet {
                raw "au BufNewFile,BufRead *.${1:$TM_FILENAME_BASE} set filetype=$TM_FILENAME_BASE",
                "",
            }
        end
        return nil
    end,
    -- complgen in config/complgen/*.usage
    usage = function()
        -- assume filename is command name
        local cmdname = vim.fn.expand('%:r')
        return snippet {
            "# || syntax allows <PATH> to take completion precidence over options (otherwise they would be suggested together)",
            cmdname .. " <PATH>... || [<OPTION>]... <PATH>...;",
            "",
            "<OPTION> ::=",
            '      (--help) "display this help and exit"',
            '    | (--version) "display PyMOL version and exit"',
            "    ;",
            "",
            "# Fix PATH not adding space after completion.",
            [[<PATH@zsh> ::= {{{ _path_files | sed '/\/$/!s/$/ /' }}};]],
            "",
        }
    end,
}

-- Disable neovim's built-in Tab/S-Tab snippet jump mappings.
-- We use <C-.>/<C-,> instead (consistent with blink.cmp).
vim.keymap.del({ "i", "s" }, "<Tab>")
vim.keymap.del({ "i", "s" }, "<S-Tab>")

-- Has to be outside here defined as variable group when clear=true.
-- Otherwise will reset for each iteration of the loop.
local grp = vim.api.nvim_create_augroup("defaultfile", { clear = true })
for ext, body in pairs(templates) do
    vim.api.nvim_create_autocmd("BufNewFile", {
        pattern = "*." .. ext,
        group = grp,
        callback = function()
            local b = body
            if type(b) == "function" then
                b = b(vim.api.nvim_buf_get_name(0))
            end
            if b then
                vim.snippet.expand(b)
                -- Stay in select mode for placeholders, otherwise return to normal mode.
                if vim.fn.mode() ~= 's' then
                    vim.cmd.stopinsert()
                end
            end
        end,
    })
end
