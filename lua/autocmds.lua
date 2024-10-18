#!/usr/bin/env lua
-- default contents in new files of specific types.
-- Either lines, or function producing lines
local defaultlines = {
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
        "#!/usr/bin/env zsh",
        "conda activate ENVIRONMENT 2> /dev/null"
    },
    ["deactivate.sh"] = {
        "#!/usr/bin/env zsh",
        "[ $CONDA_DEFAULT_ENV = base ] || conda deactivate"
    },
    -- complgen in config/complgen/*.usage
    usage = function ()
        -- assume filename is command name
        local cmdname = vim.fn.expand('%:r')
        lines = {
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
        return lines
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

