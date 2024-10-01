#!/usr/bin/env lua
-- default contents in new files of specific types
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
}

local grp = vim.api.nvim_create_augroup("defaultfile", {clear=true})

for ext, lines in pairs(defaultlines) do
    vim.api.nvim_create_autocmd("BufNewFile", {
        pattern = "*." .. ext,
        group = grp,
        callback = function ()
            vim.api.nvim_put(lines, "l", false, true)
        end
    })
end
