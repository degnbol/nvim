#!/usr/bin/env lua
-- autocmd to extend a query file by default.
local grp = vim.api.nvim_create_augroup("Extends", {clear=true})
vim.api.nvim_create_autocmd("BufNewFile", {
    pattern = "*.scm",
    group = grp,
    callback = function ()
        vim.api.nvim_put({";extends"}, "l", false, true)
    end
})

