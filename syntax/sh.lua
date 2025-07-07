#!/usr/bin/env lua

vim.api.nvim_create_autocmd("Colorscheme", {
    buffer = 0,
    group = vim.api.nvim_create_augroup("syntax", {}),
    callback = function()
        -- Appended "." to syntax specific iskeyword so that e.g. `git config alias.root ...`
        -- isn't understood to contain the shell builtin "alias"
        vim.cmd.syntax("iskeyword @,48-57,_,192-255,#,-,.")
    end
})
