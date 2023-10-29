#!/usr/bin/env lua
-- remove o, we want to continue comments while editing them only (r).
-- no t and having c+a means only comments are autoformatted. 
-- However, a made comment reformat slow, so don't use by default.
vim.opt.formatoptions = "jwcrql"
vim.opt.concealcursor = ""
vim.opt.list = false

local grp = vim.api.nvim_create_augroup("colorscheme", {clear=true})
vim.api.nvim_create_autocmd("Colorscheme", {
    buffer = 0,
    group = grp,
    callback = function ()
        vim.api.nvim_set_hl(0, "@cell", {reverse=true})
    end
})

