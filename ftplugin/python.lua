#!/usr/bin/env lua
-- remove o, we want to continue comments while editing them only (r).
-- no t and having c+a means only comments are autoformatted. 
-- However, a made comment reformat slow, so don't use by default.
vim.opt.formatoptions = "jwcrql"
vim.opt.concealcursor = ""
