#!/usr/bin/env lua
-- remove o, we want to continue comments while editing them only (r).
-- no t and having c+a means only comments are autoformatted.
vim.opt.formatoptions = "jwcrqla"
