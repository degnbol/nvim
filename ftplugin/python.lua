#!/usr/bin/env lua
local util = require "utils/init"

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

function Load_pymol_snippets()
    if not Loaded_pymol then
        Loaded_pymol = true
        require("luasnip").add_snippets("python",
            require 'luasnippets.python_pymol')
    end
end
-- manually load
vim.keymap.set('n', '<leader><leader>s', Load_pymol_snippets, { desc="Load pymol snippets." })
-- check if pymol is loaded by scanning first 10 lines
for _, line in ipairs(vim.api.nvim_buf_get_lines(0, 0, 10, false)) do
    -- might be using e.g. `from pymol_util import *`
    if line:match("import") and line:match("pymol") then
        return Load_pymol_snippets()
    end
end

