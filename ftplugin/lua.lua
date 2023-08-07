#!/usr/bin/env lua
vim.opt.list = false
-- to wrap comments using gw (default vim formatting)
vim.opt.formatoptions:append("t")

require 'luasnippets/lazy'

-- use :help instead of default Man for visual selection in lua
vim.keymap.set('x', 'K', [["hy:h <C-r>h<CR>]], { buffer=true, desc="Help" })

