#!/usr/bin/env lua
vim.opt.list = false
-- to wrap comments using gw (default vim formatting)
vim.opt.formatoptions:append("t")

require 'luasnippets/lazy'

-- use :help instead of default Man for visual selection in lua
vim.keymap.set('x', 'K', [["hy:h <C-r>h<CR>]], { buffer=true, desc="Help" })

-- when doing gf or similar obviously we should look in the lua/ folder since 
-- this is where all scripts are required from.
local rtp = vim.opt.runtimepath:get()[1]
vim.opt_local.path:append(rtp .. "/lua")

vim.cmd.iabbrev("ture", "true")
vim.cmd.iabbrev("flase", "false")

