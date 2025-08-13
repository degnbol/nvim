
local rtp = vim.opt.runtimepath:get()[1]
dofile(rtp .. "/ftplugin/sh.lua")
vim.opt_local.list = false
vim.opt_local.filetype = "sh.zsh"
