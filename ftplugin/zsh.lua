#!/usr/bin/env lua
local rtp = vim.opt.runtimepath:get()[1]
dofile(rtp .. "/ftplugin/sh.lua")
vim.opt.list = false
vim.opt_local.filetype = "sh.zsh"

