#!/usr/bin/env lua
local rtp = vim.opt.runtimepath:get()[1]
dofile(rtp .. "/ftplugin/sh.lua")
vim.opt_local.list = false
-- bash is effectively sh in neovim. Let's see if we can have zsh separate for
-- ft, and use dofile etc. to load code shared between bash and zsh.
-- This is to avoid zsh being filled with "shIf" etc. i.e. bash alternative
-- highlight groups.
-- vim.opt_local.filetype = "sh.zsh"
