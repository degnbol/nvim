local utils = require('utils')
local g = vim.g

g.mapleader = ' '
g.maplocalleader = ' '

local opt = {noremap = true, silent = true}

-- Commenting. I haven't moved them to whichkey.lua since it seems that is only for normal mode.
utils.map("n", "<leader>/", ":CommentToggle<CR>", opt)
utils.map("v", "<leader>/", ":CommentToggle<CR>", opt)



