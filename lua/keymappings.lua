local utils = require('utils')
local g = vim.g

g.mapleader = ' '
g.maplocalleader = ' '

local opt = {noremap = true, silent = true}

-- Commenting
utils.map("n", "<leader>/", ":CommentToggle<CR>", opt)
utils.map("v", "<leader>/", ":CommentToggle<CR>", opt)

