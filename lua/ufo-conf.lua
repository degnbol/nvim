#!/usr/bin/env lua
-- https://github.com/kevinhwang91/nvim-ufo

-- 1 char width column showing icons to indicate folding levels
-- vim.o.foldcolumn = '1'
-- pretty icons
vim.o.fillchars = [[eob: ,fold: ,foldopen:,foldsep: ,foldclose:]]

vim.keymap.set('n', 'zR', require('ufo').openAllFolds)
vim.keymap.set('n', 'zM', require('ufo').closeAllFolds)

require('ufo').setup()


