#!/usr/bin/env lua

local g = vim.g

g.repl_filetype_commands = {
    python = 'bpython -q'
}


-- keymaps
-- TODO make a function where if the python line starts with def then send V]M:ReplSend<CR>`>. 
-- In Julia it should instead look for function anywhere in the line and call V][:ReplSend`>
vim.api.nvim_set_keymap('n', '<localleader><CR>', ':ReplSend<CR><CR>', {expr = false, noremap = false})
-- After sending to visual the cursor jumps to the start of the selection. 
-- `> means go to mark named > which will be at the end of the previous selection.
vim.api.nvim_set_keymap('v', '<localleader><CR>', ':ReplSend<CR>`>', {expr = false, noremap = false})
vim.api.nvim_set_keymap('n', '<localleader>r', ':ReplToggle<CR>', {expr = false, noremap = false})
vim.api.nvim_set_keymap('n', '<localleader>c', ':ReplClear<CR>', {expr = false, noremap = false})
vim.api.nvim_set_keymap('n', '<localleader><tab>', '<c-w>wA', {expr = false, noremap = false})
vim.api.nvim_set_keymap('v', '<localleader><tab>', 'y<c-w>wpA', {expr = false, noremap = false})


