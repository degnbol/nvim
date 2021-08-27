#!/usr/bin/env lua

local g = vim.g

g.repl_filetype_commands = {
    python = 'bpython -q',
    r = 'radian',
    lua = 'lua',
    julia = 'julia'
}

function ReplOperator(type, ...)
    -- vim.cmd("'[,']ReplSend") -- only supports whole lines sent
    -- vim.cmd("normal `]j") -- goto end of previous selection and then down since we are sending lines
    -- hacky alternative solution that works on any motion:
    if type == "line" then
        vim.cmd('normal `["ryV`]')
    else
        vim.cmd('normal `["ryv`]')
    end
    vim.cmd('wincmd l')
    vim.cmd('let @r .= "\r"')
    vim.cmd([[normal "rpG]])
    vim.cmd('wincmd h')
    vim.cmd("normal `]w") -- goto end of previous selection and then right
end

-- keymaps
vim.api.nvim_set_keymap('n', '<CR>', 'Operator("v:lua.ReplOperator")', {expr = true, noremap = false})
vim.api.nvim_set_keymap('n', '<CR><CR>', ':ReplSend<CR>j', {expr = false, noremap = false})
-- `> means go to mark named > which will be at the end of the previous selection.
vim.api.nvim_set_keymap('v', '<CR>', ':ReplSend<CR>`>w', {expr = false, noremap = false})
vim.api.nvim_set_keymap('n', '<localleader><CR>', ':silent ReplToggle<CR>', {expr = false, noremap = false})
vim.api.nvim_set_keymap('n', '<localleader>c', ':ReplClear<CR>', {expr = false, noremap = false})
vim.api.nvim_set_keymap('n', '<localleader><tab>', '<c-w>lA', {expr = false, noremap = false})
vim.api.nvim_set_keymap('v', '<localleader><tab>', 'y<c-w>lpA', {expr = false, noremap = false})


