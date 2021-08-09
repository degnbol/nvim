local iron = require("iron")

iron.core.set_config {
  preferred = {
    python = "ipython"
  },
  -- open to the right rather than the default left
  repl_open_cmd = "rightbelow vertical split",
}


-- keymaps
-- TODO make a function where if the python line starts with def then send
-- V]M<Plug>(iron-visual-send)`>. Julia it should instead look for function
-- anywhere in the line and call V][<Plug>(iron-visual-send)`>
vim.api.nvim_set_keymap('n', '<localleader><CR>', '<Plug>(iron-send-line)j', {expr = false, noremap = false})
-- After sending to visual the cursor jumps to the start of the selection. 
-- `> means go to mark named > which will be at the end of the previous selection.
vim.api.nvim_set_keymap('v', '<localleader><CR>', '<Plug>(iron-visual-send)`>', {expr = false, noremap = false})
vim.api.nvim_set_keymap('n', '<localleader>r', '<Plug>(iron-repeat-cmd)', {expr = false, noremap = false})
vim.api.nvim_set_keymap('n', '<localleader>i', '<Plug>(iron-interrupt)', {expr = false, noremap = false})
-- Exit -> focus -> terminal mode -> newline to close window pane.
vim.api.nvim_set_keymap('n', '<localleader>q', '<Plug>(iron-exit):IronFocus<CR>A<CR>', {expr = false, noremap = false})
vim.api.nvim_set_keymap('n', '<localleader>c', '<Plug>(iron-clear)', {expr = false, noremap = false})
vim.api.nvim_set_keymap('n', '<localleader><tab>', ':IronFocus<CR>A', {expr = false, noremap = false})
vim.api.nvim_set_keymap('v', '<localleader><tab>', 'y:IronFocus<CR>pA', {expr = false, noremap = false})

-- Default behaviour for adding a newline to the console:
-- nmap c<CR> <Plug>(iron-cr)
-- other defaults on
-- https://github.com/hkupty/iron.nvim/blob/master/doc/iron.txt



