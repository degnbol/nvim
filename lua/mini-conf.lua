#!/usr/bin/env lua
require('mini.surround').setup {
    mappings = {
    add = 'ys', -- Add surrounding in Normal and Visual modes
    delete = 'ds', -- Delete surrounding
    find = '', -- Find surrounding (to the right)
    find_left = '', -- Find surrounding (to the left)
    highlight = '', -- Highlight surrounding
    replace = 'cs', -- Replace surrounding
    update_n_lines = '', -- Update `n_lines`
  },
  custom_surroundings = {
      ['d'] = { output = { left = '', right = '' } },
      -- normally used for tags such as <div> </div> which is cool but I don't use them so t for triple
      ['t'] = {
          -- figured out these patterns from the [[ ]] example in :h MiniSurround.config
          input = { find = '%"%"%".-%"%"%"', extract = '^(...).*(...)$' },
          output = { left = '"""', right = '"""' }
      },
   }
}

-- Remap adding surrounding to Visual mode selection
vim.api.nvim_del_keymap('x', 'ys')
vim.api.nvim_set_keymap('x', 'S', [[:<C-u>lua MiniSurround.add('visual')<CR>]], { noremap = true })

-- Make special mapping for "add surrounding for line"
vim.api.nvim_set_keymap('n', 'yss', 'ys_', { noremap = false })

