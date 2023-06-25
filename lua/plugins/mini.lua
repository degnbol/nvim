#!/usr/bin/env lua
return {{
    'echasnovski/mini.nvim',
    version=false,
    priority = 100, -- higher than default 50 to allow mini parts loaded in other files.
    config=function()

require('mini.bracketed').setup {
    -- ]i to go to more indented region.
    indent = { suffix = 'i', options = { change_type="more" } },
    -- compliment with using [[, ]], [], ][ to jump to less indented region and to next region
}

-- mini_indentscope = require('mini.indentscope')
-- mini_indentscope.setup {
--     draw = { animation = mini_indentscope.gen_animation.none() },
--     symbol = "▏",
-- }


require'mini.surround'.setup {
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
      -- [d]elete
      ['d'] = { output = { left = '', right = '' } },
      -- normally used for tags such as <div> </div> which is cool but I don't use them so t for triple
      ['t'] = {
          -- figured out these patterns from the [[ ]] example in :h MiniSurround.config
          input = { find = '%"%"%".-%"%"%"', extract = '^(...).*(...)$' },
          output = { left = '"""', right = '"""' }
      },
      -- Taken from `:h MiniSurround.config` for use in lua.
      -- I never really use the open versions of brackets that adds space so we override '['.
      ['['] = {
        input = { '%[%[().-()%]%]' },
        output = { left = '[[', right = ']]' },
      },
      -- for latex ``...'', ~ is chosen since it is shift+`
      ['~'] = {
        input = { "``().-()''" },
        output = { left = '``', right = "''" },
      },
   },
   -- Number of lines within which surrounding is searched
  n_lines = 100,
}

-- Remap adding surrounding to Visual mode selection
vim.api.nvim_del_keymap('x', 'ys')
vim.api.nvim_set_keymap('x', 'S', [[:<C-u>lua MiniSurround.add('visual')<CR>]], { noremap = true })

-- Make special mapping for "add surrounding for line"
vim.api.nvim_set_keymap('n', 'yss', 'ys_', { noremap = false })

end},
}
