#!/usr/bin/env lua
return {{
    'echasnovski/mini.nvim',
    version = false,
    priority = 100, -- higher than default 50 to allow mini parts loaded in other files.
    config = function ()

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

-- highlight hex colors and todos, notes etc
-- consider this or https://github.com/folke/paint.nvim if you want to highlight custom things.
local hipatterns = require'mini.hipatterns'
hipatterns.setup {
    highlighters = {
        -- Highlight standalone 'FIXME', 'HACK', 'TODO', 'NOTE'
        fixme = { pattern = '%f[%w]()FIXME()%f[%W]', group = 'MiniHipatternsFixme' },
        hack  = { pattern = '%f[%w]()HACK()%f[%W]',  group = 'MiniHipatternsHack'  },
        todo  = { pattern = '%f[%w]()TODO()%f[%W]',  group = 'MiniHipatternsTodo'  },
        note  = { pattern = '%f[%w]()NOTE()%f[%W]',  group = 'MiniHipatternsNote'  },

        -- Highlight hex color strings (`#rrggbb`) using that color
        hex_color = hipatterns.gen_highlighter.hex_color(),
    }
}


local clue = require'mini.clue'

clue.setup({
  triggers = {
    { mode = 'n', keys = '<Leader>' },
    { mode = 'x', keys = '<Leader>' },
    { mode = 'n', keys = 'g' },
    { mode = 'x', keys = 'g' },
    { mode = 'n', keys = 'y' },
    { mode = 'x', keys = 'y' },
    { mode = 'n', keys = 'z' },
    { mode = 'x', keys = 'z' },
    { mode = 'n', keys = ']' },
    { mode = 'n', keys = '[' },
    -- Built-in completion
    { mode = 'i', keys = '<C-x>' },
    -- Marks
    { mode = 'n', keys = "'" },
    { mode = 'n', keys = '`' },
    { mode = 'x', keys = "'" },
    { mode = 'x', keys = '`' },
    -- Registers
    { mode = 'n', keys = '"' },
    { mode = 'x', keys = '"' },
    { mode = 'i', keys = '<C-r>' },
    { mode = 'c', keys = '<C-r>' },
    -- Window
    { mode = 'n', keys = '<C-w>' },
  },

  clues = {
    -- use e.g. postkeys='<C-w>' to make a submode. 
    { mode = 'n', keys = '<leader>c', desc = "Choose|Color|Code" },
    { mode = 'n', keys = '<leader>d', desc = "Diagnostics|Def peek" },
    { mode = 'n', keys = '<leader>f', desc = "Find|File" },
    { mode = 'n', keys = '<leader>g', desc = "Git" },
    { mode = 'n', keys = '<leader>l', desc = "LaTeX" },
    { mode = 'n', keys = '<leader>p', desc = "Paste" },
    { mode = 'n', keys = '<leader>P', desc = "Paste before" },
    { mode = 'n', keys = '<leader>r', desc = "Rename" },
    { mode = 'n', keys = '<leader>t', desc = "Toggle|Term" },
    { mode = 'n', keys = '<leader>w', desc = "Workspace" },
    { mode = 'n', keys = '<leader>wl', desc = "List" },
    { mode = 'n', keys = 'gc', desc = "Comment" },
    { mode = 'n', keys = 'gcc', desc = "Line" },
    -- Enhance this by adding descriptions for <Leader> mapping groups
    clue.gen_clues.builtin_completion(),
    clue.gen_clues.g(),
    clue.gen_clues.marks(),
    clue.gen_clues.registers(),
    clue.gen_clues.windows{
        -- making submodes of the window stuff means we can increase hight with <C-w>++ instead of <C-w>+<C-w>+
        submode_move = true,
        submode_navigate = true,
        submode_resize = true,
    },
    clue.gen_clues.z(),
  },

  window = {
      config = {
          width = "auto",
      },
  },
})


end},
}
