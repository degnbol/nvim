#!/usr/bin/env lua
leap = require 'leap'

-- make s and S "unsafe", i.e. available immediately as a command
-- add ' and ` as safe since it would be unlikely that I would want to jump to a mark right after a leap
-- add [] as safe since it would be unlikely that I would want to jump with those after a leap
leap.opts.safe_labels = {'f','n','u','t','/','`',"'",'[',']','F','N','L','H','M','U','G','T','?','Z'}

-- mentioned in whichkey
vim.keymap.set({'n', 'x', 'o'}, '\\f', '<Plug>(leap-forward-to)')
vim.keymap.set({'n', 'x', 'o'}, '\\F', '<Plug>(leap-backward-to)')
vim.keymap.set({'n', 'x', 'o'}, '\\t', '<Plug>(leap-forward-till)')
vim.keymap.set({'n', 'x', 'o'}, '\\T', '<Plug>(leap-backward-till)')
