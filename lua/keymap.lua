#!/usr/bin/env lua
-- from https://github.com/nvim-lua/kickstart.nvim/blob/master/init.lua

-- Keymaps for better default experience
-- See `:help vim.keymap.set()`
vim.keymap.set({ 'n', 'v' }, '<Space>', '<Nop>', { silent = true })

-- Remap for dealing with word wrap
vim.keymap.set('n', 'j', "v:count == 0 ? 'gj' : 'j'", { expr = true, silent = true })
vim.keymap.set('n', 'k', "v:count == 0 ? 'gk' : 'k'", { expr = true, silent = true })
vim.keymap.set('i', '<down>', [[v:count == 0 ? '<C-\><C-O>gj' : '<down>']], { expr = true, silent = true })
vim.keymap.set('i', '<up>', [[v:count == 0 ? '<C-\><C-O>gk' : '<up>']], { expr = true, silent = true })

-- signature help is naturally bound to ctrl+k in normal mode.
-- Access it in insert mode with a shift added, since ctrl+k is already bound 
-- to jump to next snippet.
vim.keymap.set("i", "<C-S-k>", vim.lsp.buf.signature_help)

-- I don't use command window much but I often press q: when I mean :q
vim.keymap.set('n', 'q;', 'q:')
vim.keymap.set('n', 'q:', ':q')

