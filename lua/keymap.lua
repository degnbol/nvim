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

-- switch to/from Danish æøå and to insert mode, which is convenient.
vim.keymap.set("n", "<leader>td", "i<C-^>", { desc="Danish (<C-^>)" })

-- small hack to remove excess whitespace possible since iw also captures 
-- whitespace under cursor.
vim.keymap.set("n", "di ", "ciw <Esc>", { desc="Delete excess whitespace" })

-- like =p but for substitution
vim.keymap.set("n", '=ss', 'ss=`]', {remap=true, silent=true, desc="Substitute+reindent"})
vim.keymap.set("n", '<leader>h', ':noh<CR>', {silent=true, desc="Clear highlights"})
vim.keymap.set("n", "<leader>rs", "<cmd>source $XDG_CONFIG_HOME/nvim/after/plugin/luasnip.lua<CR>", {desc="Reload snippets"})

-- <leader>b can be for buffer related things, (bufdel, but :bd is already 
-- fast) or bibliography things e.g. with papis but this is not added yet and 
-- could also just fit under latex (<leader>l) since I don't really use it in 
-- markdown.
