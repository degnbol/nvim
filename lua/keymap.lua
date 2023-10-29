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

-- I don't use command window much but I often press q: when I mean :q
vim.keymap.set('n', 'q;', 'q:')
vim.keymap.set('n', 'q:', ':q')

-- switch to/from Danish æøå and to insert mode, which is convenient.
-- remap in order to utilise the remapped <C-^> which updates the cmp dictionary
vim.keymap.set("n", "<leader>td", "i<C-^>", { remap=true, desc="Danish (<C-^>)" })

-- small hack to remove excess whitespace possible since iw also captures 
-- whitespace under cursor.
vim.keymap.set("n", "di ", "ciw <Esc>", { desc="Delete excess whitespace" })

-- like =p but for substitution
vim.keymap.set("n", '=ss', 'ss=`]', {remap=true, silent=true, desc="Substitute+reindent"})
vim.keymap.set("n", "<leader>Sr", "<cmd>source $XDG_CONFIG_HOME/nvim/after/plugin/luasnip.lua<CR>", {desc="Reload snippets"})

-- <leader>b can be for buffer related things, (bufdel, but :bd is already 
-- fast) or bibliography things e.g. with papis but this is not added yet and 
-- could also just fit under latex (<leader>l) since I don't really use it in 
-- markdown.

-- use the following two commands to enable spelling
-- setlocal spell
-- set spelllang=en_us
-- ctrl+s to fix last misspelled word
-- credit: https://castel.dev/post/lecture-notes-1/
-- with git: https://github.com/gillescastel/latex-snippets
-- "i_ctrl-g u" = something about undo history, probably to keep the edit part of the insert edit?
-- [s = goto last misspelled, 1z= = replace with first spell suggestion. 
vim.keymap.set('i', '<C-s>', [[<c-g>u<Esc>[s1z=`]i<c-g>u]], { desc="Spell correct closest" })
vim.keymap.set('n', '<C-s>', [=[[s1z=``]=], { desc="Spell correct closest" })

-- not the most elegant but it works.
-- LeftMouse to move cursor to pressed location.
-- Then set @/ to the current word (\< and \> are to search strictly).
-- Then enable hlsearch. This is all a way to search without going to the next 
-- match (if we just pressed * for instance)
vim.keymap.set('n', '<RightMouse>',
    [[<LeftMouse>:let @/='\<'.expand('<cword>').'\>'|set hlsearch<CR>]],
    { silent=true, desc="Search pressed word" }
)

-- fallback search replace if both treesitter and LSP are not attached.
vim.keymap.set('n', '<leader>rn', [[:%s/<C-r><C-w>/]], { desc="Search/replace cword" })
-- use a selection that isn't a perfect cword, or just to use the simple search/replace when LSP is attached etc.
vim.keymap.set('x', '<leader>rn', [["ry:%s/<C-r>r/]], { desc="Search/replace" })

vim.keymap.set('i', '<C-6>',
    [[<C-^><C-\><C-o>:doautocmd User ToggleDansk<CR>]],
    { silent=true, desc="Toggle dansk" }
)
-- Ins key was a bit useless just doing what i does so let's make it a language switch insertion:
vim.keymap.set('n', '<Ins>', 'i<C-6>', { silent=true, desc="Insert + Toggle dansk" })
-- <sa = my keybind for enable setting autoformat
-- gww = autoformat line
-- 0   = goto column 0, so we scroll all the way back to the right
-- gi  = go to last insert location and enter insert mode. Works even with the change to the line.
vim.keymap.set('i', '<C-S-A>', "<Esc><sagww0gi", { remap=true, desc="Enable autoformat and apply it." })

vim.keymap.set('n', '<leader>bd', ":bd<CR>", { desc="Delete" })
vim.keymap.set('n', '<leader>bD', ":bd!<CR>", { desc="Delete!" })

