-- assume using typst.vim
vim.keymap.set('n', '<leader>cc', "<Cmd>TypstWatch<CR>", { desc="Compile continuously", buffer=true })
vim.keymap.set('n', '<leader>cv', "!open -a sioyek<CR>", { desc="Compile view", buffer=true, silent=true, })
