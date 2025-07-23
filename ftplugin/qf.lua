

-- quickfix specific keymaps.
-- Don't change buffer in quickfix window, first go to other win.
vim.keymap.set('n', '[b', "<C-w><C-w><Cmd>lua MiniBracketed.buffer('forward')<CR>", { desc="Change buf outside qf", buffer=true, silent=true })
vim.keymap.set('n', ']b', "<C-w><C-w><Cmd>lua MiniBracketed.buffer('forward')<CR>", { desc="Change buf outside qf", buffer=true, silent=true })

