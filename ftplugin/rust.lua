#!/usr/bin/env lua
vim.keymap.set('n', '<leader>cc',        ":Crun\n",   { desc="Cargo run" })
vim.keymap.set('n', '<leader><leader>r', ":Crun\n",   { desc="Cargo run" })
vim.keymap.set('n', '<leader><leader>b', ":Cbuild\n", { desc="Cargo build" })
vim.keymap.set('n', '<leader><leader>c', ":Ccheck\n", { desc="Cargo check" })
vim.keymap.set('n', '<leader><leader>C', ":Cclean\n", { desc="Cargo clean" })
vim.keymap.set('n', '<leader><leader>t', ":Ctest\n",  { desc="Cargo test" })
vim.keymap.set('n', '<leader><leader>u', ":Cupdate\n",{ desc="Cargo update" })
