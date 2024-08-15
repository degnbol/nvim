#!/usr/bin/env lua
vim.keymap.set('n', '<leader>cc', ":Cargo run\n", { desc="Cargo run" })
vim.keymap.set('n', '<leader>cb', ":Cargo build\n", { desc="Cargo build" })
vim.keymap.set('n', '<leader>cC', ":Cargo check\n", { desc="Cargo check" })
