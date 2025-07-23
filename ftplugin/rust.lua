
vim.keymap.set('n', '<leader>cc',     ":Crun\n",   { buffer=true, desc="Cargo run" })
vim.keymap.set('n', '<LocalLeader>r', ":Crun\n",   { buffer=true, desc="Cargo run" })
vim.keymap.set('n', '<LocalLeader>b', ":Cbuild\n", { buffer=true, desc="Cargo build" })
vim.keymap.set('n', '<LocalLeader>c', ":Ccheck\n", { buffer=true, desc="Cargo check" })
vim.keymap.set('n', '<LocalLeader>C', ":Cclean\n", { buffer=true, desc="Cargo clean" })
vim.keymap.set('n', '<LocalLeader>t', ":Ctest\n",  { buffer=true, desc="Cargo test" })
vim.keymap.set('n', '<LocalLeader>u', ":Cupdate\n",{ buffer=true, desc="Cargo update" })
