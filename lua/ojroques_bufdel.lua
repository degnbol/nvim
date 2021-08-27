-- close tab. :BufDel from https://github.com/ojroques/nvim-bufdel instead of built-in bwipeout.
vim.api.nvim_set_keymap("n", "<leader>x", [[<Cmd>BufDel<CR>]], {silent=true})
