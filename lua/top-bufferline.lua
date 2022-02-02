require "bufferline".setup {
    options = {
        numbers = "ordinal", -- can be customized to show both ordinal number and buffer id (see readme)
        offsets = {{filetype = "NvimTree"}}, -- moves bufferline to the right when tree view is open
        indicator_icon = " ", -- no line or other indicator in front of selected tab
        buffer_close_icon = "",
        modified_icon = "",
        close_icon = "",
        left_trunc_marker = "",
        right_trunc_marker = "",
        max_name_length = 25,
        max_prefix_length = 20,
        tab_size = 25,
        show_tab_indicators = false, -- lsp diagnostics etc.
        enforce_regular_tabs = false, -- enforce same size tabs
        always_show_bufferline = false, -- hide if only one file is open
        show_buffer_close_icons = false,
        show_close_icon = false,
    }
}


-- move between tabs
vim.api.nvim_set_keymap("n", "<TAB>", [[<Cmd>BufferLineCycleNext<CR>]], {silent=true})
vim.api.nvim_set_keymap("n", "<S-TAB>", [[<Cmd>BufferLineCyclePrev<CR>]], {silent=true})
vim.api.nvim_set_keymap("n", "<leader>1", [[<Cmd>BufferLineGoToBuffer 1<CR>]], {silent=true})
vim.api.nvim_set_keymap("n", "<leader>2", [[<Cmd>BufferLineGoToBuffer 2<CR>]], {silent=true})
vim.api.nvim_set_keymap("n", "<leader>3", [[<Cmd>BufferLineGoToBuffer 3<CR>]], {silent=true})
vim.api.nvim_set_keymap("n", "<leader>4", [[<Cmd>BufferLineGoToBuffer 4<CR>]], {silent=true})
vim.api.nvim_set_keymap("n", "<leader>5", [[<Cmd>BufferLineGoToBuffer 5<CR>]], {silent=true})
vim.api.nvim_set_keymap("n", "<leader>6", [[<Cmd>BufferLineGoToBuffer 6<CR>]], {silent=true})
vim.api.nvim_set_keymap("n", "<leader>7", [[<Cmd>BufferLineGoToBuffer 7<CR>]], {silent=true})
vim.api.nvim_set_keymap("n", "<leader>8", [[<Cmd>BufferLineGoToBuffer 8<CR>]], {silent=true})
vim.api.nvim_set_keymap("n", "<leader>9", [[<Cmd>BufferLineGoToBuffer 9<CR>]], {silent=true})

