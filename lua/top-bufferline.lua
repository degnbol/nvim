require "bufferline".setup {
    options = {
        offsets = {{filetype = "NvimTree", text = "", padding = 1}},
        buffer_close_icon = "",
        modified_icon = "",
        close_icon = "",
        left_trunc_marker = "",
        right_trunc_marker = "",
        max_name_length = 20,
        max_prefix_length = 15,
        tab_size = 20,
        show_tab_indicators = false,
        enforce_regular_tabs = false, -- enforce same size tabs
        view = "multiwindow",
        show_buffer_close_icons = false,
        show_close_icon = false,
        mappings = "true"
    }
}

local opt = {silent = true}
local map = vim.api.nvim_set_keymap
vim.g.mapleader = " "

-- MAPPINGS
map("n", "<leader>x", [[<Cmd>BufDel<CR>]], opt) -- close tab. :BufDel from https://github.com/ojroques/nvim-bufdel instead of built-in bwipeout.

-- move between tabs
map("n", "<TAB>", [[<Cmd>BufferLineCycleNext<CR>]], opt)
map("n", "<S-TAB>", [[<Cmd>BufferLineCyclePrev<CR>]], opt)
