require "bufferline".setup {
    options = {
        offsets = {{filetype = "NvimTree", text = "", padding = 1}},
        indicator_icon = "", -- no line or other indicator in front of selected tab
        buffer_close_icon = "",
        modified_icon = "",
        close_icon = "",
        left_trunc_marker = "",
        right_trunc_marker = "",
        max_name_length = 25,
        max_prefix_length = 20,
        tab_size = 25,
        show_tab_indicators = false,
        enforce_regular_tabs = false, -- enforce same size tabs
        view = "multiwindow",
        show_buffer_close_icons = false,
        show_close_icon = false,
    }
}


-- move between tabs
vim.api.nvim_set_keymap("n", "<TAB>", [[<Cmd>BufferLineCycleNext<CR>]], {silent=true})
vim.api.nvim_set_keymap("n", "<S-TAB>", [[<Cmd>BufferLineCyclePrev<CR>]], {silent=true})
