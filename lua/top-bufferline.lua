require "bufferline".setup {
    options = {
        -- numbers = "ordinal", -- can be customized to show both ordinal number and buffer id (see readme)
        numbers = function(opts) return string.format('%s', opts.raise(opts.ordinal)) end,
        offsets = {{filetype = "NvimTree"}}, -- moves bufferline to the right when tree view is open
        indicator_icon = " ", -- no line or other indicator in front of selected tab
        separator_style = {"", ""}, -- remove separators
        buffer_close_icon = "",
        modified_icon = "",
        close_icon = "",
        left_trunc_marker = "",
        right_trunc_marker = "",
        max_name_length = 20,
        max_prefix_length = 20,
        tab_size = 10,
        show_tab_indicators = false, -- lsp diagnostics etc.
        enforce_regular_tabs = false, -- enforce same size tabs
        always_show_bufferline = false, -- hide if only one file is open
        show_buffer_close_icons = false,
        show_close_icon = false,
    }
}

