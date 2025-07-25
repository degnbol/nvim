local hi = require "utils/highlights"
local map = require "utils/keymap"

-- add a line at the top with all the files open in the buffer
return {
    "akinsho/nvim-bufferline.lua",
    enabled = false,
    dependencies = { "nvim-tree/nvim-web-devicons" },
    init = function()
        for i = 1, 9 do
            map.n("<leader>" .. i, function() require('bufferline').go_to(i, true) end, "Buffer " .. i)
        end
        map.n("<leader>0", function() require('bufferline').go_to(-1, true) end, "Buffer last")
        map.n("<leader>bb", function() require('bufferline').go_to(vim.v.count, true) end, "Buffer N")
        map.n("<leader>b<left>", "<Cmd>BufferLineMovePrev<CR>", "Move left")
        map.n("<leader>b<right>", "<Cmd>BufferLineMoveNext<CR>", "Move right")
    end,
    opts = {
        options = {
            -- numbers = "ordinal", -- can be customized to show both ordinal number and buffer id (see readme)
            numbers = function(opts) return string.format('%s', opts.raise(opts.ordinal)) end,
            offsets = { { filetype = "NvimTree" } }, -- moves bufferline to the right when tree view is open
            indicator = { style = "none" },          -- no line or other indicator in front of selected tab
            -- separator_style = {"", ""}, -- remove separators
            separator_style = "slant",
            modified_icon = "",
            left_trunc_marker = "",
            right_trunc_marker = "",
            max_name_length = 20,
            max_prefix_length = 20,
            tab_size = 10,
            show_tab_indicators = false,    -- lsp diagnostics etc.
            enforce_regular_tabs = false,   -- enforce same size tabs
            always_show_bufferline = false, -- hide if only one file is open
            show_buffer_close_icons = false,
            show_close_icon = false,
            -- hide quickfix from bufferline:
            custom_filter = function(buf_number, buf_numbers)
                if vim.bo[buf_number].filetype ~= "qf" then
                    return true
                end
            end,
        },
        highlights = {
            buffer_selected = {
                italic = false,
            },
            numbers_selected = {
                italic = false,
            }
        },
    }
}
