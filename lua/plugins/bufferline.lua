-- add a line at the top with all the files open in the buffer
return {
    "akinsho/nvim-bufferline.lua",
    dependencies={"nvim-tree/nvim-web-devicons"},
    init = function ()
        vim.keymap.set("n", "<leader>1",  function () require('bufferline').go_to(1,true) end, {desc="Buffer 1"})
        vim.keymap.set("n", "<leader>2",  function () require('bufferline').go_to(2,true) end, {desc="Buffer 2"})
        vim.keymap.set("n", "<leader>3",  function () require('bufferline').go_to(3,true) end, {desc="Buffer 3"})
        vim.keymap.set("n", "<leader>4",  function () require('bufferline').go_to(4,true) end, {desc="Buffer 4"})
        vim.keymap.set("n", "<leader>5",  function () require('bufferline').go_to(5,true) end, {desc="Buffer 5"})
        vim.keymap.set("n", "<leader>6",  function () require('bufferline').go_to(6,true) end, {desc="Buffer 6"})
        vim.keymap.set("n", "<leader>7",  function () require('bufferline').go_to(7,true) end, {desc="Buffer 7"})
        vim.keymap.set("n", "<leader>8",  function () require('bufferline').go_to(8,true) end, {desc="Buffer 8"})
        vim.keymap.set("n", "<leader>9",  function () require('bufferline').go_to(9,true) end, {desc="Buffer 9"})
        vim.keymap.set("n", "<leader>bb", function () require('bufferline').go_to(vim.v.count,true) end, {desc="Buffer N"})
        vim.keymap.set("n", "<leader>b<left>", "<Cmd>BufferLineMovePrev<CR>", {desc="Move left"})
        vim.keymap.set("n", "<leader>b<right>","<Cmd>BufferLineMoveNext<CR>", {desc="Move right"})
    end,
    opts = {
        options = {
            -- numbers = "ordinal", -- can be customized to show both ordinal number and buffer id (see readme)
            numbers = function(opts) return string.format('%s', opts.raise(opts.ordinal)) end,
            offsets = {{filetype = "NvimTree"}}, -- moves bufferline to the right when tree view is open
            indicator = {style="none"}, -- no line or other indicator in front of selected tab
            -- separator_style = {"", ""}, -- remove separators
            separator_style = "slant",
            modified_icon = "",
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
}

