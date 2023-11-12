-- add a line at the top with all the files open in the buffer
return {
    "akinsho/nvim-bufferline.lua",
    dependencies={"DaikyXendo/nvim-material-icon"},
    init = function ()
        -- consider using [b and ]b instead. Then we could clear up tab in normal 
        -- mode for something related to completion.
        -- vim.keymap.set("n", "<TAB>", "<Cmd>BufferLineCycleNext<CR>", {desc="Next buffer"})
        -- vim.keymap.set("n", "<S-TAB>", "<Cmd>BufferLineCyclePrev<CR>", {desc="Previous buffer"})
        vim.keymap.set("n", "<leader>1",       "<Cmd>BufferLineGoToBuffer 1<CR>", {desc="Buffer 1"})
        vim.keymap.set("n", "<leader>2",       "<Cmd>BufferLineGoToBuffer 2<CR>", {desc="Buffer 2"})
        vim.keymap.set("n", "<leader>3",       "<Cmd>BufferLineGoToBuffer 3<CR>", {desc="Buffer 3"})
        vim.keymap.set("n", "<leader>4",       "<Cmd>BufferLineGoToBuffer 4<CR>", {desc="Buffer 4"})
        vim.keymap.set("n", "<leader>5",       "<Cmd>BufferLineGoToBuffer 5<CR>", {desc="Buffer 5"})
        vim.keymap.set("n", "<leader>6",       "<Cmd>BufferLineGoToBuffer 6<CR>", {desc="Buffer 6"})
        vim.keymap.set("n", "<leader>7",       "<Cmd>BufferLineGoToBuffer 7<CR>", {desc="Buffer 7"})
        vim.keymap.set("n", "<leader>8",       "<Cmd>BufferLineGoToBuffer 8<CR>", {desc="Buffer 8"})
        vim.keymap.set("n", "<leader>9",       "<Cmd>BufferLineGoToBuffer 9<CR>", {desc="Buffer 9"})
        vim.keymap.set("n", "<leader>b<left>", "<Cmd>BufferLineMovePrev<CR>",     {desc="Move left"})
        vim.keymap.set("n", "<leader>b<right>","<Cmd>BufferLineMoveNext<CR>",     {desc="Move right"})
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

