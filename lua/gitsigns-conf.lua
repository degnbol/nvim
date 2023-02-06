require("gitsigns").setup {
    signs = {
        add = {hl = "DiffAdd", text = "▌", numhl = "DiffAddNr"},
        change = {hl = "DiffChange", text = "▌", numhl = "DiffChangeNr"},
        delete = {hl = "DiffDelete", text = "_", numhl = "DiffDeleteNr"},
        topdelete = {hl = "DiffDelete", text = "‾", numhl = "DiffTopDeleteNr"},
        changedelete = {hl = "DiffChange", text = "~", numhl = "DiffChangeDeleteNr"}
    },
    -- highlight in signcolumn to the left of numbers, although :set signcolumn=no means this is suppressed.
    signcolumn = true, 
    numhl = true, -- highlight line number
    keymaps = {
        -- Default keymap options
        noremap = true,
        buffer = true,
        ["n ]h"] = {expr = true, '&diff ? \']h\' : \'<cmd>lua require"gitsigns".next_hunk()<CR>\''},
        ["n [h"] = {expr = true, '&diff ? \'[h\' : \'<cmd>lua require"gitsigns".prev_hunk()<CR>\''},
        ["n <leader>gs"] = '<cmd>lua require"gitsigns".stage_hunk()<CR>',
        ["n <leader>gu"] = '<cmd>lua require"gitsigns".undo_stage_hunk()<CR>',
        ["n <leader>gr"] = '<cmd>lua require"gitsigns".reset_hunk()<CR>',
        ["n <leader>gp"] = '<cmd>lua require"gitsigns".preview_hunk()<CR>',
        ["n <leader>gb"] = '<cmd>lua require"gitsigns".blame_line()<CR>'
    },
    watch_gitdir = { interval = 100 },
    sign_priority = 5,
    status_formatter = nil -- Use default
}
