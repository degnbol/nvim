#!/usr/bin/env lua
require("nvim-treesitter.configs").setup {
    -- for the p00f/nvim-ts-rainbow plugin
    rainbow = {
        enable = false, -- update broke this plugin
        extended_mode = true, -- Highlight also non-parentheses delimiters, boolean or table: lang -> boolean
        max_file_lines = 1000, -- Do not enable for files with more than 1000 lines, int
    },
}
