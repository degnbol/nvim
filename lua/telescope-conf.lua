#!/usr/bin/env lua
require'telescope'.setup {
    defaults = {
        layout_strategy = "vertical",
        prompt_prefix = '',
        selection_caret = '',
        entry_prefix = '',
        multi_icon = '',
    },
    extensions = {
        dash = {
            file_type_keywords = {
                python = {"python", "numpy", "scipy", "pandas"}
            }
        }
    }
}

-- load the native fzf as recommended
require'telescope'.load_extension('fzf')

