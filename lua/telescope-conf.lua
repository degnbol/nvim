
require'telescope'.setup {
    defaults = {
        layout_strategy = "vertical",
        prompt_prefix = '',
        selection_caret = '',
        entry_prefix = '',
        multi_icon = '',
    },
}

-- load the native fzf as recommended. Not working.
-- require'telescope'.load_extension('fzf')
