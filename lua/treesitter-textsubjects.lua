#!/usr/bin/env lua
-- while in visual mode these keybindings will change what is selected
require('nvim-treesitter.configs').setup {
    textsubjects = {
        enable = true,
        -- optional keymap to select the previous selection, which effectively decreases the incremental smart selection
        prev_selection = ',',
        keymaps = {
            -- this will do incremental expand select
            ['.'] = 'textsubjects-smart',
            -- treesitter based container will be selected with v;
            [';'] = 'textsubjects-container-outer',
            -- inside of treesitter based container will be selected with vi;
            ['i;'] = 'textsubjects-container-inner',
        },
    },
}
