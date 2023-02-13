#!/usr/bin/env lua
return {
    -- git
    "tpope/vim-fugitive", -- git
    {"lewis6991/gitsigns.nvim", dependencies={'nvim-lua/plenary.nvim'}, config=function() require'gitsigns-conf' end}, -- git decoration to the left
    {"rhysd/conflict-marker.vim"}, -- highlight git conflicts, jump with [x and ]x, resolve by keeping none (cn), theirs (ct), our (co), both (cb), or both reverse (cB)
    { 'sindrets/diffview.nvim', dependencies={'nvim-lua/plenary.nvim', 'nvim-tree/nvim-web-devicons'}, config=function() require'diffview-conf' end }, -- :DiffviewOpen and other commands for seeing git diff and git history for files.
}

