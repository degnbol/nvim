#!/usr/bin/env lua
return {
    -- treesitter
    {'nvim-treesitter/nvim-treesitter', build=':TSUpdate', config=function() require'treesitter' end}, -- language coloring and ensuring of installation
    {"nvim-treesitter/nvim-treesitter-refactor", dependencies='nvim-treesitter/nvim-treesitter', config=function() require'treesitter-refactor' end}, -- refactor
    {"nvim-treesitter/nvim-treesitter-textobjects", dependencies='nvim-treesitter/nvim-treesitter', config=function() require'treesitter-textobjects' end}, -- selecting, moving functions etc.
    {"RRethy/nvim-treesitter-textsubjects", dependencies='nvim-treesitter/nvim-treesitter', config=function() require'treesitter-textsubjects' end}, -- in vis mode use . , ; i; to select based on treesitter 
    -- "romgrk/nvim-treesitter-context", -- show the "context" at the top line, i.e. function name when in a function
    -- error for julia tree-sitter:
    -- {"andymass/vim-matchup", dependencies='nvim-treesitter/nvim-treesitter', config=function() require'matchup' end}, -- % jumps between matching coding blocks, not just single chars.
    {"p00f/nvim-ts-rainbow", dependencies='nvim-treesitter/nvim-treesitter', config=function() require'treesitter-rainbow' end}, -- tree sitter based rainbow color parenthesis to easily see the matching
    {"nvim-treesitter/playground"}, -- provides :TSHighlightCapturesUnderCursor to see highlight groups for a word under the cursor, TSPlaygroundToggle. "a" for hidden, "o" for scratch edit scheme.
}

