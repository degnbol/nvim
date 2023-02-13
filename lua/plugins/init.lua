#!/usr/bin/env lua
return {
    -- core behaviour
    "tpope/vim-repeat", -- change . to repeat last native command to last "full" command, which feels more natural.
    -- problems: it doesn't have simple ds working at all.
    -- {"kylechui/nvim-surround", config=function() require"surround-conf" end}, -- press cs'" to change surrounding ' with ", ds' to delete surrounding ', ysiw) to surround word with ) and yss[ to surround line with [ ... ] (incl. spaces)
    -- "tpope/vim-sensible", -- Y should yank to end of line which is consistent with other uppercase use, rather than yank whole line like yy which is for ancient vi compatibility.
    "svermeulen/vim-subversive", -- add substitution functions to e.g. replace a word with clipboard content by writing siw
    {"gbprod/cutlass.nvim", config=function() require'cutlass-conf' end}, -- c(hange), d(elete) no longer copies, remapped in keymapping file so x will cut. Since we have added backspace and delete button support in normal mode there is no need for default x behavior
    "svermeulen/vim-yoink", -- yank history that you can cycle with c-n and c-p
    { 'ibhagwan/smartyank.nvim', config=function() require'smartyank-conf' end }, -- yank in tmux and over ssh
    -- "mg979/vim-visual-multi", -- multi cursor TODO https://github.com/mg979/vim-visual-multi/wiki/Quick-start
    -- {"moll/vim-bbye", config=function() require'bbye' end},
    "farmergreg/vim-lastplace", -- open file in last edited location
    {"inkarkat/vim-UnconditionalPaste", dependencies={'inkarkat/vim-ingo-library'}, config=function() require'unconditionalPaste' end}, -- lots of ways to paste using g{c,C,l,b}{,i}{p,P} and may others 
    {"AndrewRadev/whitespaste.vim", config=function() require'whitespaste' end}, -- paste without empty newlines
    "google/vim-searchindex", -- let's search result box show number of matches when there's >99 matches
    "haya14busa/vim-asterisk", -- improvements to z* and visual *. See git for uses https://github.com/haya14busa/vim-asterisk
    -- "kana/vim-textobj-user", -- easily define custom textobjects such as i( and a( to select in/an \left( \right) block in latex
    -- TODO add from https://github.com/kana/vim-textobj-user and https://github.com/kana/vim-textobj-user/wiki
    {"glts/vim-textobj-comment", dependencies={"kana/vim-textobj-user"}}, -- not working?
    {"AckslD/nvim-trevJ.lua", config=function() require'trevj-conf' end}, -- if it fails (gj), try revj (<leader>j[j] or motion INSIDE brackets)
    {"AckslD/nvim-revJ.lua", config=function() require'revj-conf' end, dependencies={'kana/vim-textobj-user', 'sgur/vim-textobj-parameter'}},
    {"monaqa/dial.nvim", config=function() require'dial-conf' end}, -- increment and decrement numbers, dates, color hex, even bool
    "monkoose/matchparen.nvim", -- supposedly faster and less buggy version of neovim builtin (:h )matchparen which highlights matching parenthesis etc.
    {'ggandor/leap.nvim', config=function() require"leap-conf" end}, -- jump to anywhere with \ + f or F or t or T
    
    -- color
    {"NvChad/nvim-colorizer.lua", config=function() require'colorizer-conf' end}, -- when a hex or other color is defined, highlight the text with its color
    {"norcalli/nvim-base16.lua", dependencies={"norcalli/nvim.lua"}},
    -- {"unblevable/quick-scope", config=function() require'quick-scope' end}, -- highlight letters for jumping with f/F/t/T
    
    {"DaikyXendo/nvim-material-icon", config=function() require'icons' end},
    
    "sakshamgupta05/vim-todo-highlight", -- highlight todos
    -- {"folke/twilight.nvim", config=function() require'twilight'.setup{dimming={alpha=0.5}, context=30} end}, -- dim code that isn't currently being edited with :Twilight.
}

