#!/usr/bin/env lua
return {
    -- core behaviour
    "tpope/vim-repeat", -- change . to repeat last native command to last "full" command, which feels more natural.
    -- Y should yank to end of line which is consistent with other uppercase use, rather than yank whole line like yy which is for ancient vi compatibility.
    -- "tpope/vim-sensible",
    -- add substitution functions to e.g. replace a word with clipboard content by writing siw
    "svermeulen/vim-subversive",
    -- "mg979/vim-visual-multi", -- multi cursor TODO https://github.com/mg979/vim-visual-multi/wiki/Quick-start
    "farmergreg/vim-lastplace", -- open file in last edited location
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
    -- when a hex or other color is defined, highlight the text with its color
    {"NvChad/nvim-colorizer.lua", opts={
        filetypes = {
            -- #RGB hex codes are too simple and may show up in julia in rare cases unnamed variables
            julia = { RGB = false, }
        }
    }}, 
    {"norcalli/nvim-base16.lua", dependencies={"norcalli/nvim.lua"}},
    {"unblevable/quick-scope", enabled=false, config=function()
        -- Trigger a highlight in the appropriate direction when pressing these keys:
        vim.g.qs_highlight_on_keys = {'f', 'F', 't', 'T'}
    end}, -- highlight letters for jumping with f/F/t/T

    "sakshamgupta05/vim-todo-highlight", -- highlight todos
    -- {"folke/twilight.nvim", config=function() require'twilight'.setup{dimming={alpha=0.5}, context=30} end}, -- dim code that isn't currently being edited with :Twilight.
}

