#!/usr/bin/env lua
return {
    -- core behaviour
    "tpope/vim-repeat", -- change . to repeat last native command to last "full" command, which feels more natural.
    -- Y should yank to end of line which is consistent with other uppercase use, rather than yank whole line like yy which is for ancient vi compatibility.
    -- "tpope/vim-sensible",
    "tpope/vim-unimpaired", -- [e, ]e exchange line with above/below, ]<space> add newlines, more: https://github.com/tpope/
    "tpope/vim-characterize", -- unicode shown for ga
    -- add substitution functions to e.g. replace a word with clipboard content by writing siw
    "svermeulen/vim-subversive",
    -- "mg979/vim-visual-multi", -- multi cursor TODO https://github.com/mg979/vim-visual-multi/wiki/Quick-start
    "farmergreg/vim-lastplace", -- open file in last edited location
    "haya14busa/vim-asterisk", -- improvements to z* and visual *. See git for uses https://github.com/haya14busa/vim-asterisk
    -- "kana/vim-textobj-user", -- easily define custom textobjects such as i( and a( to select in/an \left( \right) block in latex
    -- TODO add from https://github.com/kana/vim-textobj-user and https://github.com/kana/vim-textobj-user/wiki
    {"glts/vim-textobj-comment", dependencies={"kana/vim-textobj-user"}}, -- not working?
    -- increment and decrement numbers, dates, color hex, even bool
    -- https://github.com/monaqa/dial.nvim
    {"monaqa/dial.nvim", config=function()
        local augend = require("dial.augend")
        require("dial.config").augends:register_group {
            default = {
                augend.integer.alias.decimal_int,
                augend.constant.alias.bool,    -- boolean value (true <-> false)
                augend.constant.new{ elements={"True", "False"}, word=true, cyclic=true, }, -- python
                augend.hexcolor.new{ case="lower", },
                augend.date.alias["%Y/%m/%d"],
                augend.date.alias["%Y-%m-%d"],
                augend.date.alias["%m/%d"],
                augend.date.alias["%H:%M"],
            },
        }

        vim.api.nvim_set_keymap("n", "<C-a>", require("dial.map").inc_normal(), {noremap = true})
        vim.api.nvim_set_keymap("n", "<C-x>", require("dial.map").dec_normal(), {noremap = true})
        vim.api.nvim_set_keymap("v", "<C-a>", require("dial.map").inc_visual(), {noremap = true})
        vim.api.nvim_set_keymap("v", "<C-x>", require("dial.map").dec_visual(), {noremap = true})
        vim.api.nvim_set_keymap("v", "g<C-a>", require("dial.map").inc_gvisual(), {noremap = true})
        vim.api.nvim_set_keymap("v", "g<C-x>", require("dial.map").dec_gvisual(), {noremap = true})
    end},
    
    "monkoose/matchparen.nvim", -- supposedly faster and less buggy version of neovim builtin (:h )matchparen which highlights matching parenthesis etc.
    
    -- jump to anywhere with \ + f or F or t or T
    {'ggandor/leap.nvim', config=function()
        leap = require 'leap'

        -- make s and S "unsafe", i.e. available immediately as a command
        -- add ' and ` as safe since it would be unlikely that I would want to jump to a mark right after a leap
        -- add [] as safe since it would be unlikely that I would want to jump with those after a leap
        leap.opts.safe_labels = {'f','n','u','t','/','`',"'",'[',']','F','N','L','H','M','U','G','T','?','Z'}

        -- mentioned in whichkey
        vim.keymap.set({'n', 'x', 'o'}, '\\f', '<Plug>(leap-forward-to)')
        vim.keymap.set({'n', 'x', 'o'}, '\\F', '<Plug>(leap-backward-to)')
        vim.keymap.set({'n', 'x', 'o'}, '\\t', '<Plug>(leap-forward-till)')
        vim.keymap.set({'n', 'x', 'o'}, '\\T', '<Plug>(leap-backward-till)')
    end},

    -- color
    -- when a hex or other color is defined, highlight the text with its color
    {"NvChad/nvim-colorizer.lua", opts={
        filetypes = {
            -- #RGB hex codes are too simple and may show up in julia in rare cases unnamed variables
            julia = { RGB = false, }
        }
    }}, 
    {"norcalli/nvim-base16.lua", dependencies={"norcalli/nvim.lua"}, config=function ()
        base16 = require'base16'
        theme = base16.theme_from_array(require("themes/gigavoltArray"))
        -- theme = base16.themes["unikitty-dark"]
        base16(theme, true)
    end},
    
    -- highlight letters for jumping with f/F/t/T
    {"unblevable/quick-scope", enabled=false, config=function()
        -- Trigger a highlight in the appropriate direction when pressing these keys:
        vim.g.qs_highlight_on_keys = {'f', 'F', 't', 'T'}
    end},

    "sakshamgupta05/vim-todo-highlight", -- highlight todos
    -- {"folke/twilight.nvim", config=function() require'twilight'.setup{dimming={alpha=0.5}, context=30} end}, -- dim code that isn't currently being edited with :Twilight.

    {
        "chrisgrieser/nvim-various-textobjs",
        opts = { useDefaultKeymaps = true },
    },

}

