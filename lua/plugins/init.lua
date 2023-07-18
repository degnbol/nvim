#!/usr/bin/env lua
return {
    -- core behaviour
    "tpope/vim-repeat", -- change . to repeat last native command to last "full" command, which feels more natural.
    -- Y should yank to end of line which is consistent with other uppercase use, rather than yank whole line like yy which is for ancient vi compatibility.
    -- "tpope/vim-sensible",
    -- TODO: list these with whichkey?
    "tpope/vim-unimpaired", -- [e, ]e exchange line with above/below, ]<space> add newlines, and more
    {
        -- unicode shown for ga, we have go-align on that so have to use gA
        "tpope/vim-characterize",
        keys="gA",
        init = function () vim.keymap.set("n", "gA", "ga") end
    },
    -- add substitution functions to e.g. replace a word with clipboard content by writing siw
    "svermeulen/vim-subversive",
    -- "mg979/vim-visual-multi", -- multi cursor TODO https://github.com/mg979/vim-visual-multi/wiki/Quick-start
    "farmergreg/vim-lastplace", -- open file in last edited location
    "haya14busa/vim-asterisk", -- improvements to z* and visual *. See git for uses https://github.com/haya14busa/vim-asterisk
    -- "gioele/vim-autoswap",
    -- easily define custom textobjects
    -- some config in ftplugin/tex.lua
    "kana/vim-textobj-user", 
    -- TODO: add from https://github.com/kana/vim-textobj-user and https://github.com/kana/vim-textobj-user/wiki
    {"glts/vim-textobj-comment", dependencies={"kana/vim-textobj-user"}}, -- not working?
    {
        -- https://github.com/chrisgrieser/nvim-various-textobjs
        "chrisgrieser/nvim-various-textobjs",
        -- opts = { useDefaultKeymaps = true },
        config = function ()
            
            require("various-textobjs").setup {
                useDefaultKeymaps = true,
                -- disable some default keymaps, e.g. { "ai", "ii" }
                disabledKeymaps = {"ak", "ik", "av", "iv"},
            }
            
            local open = require("utils/init").open
            
            vim.keymap.set("n", "gx", function()
                -- go to github for plugin easily.
                -- First check if we are editing a file read by lazy.nvim
                if vim.api.nvim_buf_get_name(0):find('nvim/lua/plugins/') then
                    -- get first string on the line, assumes we don't list multiple plugins on one line.
                    local line = vim.api.nvim_get_current_line()
                    local repo = line:match([["([%w%p]+/[%w%p]+)"]])
                    repo = repo or line:match([['([%w%p]+/[%w%p]+)']])
                    if repo then
                        return open("https://github.com/" .. repo)
                    end
                end
                
                -- visually select URL
                require("various-textobjs").url()

                -- plugin only switches to visual mode when textobj found
                local foundURL = vim.fn.mode():find("v")

                -- retrieve URL with the z-register as intermediary
                vim.cmd.normal { '"zy', bang = true }
                local url = vim.fn.getreg("z")
                open(url)

            end, { desc = "Smart URL Opener" })
        end
    },
    -- increment and decrement numbers, dates, color hex, even bool
    {"monaqa/dial.nvim",
    keys = {"<C-a>", "<C-x>", "g<C-a>", "g<C-x>"},
    config=function()
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

        vim.api.nvim_set_keymap("n", "<C-a>",  require("dial.map").inc_normal(), {})
        vim.api.nvim_set_keymap("n", "<C-x>",  require("dial.map").dec_normal(), {})
        vim.api.nvim_set_keymap("v", "<C-a>",  require("dial.map").inc_visual(), {})
        vim.api.nvim_set_keymap("v", "<C-x>",  require("dial.map").dec_visual(), {})
        vim.api.nvim_set_keymap("v", "g<C-a>", require("dial.map").inc_gvisual(), {})
        vim.api.nvim_set_keymap("v", "g<C-x>", require("dial.map").dec_gvisual(), {})
    end},
    
    "lervag/file-line", -- open a file on a line with vi filepath:linenumber

    -- supposedly faster and less buggy version of neovim builtin (:h )matchparen which highlights matching parenthesis etc.
    "monkoose/matchparen.nvim",

    -- jump to anywhere with \ + f or F or t or T (set in whichkey)
    {'ggandor/leap.nvim',
    keys = "\\",
    config=function()
        -- make s and S "unsafe", i.e. available immediately as a command
        -- add ' and ` as safe since it would be unlikely that I would want to jump to a mark right after a leap
        -- add [] as safe since it would be unlikely that I would want to jump with those after a leap
        require'leap'.opts.safe_labels = {'f','n','u','t','/','`',"'",'[',']','F','N','L','H','M','U','G','T','?','Z'}
    end},

    -- highlight letters for jumping with f/F/t/T
    {"unblevable/quick-scope", enabled=false, config=function()
        -- Trigger a highlight in the appropriate direction when pressing these keys:
        vim.g.qs_highlight_on_keys = {'f', 'F', 't', 'T'}
    end},

    "mg979/vim-visual-multi",

    -- when a hex or other color is defined, highlight the text with its color
    {"NvChad/nvim-colorizer.lua", event = "VeryLazy", opts={
        -- write filetype to auto-attach to its buffer
        filetypes = {
            'lua',
            'R',
            'python',
            -- #RGB hex codes are too simple and may show up in julia in rare cases unnamed variables
            julia = { RGB = false, }
        }
    }},

    -- dim code that isn't currently being edited with :Twilight.
    {"folke/twilight.nvim", cmd={"Twilight", "TwilightEnable"}, opts = {dimming={alpha=0.5}, context=20}},
    
    {"ThePrimeagen/vim-be-good", cmd="VimBeGood"},

}

