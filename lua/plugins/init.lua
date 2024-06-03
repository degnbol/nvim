#!/usr/bin/env lua
return {
    -- core behaviour
    "tpope/vim-repeat", -- change . to repeat last native command to last "full" command, which feels more natural.
    -- Y should yank to end of line which is consistent with other uppercase use, rather than yank whole line like yy which is for ancient vi compatibility.
    -- "tpope/vim-sensible",
    -- :Abolish is iabbrev which adds beginning Uppercase, whole word UPPERCASE and with {a,b} allowed for multiple combinations
    -- crs: snake_case, crm: MixedCase, crc: camelCase, cru: UPPER_CASE, cr-: dash-case, cr.: dot.case, cr<SPACE>: space case, crt: Title Case
    "tpope/vim-abolish",
    {
        "tpope/vim-unimpaired",
        config = function ()
            require "utils/keymap"
            set_keymap_desc('n', 'yo', "Option toggle")
            -- duplicate mappings
            set_keymap_desc('n', '=s', "Setting toggle | Substitute+reindent")
            set_keymap_desc('n', '<s', "Setting enable")
            set_keymap_desc('n', '>s', "Setting disable")
            -- argument list, i.e. the original list of files given to the nvim cmd.
            -- Seems more intuitive to just use ]b for buffers.
            -- We replace it with moving to next argument instead.
            set_keymap_desc('n', '[a', ":previous")
            set_keymap_desc('n', ']a', ":next")
            set_keymap_desc('n', '[A', ":first")
            set_keymap_desc('n', ']A', ":last")
            -- next/prev file in the same dir. Not that useful so we replace it 
            -- with function by treesitter (which defaults to m but we want 
            -- that for math in latex)
            set_keymap_desc('n', '[f', "File prev")
            set_keymap_desc('n', ']f', "File next")
            set_keymap_desc('n', '[e', "Exchange line before")
            set_keymap_desc('n', ']e', "Exchange line after")
            set_keymap_desc('n', ']n', "Conflict")
            set_keymap_desc('n', ']n', "Conflict")
            set_keymap_desc('n', '[y', "C style escape motion")
            set_keymap_desc('n', ']y', "C style unescape motion")
            set_keymap_desc('n', '[ ', "Add empty before")
            set_keymap_desc('n', '] ', "Add empty after")
            set_keymap_desc('n', '[p', "Put above, same indent")
            set_keymap_desc('n', '[P', "Put above, same indent")
            set_keymap_desc('n', ']p', "Put below, same indent")
            set_keymap_desc('n', ']P', "Put below, same indent")
            set_keymap_desc('n', '>P', "Put above, incr indent")
            set_keymap_desc('n', '>p', "Put below, incr indent")
            set_keymap_desc('n', '<P', "Put above, decr indent")
            set_keymap_desc('n', '<p', "Put below, decr indent")
            set_keymap_desc('n', '=P', "Put above, reindent")
            set_keymap_desc('n', '=p', "Put below, reindent")

            -- some extra functions that I felt were missing
            -- To conceal or not, set by changing the conceallevel between 0 
            -- and a non-zero value.
            local toggle_conceal = function ()
                if vim.opt.conceallevel:get() > 0 then
                    nzConcealLvl = vim.opt.conceallevel:get()
                    vim.opt.conceallevel = 0
                else
                    vim.opt.conceallevel = nzConcealLvl or 1
                end
            end
            local toggle_colcealcursor = function ()
                if vim.opt.concealcursor:get():match('n') then
                    vim.cmd 'setlocal concealcursor-=n' -- lua version not simple
                else
                    vim.opt.concealcursor:append('n')
                end
            end
            local enable_autoformat = function ()
                vim.opt.formatoptions:append('a')
                -- also remove sidescroll offset since there should be enough space on the screen
                -- Keep record of the original value
                sidescrolloff = vim.opt.sidescrolloff:get()
                vim.opt.sidescrolloff = 0
                notify = "fo+=a | sidescrolloff=0"
                print(notify)
                return notify
            end
            local disable_autoformat = function ()
                vim.opt.formatoptions:remove('a')
                -- reset sidescrolloff.
                -- if the script local var 'sidescrolloff' hasn't been defined 
                -- in a call to enable_autoformat, we set it to what is 
                -- currently the default in lua/options.lua
                vim.opt.sidescrolloff = sidescrolloff or 12
                notify = "fo-=a | sidescrolloff=" .. vim.opt.sidescrolloff:get()
                print(notify)
                return notify
            end
            local enable_wrap = function ()
                vim.opt.wrap = true
                -- also disable autoformat when wrapping
                -- Keep record of the original value
                autoformat = vim.opt.formatoptions:get()['a']
                notify = disable_autoformat() .. " | wrap"
                print(notify)
                return notify
            end
            local disable_wrap = function ()
                vim.opt.wrap = false
                -- reset autoformat.
                -- if the script local var 'autoformat' hasn't been defined 
                -- in a call to enable_wrap, we default to false
                if autoformat then
                    notify = enable_autoformat() .. " | nowrap"
                else
                    notify = "nowrap"
                end
                print(notify)
                return notify
            end
            local toggle_autoformat = function ()
                if vim.opt.formatoptions:get()['a'] then
                    disable_autoformat()
                else
                    enable_autoformat()
                end
            end
            local function toggle_wrap()
                if vim.opt.wrap:get() then
                    disable_wrap()
                else
                    enable_wrap()
                end
            end

            -- default for yoc is another binding for cursorline which is a lot 
            -- less useful than conceal
            vim.keymap.set('n', 'yoc', toggle_conceal, { desc="conceal" })
            vim.keymap.set('n', '=sc', toggle_conceal, { desc="conceal" })
            vim.keymap.set('n', 'yoC', toggle_colcealcursor, { desc="concealcursor" })
            vim.keymap.set('n', '=sC', toggle_colcealcursor, { desc="concealcursor" })
            vim.keymap.set('n', 'yoa', toggle_autoformat, { desc="formatoptions auto" })
            vim.keymap.set('n', '=sa', toggle_autoformat, { desc="formatoptions auto" })
            vim.keymap.set('n', '<sa', enable_autoformat, { desc="formatoptions auto" })
            vim.keymap.set('n', '>sa', disable_autoformat, { desc="formatoptions auto" })
            vim.keymap.set('n', 'yow', toggle_wrap, { desc="wrap" })
            vim.keymap.set('n', '=sw', toggle_wrap, { desc="wrap" })
            vim.keymap.set('n', '<sw', enable_wrap, { desc="wrap" })
            vim.keymap.set('n', '>sw', disable_wrap, { desc="wrap" })
        end,
    },
    {
        -- unicode shown for ga, we have go-align on that so have to use gA
        "tpope/vim-characterize",
        keys="gA",
        init = function ()
            vim.keymap.set("n", "gA", "ga", { desc="Char info" })
        end
    },

    -- multi cursor
    {
        "mg979/vim-visual-multi",
        config = function ()
            -- use alt instead of ctrl since ctrl arrows moves mac os windows.
            -- Other keybindings (e.g. [] and y) cannot be set by multicursor since they are set elsewhere, hence startup warnings.
            -- Ignore for now since we can live without them.
            vim.keymap.set('n', '<A-Down>', "<C-Down>", { remap=true, desc="Multi cursor down" })
            vim.keymap.set('n', '<A-Up>', "<C-Up>", { remap=true, desc="Multi cursor up" })
        end,
    },
    {
        "smoka7/multicursors.nvim",
        event = "VeryLazy",
        dependencies = {
            'nvim-treesitter/nvim-treesitter',
            'smoka7/hydra.nvim',
        },
        opts = {
            generate_hints = {
                normal = true,
                insert = true,
                extend = true,
            },
        },
        cmd = { 'MCstart', 'MCvisual', 'MCclear', 'MCpattern', 'MCvisualPattern', 'MCunderCursor' },
        keys = {
            {
                mode = { 'v', 'n' },
                '<Leader>m',
                '<Cmd>MCstart<CR>',
                desc = 'Create a selection for selected text or word under the cursor',
            },
            {
                mode = { 'v', 'n' },
                '<Leader>M',
                '<Cmd>MCunderCursor<CR>',
                desc = 'Select the char under the cursor and start listening for the actions.',
            },
        },
    },
    "farmergreg/vim-lastplace", -- open file in last edited location
    "haya14busa/vim-asterisk", -- improvements to z* and visual *. See git for uses https://github.com/haya14busa/vim-asterisk
    -- "gioele/vim-autoswap",
    -- easily define custom textobjects
    -- some config in ftplugin/tex.lua
    "kana/vim-textobj-user", 
    -- List of uses of the kana plugin:
    -- https://github.com/kana/vim-textobj-user/wiki
    -- provides iC and aC for multi line comments.
    -- Also provides ic which will catch multiline so not ideal but the 
    -- treesitter one doesn't seem to work.
    {"glts/vim-textobj-comment", dependencies={"kana/vim-textobj-user"}},
    -- nvim version of the kana plugin
    {
        "chrisgrieser/nvim-various-textobjs",
        -- opts = { useDefaultKeymaps = true },
        config = function ()
            require("various-textobjs").setup {
                useDefaultKeymaps = true,
                disabledKeymaps = {"gc"}, -- breaks go comment in visual
            }
            -- custom gx function that opens github repos given the short 
            -- version written in these config files.
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
            end, { desc = "Smart URL opener" })
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

        vim.api.nvim_set_keymap("n", "<C-a>",  require("dial.map").inc_normal(),  {desc="Value incr"})
        vim.api.nvim_set_keymap("n", "<C-x>",  require("dial.map").dec_normal(),  {desc="Value decr"})
        vim.api.nvim_set_keymap("v", "<C-a>",  require("dial.map").inc_visual(),  {desc="Value incr"})
        vim.api.nvim_set_keymap("v", "<C-x>",  require("dial.map").dec_visual(),  {desc="Value decr"})
        vim.api.nvim_set_keymap("v", "g<C-a>", require("dial.map").inc_gvisual(), {desc="Value incr"})
        vim.api.nvim_set_keymap("v", "g<C-x>", require("dial.map").dec_gvisual(), {desc="Value decr"})
    end},
    
    "lervag/file-line", -- open a file on a line with vi filepath:linenumber

    {
        "ludovicchabant/vim-gutentags",
        -- NOTE: can't use ft filter with this plugin
        init = function ()
            -- hide tag files where logs and other caches for nvim are stored
            vim.g.gutentags_cache_dir = "~/.local/state/nvim/tags/"
            -- exclude filetypes as default in order to opt-in.
            -- I added all that I use from `ctags --list-languages`, especially 
            -- those that has LSP and/or treesitter.
            -- Using tags when there is LSP is useful for using c-] inplace of gd when I import using ROOT
            vim.g.gutentags_exclude_filetypes = {
                "python",
                "matlab",
                -- "julia",
                -- "java",
                "vim",
                "lua",
                "markdown",
                "quarto",
                "r",
                "rust",
                "sql",
                "tex",
                "zsh",
                "sh",
                "json",
                "bibtex",
            }
            -- to be extra explicit I specify a comma separated list of case 
            -- insensitive enabled languages
            vim.g.gutentags_ctags_extra_args = {
                "--languages=asciidoc,julia,java",
            }
        end,
    },

    -- supposedly faster and less buggy version of neovim builtin (:h )matchparen which highlights matching parenthesis etc.
    "monkoose/matchparen.nvim",

    {
        'ggandor/leap.nvim',
        enabled = false,
        config=function()
            -- make s and S "unsafe", i.e. available immediately as a command
            -- add ' and ` as safe since it would be unlikely that I would want to jump to a mark right after a leap
            -- add [] as safe since it would be unlikely that I would want to jump with those after a leap
            require'leap'.opts.safe_labels = {'f','n','u','t','/','`',"'",'[',']','F','N','L','H','M','U','G','T','?','Z'}

            vim.keymap.set({'n', 'x', 'o'}, "\\", '<Plug>(leap-forward-to)', { desc="forward to (leap)" })
            vim.keymap.set({'n', 'x', 'o'}, "|", '<Plug>(leap-backward-to)', { desc="backward to (leap)" })
        end
    },

    -- when a hex or other color is defined, highlight the text with its color
    -- trying out mini hipatterns instead.
    {"NvChad/nvim-colorizer.lua", enabled = false, event = "VeryLazy", opts={
        -- write filetype to auto-attach to its buffer
        filetypes = {
            'lua', 'R', 'python',
            -- #RGB hex codes are too simple and may show up in julia in rare cases unnamed variables
            julia = { RGB = false, }
        }
    }},

    -- dim code that isn't currently being edited with :Twilight.
    {"folke/twilight.nvim", cmd={"Twilight", "TwilightEnable"}, opts = {dimming={alpha=0.5}, context=20}},

    {"ThePrimeagen/vim-be-good", cmd="VimBeGood"},

}

