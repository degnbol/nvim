local map = require "utils/keymap"

return {
    -- core behaviour
    "tpope/vim-repeat", -- change . to repeat last native command to last "full" command, which feels more natural.
    -- :Abolish is iabbrev which adds beginning Uppercase, whole word UPPERCASE and with {a,b} allowed for multiple combinations
    -- coerce:
    -- crs: snake_case, crm: MixedCase, crc: camelCase, cru: UPPER_CASE, cr-: dash-case, cr.: dot.case, cr<SPACE>: space case, crt: Title Case
    "tpope/vim-abolish",
    {
        "tpope/vim-unimpaired",
        config = function()
            -- next/prev file in the same dir. Not that useful so we replace it
            -- with function by treesitter (which defaults to m but we want
            -- that for math in latex)
            map.desc('n', '[f', "File prev")
            map.desc('n', ']f', "File next")
            map.desc('n', '[e', "Exchange line before")
            map.desc('n', ']e', "Exchange line after")
            map.desc('n', ']n', "Conflict")
            map.desc('n', ']n', "Conflict")
            -- TODO: move escape motions to more useful keymap. [y and ]y is used by yoink to swap yanks.
            -- map.desc('n', '[y', "C style escape motion")
            -- map.desc('n', ']y', "C style unescape motion")
            map.desc('n', '[p', "Put above, same indent")
            map.desc('n', '[P', "Put above, same indent")
            map.desc('n', ']p', "Put below, same indent")
            map.desc('n', ']P', "Put below, same indent")
            map.desc('n', '>P', "Put above, incr indent")
            map.desc('n', '>p', "Put below, incr indent")
            map.desc('n', '<P', "Put above, decr indent")
            map.desc('n', '<p', "Put below, decr indent")
            map.desc('n', '=P', "Put above, reindent")
            map.desc('n', '=p', "Put below, reindent")
        end,
    },
    {
        -- unicode shown for ga, we have go-align on that so have to use gA
        "tpope/vim-characterize",
        lazy = true,
        keys = "gA",
        init = function()
            map.n("gA", "ga", "Char info")
        end
    },
    "farmergreg/vim-lastplace", -- open file in last edited location
    -- Open a file on a line with `nvim FILEPATH:LINENUMBER`.
    -- Also supports specific column.
    "lewis6991/fileline.nvim",
    "haya14busa/vim-asterisk", -- improvements to z* and visual *. See git for uses https://github.com/haya14busa/vim-asterisk
    -- easily define custom textobjects
    -- some config in ftplugin/tex.lua
    "kana/vim-textobj-user",
    -- List of uses of the kana plugin:
    -- https://github.com/kana/vim-textobj-user/wiki
    -- provides iC and aC for multi line comments.
    -- Also provides ic which will catch multiline so not ideal but the
    -- treesitter one doesn't seem to work.
    {
        "glts/vim-textobj-comment",
        dependencies = { "kana/vim-textobj-user" },
    },
    -- nvim version of the kana plugin
    {
        "chrisgrieser/nvim-various-textobjs",
        config = function()
            require("various-textobjs").setup {
                keymaps = {
                    useDefaults = true,
                    disabledDefaults = { "gc" }, -- breaks go comment in visual
                },
            }
        end
    },
    -- increment and decrement numbers, dates, color hex, even bool
    {
        "monaqa/dial.nvim",
        lazy = true,
        keys = { "<C-a>", "<C-x>", "g<C-a>", "g<C-x>" },
        config = function()
            local augend = require("dial.augend")
            require("dial.config").augends:register_group {
                -- you can make other groups than "default" but it's not about filetype, so it's more effort than just adding a new key below.
                default = {
                    augend.integer.alias.decimal_int,
                    augend.constant.alias.bool, -- boolean value (true <-> false)
                    augend.hexcolor.new { case = "lower", },
                    augend.date.alias["%Y/%m/%d"],
                    augend.date.alias["%Y-%m-%d"],
                    augend.date.alias["%m/%d"],
                    augend.date.alias["%H:%M"],
                    -- tex
                    augend.constant.new { elements = {
                        "tiny", "scriptsize", "footnotesize", "small", "normalsize", "large", "Large", "LARGE", "huge", "Huge" },
                        word = true, cyclic = false, },
                    -- python
                    augend.constant.new { elements = { "True", "False" }, word = true, cyclic = true, },
                },
            }
        end,
        init = function()
            local function map_dial(mode, mode_name, lhs, incr_or_decr)
                map(mode, lhs, function()
                    require("dial.map").manipulate(incr_or_decr, mode_name)
                end, { desc = incr_or_decr:sub(1, 4) .. " value" })
            end
            map_dial('n', "normal", '<C-a>', "increment")
            map_dial('n', "normal", '<C-x>', "decrement")
            map_dial('n', "gnormal", 'g<C-a>', "increment")
            map_dial('n', "gnormal", 'g<C-x>', "decrement")
            map_dial('v', "visual", '<C-a>', "increment")
            map_dial('v', "visual", '<C-x>', "decrement")
            map_dial('v', "gvisual", 'g<C-a>', "increment")
            map_dial('v', "gvisual", 'g<C-x>', "decrement")
        end
    },
    {
        "ludovicchabant/vim-gutentags",
        -- NOTE: can't use ft filter with this plugin
        init = function()
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
                -- "lua",
                "markdown",
                "quarto",
                "r",
                "rust",
                "sql",
                -- "tex",
                "zsh",
                "sh",
                "json",
                "bibtex",
            }
            -- to be extra explicit I specify a comma separated list of case
            -- insensitive enabled languages
            vim.g.gutentags_ctags_extra_args = {
                "--languages=asciidoc,julia,java,lua,tex",
            }
            vim.g.gutentags_ctags_extra_args = {
                -- tried, almost worked. Was slow when worked and only inside buffer.
                -- [[--langdef=tex]],
                -- [[--langmap=tex:.tex]],
                -- [[--regex-tex=/^\\newglossaryentry\{([A-Za-z0-9 ]+)/\1/g,glossary/i]],
            }
        end,
    },

    -- supposedly faster and less buggy version of neovim builtin
    -- (:h matchparen) which highlights matching parenthesis etc.
    "monkoose/matchparen.nvim",

    -- dim code that isn't currently being edited with :Twilight.
    {
        "folke/twilight.nvim",
        lazy = true,
        cmd = { "Twilight", "TwilightEnable" },
        opts = { dimming = { alpha = 0.5 }, context = 20 },
    },

}
