local map = require "utils/keymap"

return {
    -- core behaviour
    "tpope/vim-repeat", -- change . to repeat last native command to last "full" command, which feels more natural.
    -- :Abolish is iabbrev which adds beginning Uppercase, whole word UPPERCASE and with {a,b} allowed for multiple combinations
    -- coerce:
    -- crs: snake_case, crm: MixedCase, crc: camelCase, cru: UPPER_CASE, cr-: dash-case, cr.: dot.case, cr<SPACE>: space case, crt: Title Case
    "tpope/vim-abolish",
    -- monkoose/matchparen.nvim is supposedly faster and less buggy version of neovim default plugin matchparen,
    -- which detects matching parenthesis etc., highlights them with `MatchParen` and allows jumping to them with %.
    -- { "monkoose/matchparen.nvim", config = true, },
    -- We use the popular vim-matchup, even though the author monkoose once said it is really slow.
    -- I'm not seeing the slowness, and it has matching for quotation marks around strings.
    "andymass/vim-matchup",
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
    -- Open a file on a line with `nvim FILEPATH:LINENUMBER`.
    -- Also supports specific column.
    "lewis6991/fileline.nvim",
    -- easily define custom textobjects
    -- some config in ftplugin/tex.lua
    "kana/vim-textobj-user",
    -- List of uses of the kana plugin:
    -- https://github.com/kana/vim-textobj-user/wiki
    -- provides iC and aC for multi line comments.
    -- Also provides ic which will catch multiline so not ideal but the
    -- treesitter one doesn't seem to work.
    {
        "kana/vim-arpeggio",
        config = function ()
            vim.fn["arpeggio#map"]('i', '', 0, 'jk', '<Esc>')
        end,
    },
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
        -- Lazy-loading seems to work immediately but stops working after something else loads.
        -- keys = { "<C-a>", "<C-x>", "g<C-a>", "g<C-x>" },
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
        "linrongbin16/gentags.nvim",
        opts = {
            -- Auto cleaning up tag cache.
            gc = {
                -- Remove old tags when cache gets too big.
                trigger = "20MB",
            },
            -- Excluded filetypes.
            -- Exclude filetypes as default in order to opt-in.
            -- Here listing all that I use from `ctags --list-languages`,
            -- especially those that has LSP and/or treesitter.
            -- Using tags when there is LSP is useful for using c-] inplace of
            -- gd when I import using ROOT.
            disabled_filetypes = {
                -- "asciidoc",
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
            },
            -- ctags options
            ctags = {
                "--tag-relative=never",

                -- Recommended Options:

                -- exclude logs
                "--exclude=*.log",

                -- exclude vcs
                "--exclude=*.git",
                "--exclude=*.svg",
                "--exclude=*.hg",

                -- exclude nodejs
                "--exclude=node_modules",

                -- exclude tags/cscope
                "--exclude=*tags*",
                "--exclude=*cscope.*",

                -- exclude python
                "--exclude=*.pyc",

                -- exclude jvm class
                "--exclude=*.class",

                -- exclude VS project generated
                "--exclude=*.pdb",
                "--exclude=*.sln",
                "--exclude=*.csproj",
                "--exclude=*.csproj.user",

                -- exclude blobs
                "--exclude=*.exe",
                "--exclude=*.dll",
                "--exclude=*.mp3",
                "--exclude=*.ogg",
                "--exclude=*.flac",
                "--exclude=*.swp",
                "--exclude=*.swo",
                "--exclude=*.bmp",
                "--exclude=*.gif",
                "--exclude=*.ico",
                "--exclude=*.jpg",
                "--exclude=*.png",
                "--exclude=*.rar",
                "--exclude=*.zip",
                "--exclude=*.tar",
                "--exclude=*.tar.gz",
                "--exclude=*.tar.xz",
                "--exclude=*.tar.bz2",
                "--exclude=*.pdf",
                "--exclude=*.doc",
                "--exclude=*.docx",
                "--exclude=*.ppt",
                "--exclude=*.pptx",
            },
        },
    },
}
