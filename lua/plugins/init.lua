local util = require "utils/init"
local map = require "utils/keymap"

return {
    -- core behaviour
    "tpope/vim-repeat", -- change . to repeat last native command to last "full" command, which feels more natural.
    -- Y should yank to end of line which is consistent with other uppercase use, rather than yank whole line like yy which is for ancient vi compatibility.
    -- "tpope/vim-sensible",
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
        keys = "gA",
        init = function()
            vim.keymap.set("n", "gA", "ga", { desc = "Char info" })
        end
    },
    "farmergreg/vim-lastplace", -- open file in last edited location
    "haya14busa/vim-asterisk",  -- improvements to z* and visual *. See git for uses https://github.com/haya14busa/vim-asterisk
    -- easily define custom textobjects
    -- some config in ftplugin/tex.lua
    "kana/vim-textobj-user",
    -- List of uses of the kana plugin:
    -- https://github.com/kana/vim-textobj-user/wiki
    -- provides iC and aC for multi line comments.
    -- Also provides ic which will catch multiline so not ideal but the
    -- treesitter one doesn't seem to work.
    { "glts/vim-textobj-comment", dependencies = { "kana/vim-textobj-user" } },
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
            -- custom gx function that opens github repos given the short
            -- version written in these config files.
            vim.keymap.set("n", "gx", function()
                -- go to github for plugin easily.
                -- First check if we are editing a file read by lazy.nvim
                if vim.api.nvim_buf_get_name(0):find('nvim/lua/plugins/') then
                    -- get first string on the line, assumes we don't list multiple plugins on one line.
                    local line = vim.api.nvim_get_current_line()
                    local repo = line:match([["([%w%p]+/[%w%p]+)"]])
                    repo = repo or line:match([['([%w%p]+/[%w%p]+)']])
                    if repo then
                        return util.open("https://github.com/" .. repo)
                    end
                end

                -- for latex packages...
                if vim.bo.filetype == "tex" then
                    local line = vim.api.nvim_get_current_line()
                    local pac = line:match("\\usepackage.*{([%w_-]+)}")
                    if pac ~= nil then
                        local ctan = "https://ctan.org/pkg/" .. pac .. "?lang=en"
                        return util.open(ctan)
                    end
                end

                -- visually select URL
                require("various-textobjs").url()
                -- plugin only switches to visual mode when textobj found
                local foundURL = vim.fn.mode():find("v")
                -- retrieve URL with the z-register as intermediary
                vim.cmd.normal { '"zy', bang = true }
                local url = vim.fn.getreg("z")
                util.open(url)
            end, { desc = "Smart URL opener" })
        end
    },
    -- increment and decrement numbers, dates, color hex, even bool
    {
        "monaqa/dial.nvim",
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

            vim.api.nvim_set_keymap("n", "<C-a>", require("dial.map").inc_normal(), { desc = "Value incr" })
            vim.api.nvim_set_keymap("n", "<C-x>", require("dial.map").dec_normal(), { desc = "Value decr" })
            vim.api.nvim_set_keymap("v", "<C-a>", require("dial.map").inc_visual(), { desc = "Value incr" })
            vim.api.nvim_set_keymap("v", "<C-x>", require("dial.map").dec_visual(), { desc = "Value decr" })
            vim.api.nvim_set_keymap("v", "g<C-a>", require("dial.map").inc_gvisual(), { desc = "Value incr" })
            vim.api.nvim_set_keymap("v", "g<C-x>", require("dial.map").dec_gvisual(), { desc = "Value decr" })
        end
    },

    "lervag/file-line", -- open a file on a line with vi filepath:linenumber

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

    -- supposedly faster and less buggy version of neovim builtin (:h )matchparen which highlights matching parenthesis etc.
    "monkoose/matchparen.nvim",

    {
        'ggandor/leap.nvim',
        enabled = false,
        config = function()
            -- make s and S "unsafe", i.e. available immediately as a command
            -- add ' and ` as safe since it would be unlikely that I would want to jump to a mark right after a leap
            -- add [] as safe since it would be unlikely that I would want to jump with those after a leap
            require 'leap'.opts.safe_labels = { 'f', 'n', 'u', 't', '/', '`', "'", '[', ']', 'F', 'N', 'L', 'H', 'M', 'U',
                'G', 'T', '?', 'Z' }

            vim.keymap.set({ 'n', 'x', 'o' }, "\\", '<Plug>(leap-forward-to)', { desc = "forward to (leap)" })
            vim.keymap.set({ 'n', 'x', 'o' }, "|", '<Plug>(leap-backward-to)', { desc = "backward to (leap)" })
        end
    },

    -- when a hex or other color is defined, highlight the text with its color
    -- trying out mini hipatterns instead.
    {
        "NvChad/nvim-colorizer.lua",
        enabled = false,
        event = "VeryLazy",
        opts = {
            -- write filetype to auto-attach to its buffer
            filetypes = {
                'lua',
                'R',
                'python',
                -- #RGB hex codes are too simple and may show up in julia in rare cases unnamed variables
                julia = { RGB = false, }
            }
        }
    },

    -- dim code that isn't currently being edited with :Twilight.
    { "folke/twilight.nvim",      cmd = { "Twilight", "TwilightEnable" },    opts = { dimming = { alpha = 0.5 }, context = 20 } },

    { "ThePrimeagen/vim-be-good", cmd = "VimBeGood" },

}
