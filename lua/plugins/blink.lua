local map = require "utils/keymap"
-- Using blink since it seems to be more responsive,
-- and allows typo, frecency, etc.

-- TODO:
-- Consider other sources. Mentioned on blink's website:
-- https://github.com/netmute/ctags-lsp.nvim
-- https://github.com/kristijanhusak/vim-dadbod-completion

-- TODO:
-- Show nth item.
-- Select nth item with e.g. 3<C-space>.
-- https://cmp.saghen.dev/recipes#select-nth-item-from-the-list

local source_icon = {
    buffer         = " ",
    omni           = " ",
    lsp            = " ",
    snippets       = " ",
    nvim_lua       = " ",
    nerdfonts      = "󰊪 ",
    latex_symbols  = " ",
    vimtex         = " ",
    calc           = " ",
    path           = "/ ",
    dictionary     = " ",
    spell          = " ",
    treesitter     = " ",
    zsh            = "󰞷 ",
    plugins        = " ",
    pymol_settings = " ",
    mlr            = "󰓫 ",
    pymol_select   = " ",
    plotly         = " ",
    kitty          = "󰄛 ",
    asciidoc       = " ",
    emoji          = "😃",
}

local zsh_sources = { "mlr_columns", "zsh", "lsp", "path", "snippets", "buffer" }

local only_snippets = false

local keymap = {
    preset = 'none',
    ['<C-p>'] = { 'select_prev', 'fallback' },
    ['<C-n>'] = { 'show', 'select_next', 'fallback' },
    ['<C-e>'] = { 'hide' },
    ['<C-y>'] = { 'select_and_accept' },
    ['<C-space>'] = { 'show', 'select_and_accept', 'fallback' },
    -- cancel = revert auto_insert and hide completion menu (but stay in insert mode)
    ['<C-c>'] = { 'cancel', 'fallback' },
    ['<C-b>'] = { 'show_documentation', 'scroll_documentation_up', 'fallback' },
    ['<C-f>'] = { 'show_documentation', 'scroll_documentation_down', 'fallback' },
    -- no fallbacks, since that is just inserting . and ,
    ['<C-,>'] = { 'snippet_backward', 'fallback' },
    ['<C-.>'] = { 'snippet_forward', 'fallback' },
    -- Using builtin.
    -- ['<C-s>'] = { 'show_signature' },
}
-- Same keymaps as opts.keymap. Consistency. <C-space> does the same etc.
local keymap_cmdline = vim.tbl_extend('force', keymap, {
    -- Removing fallback seems to fix mappings doing nothing.
    ['<C-p>'] = { 'select_prev', },
    ['<C-n>'] = { 'show', 'select_next', },
})


return {
    {
        'blink.compat',
        lazy = true,
        after = function()
            require("blink.compat").setup { impersonate_nvim_cmp = true }
        end,
    },
    {
        "lazydev.nvim",
        ft = "lua",
        after = function()
            require("lazydev").setup {
                library = {
                    { path = "nvim-treesitter-context", words = { "TSContext" } },
                    { path = "blink.cmp", words = { "blink.cmp" } },
                    { path = "blink-cmp-dictionary", words = { "blink%-cmp%-dictionary" } },
                },
            }
        end,
    },
    {
        "blink-cmp-dictionary",
        lazy = true,
    },
    {
        'blink.cmp',
        enabled = true,
        event = "DeferredUIEnter",
        before = function()
            -- Also jump between snippets with same mappings in normal mode.
            map.n('<C-,>', function() require "blink.cmp".snippet_backward() end, "snippet_backward")
            map.n('<C-.>', function() require "blink.cmp".snippet_forward() end, "snippet_forward")

            -- Complete only snippets, or if snippet is already active cycle choice node.
            -- Also works when completion menu is visible.
            -- Like pressing ? for choices.
            map({ "i", "s", "n" }, "<C-/>", function()
                local luasnip = require "luasnip"
                if luasnip.choice_active() then
                    luasnip.change_choice(1)
                else
                    local cmp = require 'blink.cmp'
                    -- This will also change an active popup menu listing to only show snippets.
                    -- Include LSP since LSP includes some snippets, then filter with transform_items to only get those.
                    only_snippets = true
                    vim.api.nvim_create_autocmd('User', {
                        pattern = 'BlinkCmpHide',
                        callback = function(event)
                            only_snippets = false
                        end,
                        once = true,
                    })
                    cmp.show({ providers = { 'snippets', 'lsp', } })
                end
            end)
            -- with shift to go backwards
            map({ "i", "s", "n" }, "<C-S-/>", function()
                local luasnip = require "luasnip"
                if luasnip.choice_active() then
                    luasnip.change_choice(-1)
                end
            end)

            -- Spell completion. We have the standard <C-x><C-s> and now <C-s> defaults to signature help.
            map({ "i", "s", "n" }, "<C-S-s>", function()
                local cmp = require 'blink.cmp'
                -- This will also change an active popup menu listing to only show snippets
                cmp.show({ providers = { 'dictionary' } })
            end)
        end,

        after = function()
            -- Deps with lazy=true but no trigger: load inline
            vim.cmd.packadd("blink.compat")
            require("blink.compat").setup { impersonate_nvim_cmp = true }
            vim.cmd.packadd("blink-cmp-dictionary")

            ---@module 'blink.cmp'
            ---@type blink.cmp.Config
            local opts = {
                enabled = function()
                    -- Auto completes aggresively in DAP REPL when no asked.
                    -- I tried setting all the auto-ish settings to false and it still does it.
                    return vim.bo.buftype ~= 'prompt'
                end,
                keymap = keymap,
                cmdline = { keymap = keymap_cmdline, },

                -- Default list of enabled providers defined so that you can extend it
                -- elsewhere in your config, without redefining it, due to `opts_extend`
                sources = {
                    default = { "lsp", "omni", "path", "snippets", "buffer" },
                    per_filetype = {
                        lua = { "snippets", "lazydev", "nvim_lua", "lsp", "path", "buffer", "nerdfonts", },
                        sh = zsh_sources,
                        zsh = zsh_sources,
                        ["sh.zsh"] = zsh_sources,
                        python = { "pymol_settings", "pymol_select", "lsp", "omni", "path", "snippets", "buffer" },
                        julia = { "plotly", "lsp", "omni", "path", "snippets", "buffer" },
                        kitty = { "kitty", "snippets", "buffer" },
                        -- lacks LSP, hence the custom asciidoc provider
                        asciidoc = { "asciidoc", "lsp", "omni", "snippets", "buffer", "dictionary" },
                        tex = {
                            "lsp",
                            "omni",
                            "path",
                            "snippets",
                            "buffer",
                            "dictionary",
                        },
                    },
                    -- don't lower snippet scores since we use so many custom ones
                    transform_items = function(_, items) return items end,
                    providers = {
                        lsp = {
                            -- Buffer is fallback of LSP by default.
                            -- Change this so we don't have to wait for zero LSP results before getting buffer.
                            fallbacks = {},
                            transform_items = function(_, items)
                                if only_snippets then
                                    return vim.tbl_filter(function(item)
                                        return item.kind == require('blink.cmp.types').CompletionItemKind.Snippet
                                    end, items)
                                end
                                return items
                            end
                        },
                        path = {
                            -- When a path is relevant to complete, it usually should.
                            score_offset = 4,
                        },
                        snippets = {
                            -- from default -3. Make sure path is stil higher.
                            score_offset = -1,
                        },
                        lazydev = {
                            name = "LazyDev",
                            module = "lazydev.integrations.blink",
                            -- make lazydev completions top priority (see `:h blink.cmp`)
                            score_offset = 100,
                        },
                        nvim_lua = {
                            name = "nvim_lua",
                            module = 'blink.compat.source',
                        },
                        zsh = {
                            name = "zsh",
                            module = 'blink.compat.source',
                        },
                        nerdfonts = {
                            name = "nerdfonts",
                            module = 'blink.compat.source',
                        },
                        pymol_settings = {
                            name = "pymol_settings",
                            module = "completion.pymol.blink_pymol_settings",
                            enabled = function() return vim.g.loaded_pymol end
                        },
                        pymol_select = {
                            name = "pymol_select",
                            module = "completion.pymol.blink_pymol_select",
                            enabled = function() return vim.g.loaded_pymol end
                        },
                        mlr_columns = {
                            name = "mlr",
                            module = "completion.mlr.blink_mlr",
                        },
                        plotly = {
                            name = "plotly",
                            module = "completion.plotlyjs.blink_plotlyjs",
                            enabled = function() return vim.g.loaded_plotly end
                        },
                        kitty = {
                            name = "kitty",
                            module = "kitty-conf",
                        },
                        asciidoc = {
                            name = "asciidoc",
                            module = "completion.asciidoc.blink_asciidoc",
                        },
                        -- too slow, use the dictionary plugin instead, it's amazing.
                        spell = {
                            name = "spell",
                            module = "completion.spell",
                            async = true,
                            max_items = 10,
                            min_keyword_length = 2,
                            -- timeout_ms = 500,
                        },
                        dictionary = {
                            module = 'blink-cmp-dictionary',
                            name = 'Dict',
                            score_offset = -100,
                            max_items = 10,
                            min_keyword_length = 3,
                            --- @module 'blink-cmp-dictionary'
                            --- @type blink-cmp-dictionary.Options
                            opts = {
                                dictionary_files = {
                                    vim.fn.expand("~/.config/nvim/spell/en.dic"),
                                    vim.fn.expand("~/.config/nvim/spell/custom.utf8.add"),
                                }
                            }
                        },
                        buffer = {
                            score_offset = -10,
                            -- keep case of first char
                            transform_items = function(a, items)
                                local keyword = a.get_keyword()
                                local correct, case
                                if keyword:match('^%l') then
                                    correct = '^%u%l+$'
                                    case = string.lower
                                elseif keyword:match('^%u') then
                                    correct = '^%l+$'
                                    case = string.upper
                                else
                                    return items
                                end
                                -- avoid duplicates from the corrections
                                local seen = {}
                                local out = {}
                                for _, item in ipairs(items) do
                                    local raw = item.insertText
                                    if raw:match(correct) then
                                        local text = case(raw:sub(1, 1)) .. raw:sub(2)
                                        item.insertText = text
                                        item.label = text
                                    end
                                    if not seen[item.insertText] then
                                        seen[item.insertText] = true
                                        table.insert(out, item)
                                    end
                                end
                                return out
                            end
                        },
                        -- not currently working
                        plugins = {
                            name = "plugins",
                            module = 'blink.compat.source',
                            opts = { files = { "nvim/lua/plugins/" } },
                        },
                        -- not currently working
                        treesitter = {
                            name = "treesitter",
                            module = 'blink.compat.source',
                        },
                        -- not currently working
                        vimtex = {
                            name = "vimtex",
                            module = 'blink.compat.source',
                        },
                        -- emoji = {
                        --     module = "blink-emoji",
                        --     name = "Emoji",
                        --     score_offset = 15, -- Tune by preference
                        --     opts = { insert = true }, -- Insert emoji (default) or complete its name
                        -- },
                    },
                },
                completion = {
                    list = {
                        selection = {
                            preselect = false,
                            auto_insert = function(ctx)
                                return ctx.mode ~= 'cmdline' and vim.bo.buftype ~= "prompt"
                            end,
                        },
                    },
                    accept = {
                        auto_brackets = {
                            override_brackets_for_filetypes = {
                                tex = { "{", "}" },
                            },
                            -- Synchronously use the kind of the item to determine if brackets should be added
                            kind_resolution = {
                                enabled = true,
                                blocked_filetypes = {
                                    'tex', -- incorrectly inserting brackets after e.g. label completion \cref{fig:...|}
                                },
                            },
                            -- Asynchronously use semantic token to determine if brackets should be added
                            semantic_token_resolution = {
                                enabled = true,
                                -- blocked_filetypes = { },
                                -- How long to wait for semantic tokens to return before assuming no brackets should be added
                                timeout_ms = 400,
                            },
                        },
                    },
                    menu = {
                        -- don't autoshow completion in cmdline
                        auto_show = function(ctx)
                            return ctx.mode ~= 'cmdline' and vim.bo.buftype ~= "prompt"
                        end,
                        draw = {
                            components = {
                                source_icon = {
                                    ellipsis = false,
                                    text = function(ctx)
                                        return source_icon[ctx.source_id] or " "
                                    end,
                                    highlight = 'BlinkCmpSource',
                                },
                            },
                            -- icon at end instead of before word
                            columns = {
                                { "label", gap = 1, "label_description" },
                                {
                                    "kind_icon",
                                    gap = 1,
                                    "source_icon",
                                    -- "source_name"
                                }
                            },
                            -- padding = {0,1},
                            padding = 0,
                            -- Use treesitter to highlight the label text for the given list of sources
                            -- Broken in sh for some reason
                            -- treesitter = { 'lsp' },
                        }
                    },
                    -- Ghost text clashes with auto_insert
                    -- ghost_text = { enabled = true },
                    trigger = {
                        -- Block ':' as a trigger character in typst so it doesn't reset
                        -- the completion context when typing @fig:label references.
                        show_on_blocked_trigger_characters = function()
                            if vim.bo.filetype == 'typst' then
                                return { ' ', '\n', '\t', ':' }
                            end
                            return { ' ', '\n', '\t' }
                        end,
                    },
                    keyword = {
                        -- 'prefix' will fuzzy match on the text before the cursor
                        -- 'full' will fuzzy match on the text before *and* after the cursor
                        -- example: 'foo_|_bar' will match 'foo_' for 'prefix' and 'foo__bar' for 'full'
                        -- range = 'prefix',
                        -- Regex used to get the text when fuzzy matching
                        -- regex = '[-_]\\|\\k',
                        -- After matching with regex, any characters matching this regex at the prefix will be excluded
                        -- exclude_from_prefix_regex = '[\\-]',
                    },
                    documentation = {
                        auto_show = true,
                        auto_show_delay_ms = 0,
                    },
                },
                -- Experimental signature help support.
                -- https://cmp.saghen.dev/configuration/reference.html#signature
                signature = {
                    -- Work on the builtin signature_help.
                    enabled = false,
                    trigger = {
                        -- Don't auto-open signature help. It's too distracting in e.g. julia.
                        -- Actually have to be enabled to have the float update.
                        enabled = true,
                        -- Show the signature help window when the cursor comes after a trigger character when entering insert mode.
                        -- The line above describes the intended purpose, but actually we need this on to update the float.
                        show_on_insert_on_trigger_character = true,
                    },
                },
                appearance = {
                    -- Defaults commented out.
                    -- https://cmp.saghen.dev/configuration/reference.html#completion-menu
                    kind_icons = {
                        -- Text = '󰉿',
                        -- Method = '󰊕',
                        -- Function = '󰊕',
                        -- Constructor = '󰒓',
                        --
                        -- Field = '󰜢',
                        -- Variable = '󰆦',
                        -- Property = '󰖷',
                        --
                        -- Class = '󱡠',
                        -- Interface = '󱡠',
                        -- Struct = '󱡠',
                        -- Module = '󰅩',
                        --
                        -- Unit = '󰪚',
                        -- Value = '󰦨',
                        -- Enum = '󰦨',
                        -- EnumMember = '󰦨',
                        --
                        -- Keyword = '󰻾',
                        -- Constant = '󰏿',
                        --
                        -- Snippet = '󱄽',
                        -- Color = '󰏘',
                        -- File = '󰈔',
                        -- Reference = '󰬲',
                        -- Folder = '󰉋',
                        -- Event = '󱐋',
                        Operator = '±',
                        -- TypeParameter = '󰬛',
                    },
                },
                snippets = {
                    preset = "luasnip",
                },
            }

            require('blink.cmp').setup(opts)

            -- In typst, treat ':' as part of keywords so completion stays open
            -- for @fig:label, @sec:name references. Blink's keyword chars are
            -- hardcoded (letters, digits, underscore, hyphen) in both the Rust
            -- and Lua implementations. This extends them for typst buffers.
            local fuzzy = require('blink.cmp.fuzzy')

            local orig_is_kw = fuzzy.is_keyword_character
            function fuzzy.is_keyword_character(char)
                if char == ':' and vim.bo.filetype == 'typst' then return true end
                return orig_is_kw(char)
            end

            local orig_range = fuzzy.get_keyword_range
            function fuzzy.get_keyword_range(line, col, range)
                local s, e = orig_range(line, col, range)
                if vim.bo.filetype == 'typst' then
                    -- Extend keyword range backwards through ':' segments.
                    -- e.g. for @fig:my| → extends from "my" to "fig:my"
                    while s > 0 do
                        -- line:byte(s) is Lua 1-indexed = 0-indexed s-1, the char before keyword
                        if line:byte(s) ~= 58 then break end -- 58 = ':'
                        local ps = orig_range(line, s - 1, range)
                        if ps >= s - 1 then break end -- no keyword chars before ':'
                        s = ps
                    end
                end
                return s, e
            end
        end,
    },
}
