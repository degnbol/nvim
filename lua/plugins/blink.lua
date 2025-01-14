-- Using blink since it seems to be more responsive,
-- and allows typo, frecency, etc.

-- TODO:
-- Consider other sources. Mentioned on blink's website:
-- https://github.com/netmute/ctags-lsp.nvim
-- https://github.com/kristijanhusak/vim-dadbod-completion

local source_icon = {
    buffer        = " ",
    omni          = " ",
    lsp           = " ",
    luasnip       = " ",
    nvim_lua      = " ",
    nerdfonts     = "󰊪 ",
    latex_symbols = " ",
    vimtex        = " ",
    calc          = " ",
    path          = "/ ",
    dictionary    = " ",
    spell         = " ",
    treesitter    = " ",
    zsh           = "󰞷 ",
    plugins       = " ",
    pymol_settings= " ",
    plotly        = " ",
    asciidoc      = " ",
    emoji         = "😃",
}

local zsh_sources = { "zsh", "lsp", "path", "snippets", "buffer" }

local only_snippets = false

return {
    {
        'saghen/blink.compat',
        -- use the latest release, via version = '*', if you also use the latest release for blink.cmp
        version = '*',
        lazy = true,
        opts = {impersonate_nvim_cmp = true,},
    },
    {
        "folke/lazydev.nvim",
        ft = "lua",
    },
    {
        "Kaiser-Yang/blink-cmp-dictionary",
        lazy = true,
        build = "brew install wordnet",
    },
    {
        'saghen/blink.cmp',

        dependencies = {
           {"L3MON4D3/LuaSnip", version = 'v2.*' },
            -- sources
            "folke/lazydev.nvim",
            "hrsh7th/cmp-nvim-lua", -- neovim Lua API
            "hrsh7th/cmp-omni", -- useful for vimscript syntax hl if `vim.opt_local.omnifunc = "syntaxcomplete#Complete"`
            -- "L3MON4D3/cmp-luasnip-choice", -- show choice node choices
            "tamago324/cmp-zsh",         -- neovim zsh completion
            -- "hrsh7th/cmp-calc",          -- quick math in completion
            "ray-x/cmp-treesitter",      -- treesitter nodes
            -- "jmbuhr/otter.nvim",         -- TODO: use this for code injected in markdown
            -- "chrisgrieser/cmp-nerdfont", -- :<search string> to get icons
            "KadoBOT/cmp-plugins",
            "Kaiser-Yang/blink-cmp-dictionary",
            -- "moyiz/blink-emoji.nvim",
            -- "micangl/cmp-vimtex",
            "davidmh/cmp-nerdfonts",
        },

        -- use a release tag to download pre-built binaries
        version = '*',
        -- AND/OR build from source, requires nightly: https://rust-lang.github.io/rustup/concepts/channels.html#working-with-nightly-rust
        -- build = 'cargo build --release',
        -- If you use nix, you can build from source using latest nightly rust with:
        -- build = 'nix run .#build-plugin',

        init = function ()
            -- We are managing snippets with luasnip instead of default currently, since we already wrote many snippets with luasnip and it allows for more complexity.
            -- However, we have also written some quick ones in VS code style that we should have available.
            require("luasnip.loaders.from_vscode").lazy_load({ paths = "~/.config/nvim/snippets" })

            -- Also jump between snippets with same mappings in normal mode.
            vim.keymap.set('n', '<C-,>', function ()
                require"blink.cmp".snippet_backward()
            end, { desc="snippet_backward" })
            vim.keymap.set('n', '<C-.>', function ()
                require"blink.cmp".snippet_forward()
            end, { desc="snippet_forward" })

            -- Complete only snippets, or if snippet is already active cycle choice node.
            -- Also works when completion menu is visible.
            -- Like pressing ? for choices.
            vim.keymap.set({ "i", "s", "n" }, "<C-/>", function()
                local luasnip = require"luasnip"
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
                    cmp.show({ providers = { 'snippets', 'lsp', }})
                end
            end)
            -- with shift to go backwards
            vim.keymap.set({ "i", "s", "n" }, "<C-S-/>", function()
                local luasnip = require"luasnip"
                if luasnip.choice_active() then
                    luasnip.change_choice(-1)
                end
            end)

            -- Spell completion.
            vim.keymap.set({ "i", "s", "n" }, "<C-s>", function()
                local cmp = require 'blink.cmp'
                -- This will also change an active popup menu listing to only show snippets
                cmp.show({ providers = { 'dictionary' } })
            end)

            -- underline active parameter in signature help rather than colour it in some pale unhelpful colour
            vim.api.nvim_set_hl(0, "BlinkCmpSignatureHelpActiveParameter", {link="Underlined"})
        end,

        ---@module 'blink.cmp'
        ---@type blink.cmp.Config
        opts = {
            -- 'default' for mappings similar to built-in completion
            -- 'super-tab' for mappings similar to vscode (tab to accept, arrow keys to navigate)
            -- 'enter' for mappings similar to 'super-tab' but with 'enter' to accept
            -- See the full "keymap" documentation for information on defining your own keymap.
            keymap = {
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
                ['<C-,>'] = { 'snippet_backward' },
                ['<C-.>'] = { 'snippet_forward' },
            },

            -- Default list of enabled providers defined so that you can extend it
            -- elsewhere in your config, without redefining it, due to `opts_extend`
            sources = {
                default = { "lsp", "omni", "path", "snippets", "buffer" },
                per_filetype = {
                    lua = { "snippets", "lazydev", "nvim_lua", "lsp", "path", "buffer", "nerdfonts", },
                    sh = zsh_sources,
                    zsh = zsh_sources,
                    ["sh.zsh"] = zsh_sources,
                    python = { "pymol_settings", "lsp", "omni", "path", "snippets", "buffer" },
                    julia = { "plotly", "lsp", "omni", "path", "snippets", "buffer" },
                    -- lacks LSP, hence the custom asciidoc provider
                    asciidoc = {"asciidoc", "lsp", "omni", "snippets", "buffer", "dictionary" },
                    tex = { "lsp", "omni", "path", "snippets", "buffer", "dictionary", },
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
                    omni = {
                        name = "omni",
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
                        enabled = function () return vim.g.loaded_pymol end
                    },
                    plotly = {
                        name = "plotly",
                        module = "completion.plotlyjs.blink_plotlyjs",
                        enabled = function () return vim.g.loaded_plotly end
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
                            get_command = {
                                'rg', -- make sure this command is available in your system
                                '--color=never',
                                '--no-line-number',
                                '--no-messages',
                                '--no-filename',
                                '--ignore-case',
                                '--',
                                '${prefix}', -- this will be replaced by the result of 'get_prefix' function
                                vim.fn.expand("~/.config/nvim/spell/en.dic"),
                                vim.fn.expand("~/.config/nvim/spell/custom.utf8.add"),
                            },
                            documentation = {
                                enable = true, -- enable documentation to show the definition of the word
                                get_command = {
                                    'wn', -- make sure this command is available in your system
                                    '${word}', -- this will be replaced by the word to search
                                    '-over'
                                }
                            }
                        }
                    },
                    buffer = {
                        -- keep case of first char
                        transform_items = function (a, items)
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
                                    local text = case(raw:sub(1,1)) .. raw:sub(2)
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
                list = { selection = { preselect = false, auto_insert = true },},
                accept = {
                    auto_brackets = {
                        override_brackets_for_filetypes = {
                            tex = {"{", "}"},
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
                    auto_show = function(ctx) return ctx.mode ~= 'cmdline' end,
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
                            { "label", "label_description" },
                            { "kind_icon", gap=1, "source_icon",
                                -- "source_name"
                            }
                        },
                        -- padding = {0,1},
                        padding = 0,
                        -- Use treesitter to highlight the label text for the given list of sources
                        treesitter = { 'lsp' },
                    }
                },
                -- Ghost text clashes with auto_insert
                -- ghost_text = { enabled = true },
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
            -- Experimental signature help support
            signature = { enabled = true },
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
        },
        opts_extend = { "sources.default" }
    },
}
