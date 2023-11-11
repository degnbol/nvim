#!/usr/bin/env lua
return {
    {"hrsh7th/nvim-cmp",
    event = "InsertEnter",
    dependencies = {
        'onsails/lspkind.nvim', -- pretty pictograms
        -- putting completion sources as dependencies so they only load when cmp is loaded.
        'L3MON4D3/LuaSnip',
        'hrsh7th/cmp-nvim-lsp',
        'hrsh7th/cmp-buffer',
        'hrsh7th/cmp-path',
        'hrsh7th/cmp-nvim-lsp-signature-help',
        'hrsh7th/cmp-nvim-lua', -- neovim Lua API
        -- 'L3MON4D3/cmp-luasnip-choice', -- show choice node choices
        'tamago324/cmp-zsh', -- neovim zsh completion
        'hrsh7th/cmp-calc', -- quick math in completion
        'ray-x/cmp-treesitter', -- treesitter nodes
        'jmbuhr/otter.nvim', -- TODO: use this for code injected in markdown
        'chrisgrieser/cmp-nerdfont', -- :<search string> to get icons
        'KadoBOT/cmp-plugins',
        'uga-rosa/cmp-dictionary',
        'saadparwaiz1/cmp_luasnip',
        'honza/vim-snippets',
        'rafamadriz/friendly-snippets',
    },
    config=function()
        -- inspiration from https://vonheikemen.github.io/devlog/tools/setup-nvim-lspconfig-plus-nvim-cmp/
        local cmp = require "cmp"
        -- https://github.com/onsails/lspkind.nvim
        local lspkind = require "lspkind"

        -- menu=show completion menu. menuone=also when only one option. noselect=don't select automatically.
        vim.opt.completeopt = {"menu", "menuone", "noselect"}

        -- for tab completion
        local has_words_before = function()
            local line, col = unpack(vim.api.nvim_win_get_cursor(0))
            return col ~= 0 and vim.api.nvim_buf_get_lines(0, line - 1, line, true)[1]:sub(col, col):match("%s") == nil
        end

        -- add a version of cmp.mapping.close that always does the fallback.
        -- Modified based on https://github.com/hrsh7th/nvim-cmp/blob/main/lua/cmp/config/mapping.lua#L115
        -- NOTE: changing this was highly weird. I had different result in a test file by opening and closing it without making changes to config.
        cmp.mapping.closeFallback = function()
            return function(fallback)
                require('cmp').close()
                fallback()
            end
        end

        -- for tab support, code copied from https://github.com/hrsh7th/nvim-cmp/wiki/Example-mappings#luasnip
        mappings = {
            -- the Down and Up calls means we don't move in the list (default) but rather ignore the menu and move the cursor in the file.
            ['<up>'] = cmp.mapping.closeFallback(),
            ['<down>'] = cmp.mapping.closeFallback(),
            ['<C-b>'] = cmp.mapping.scroll_docs(-4),
            ['<C-f>'] = cmp.mapping.scroll_docs(4),
            ['<C-space>'] = cmp.mapping.complete(),
            ['<C-e>'] = cmp.mapping.abort(),
            ['<C-c>'] = cmp.mapping.abort(),
            -- Accept currently selected item. Set `select` to `false` to only confirm explicitly selected items.
            ['<CR>'] = cmp.mapping.confirm({select=false}),
            ["<Tab>"] = cmp.mapping(function(fallback)
                if vim.bo.filetype == "tsv" then
                    fallback()
                elseif cmp.visible() then
                    -- You could get the autosnippets to show up in the menu, 
                    -- but then they would autocomplete when you move through 
                    -- the completion items. Currently the only way I can think 
                    -- of to fix this is to either never auto complete the item 
                    -- when picking or writing something more intelligent to do 
                    -- this only when the next item is an (auto)snippet.
                    -- cmp.select_next_item({behavior=cmp.SelectBehavior.Select})
                    cmp.select_next_item()
                elseif has_words_before() then
                    cmp.complete()
                else
                    fallback()
                end
            end, { "i", "s" }),
            ["<S-Tab>"] = cmp.mapping(function(fallback)
                if vim.bo.filetype == "tsv" then
                    fallback()
                elseif cmp.visible() then
                    cmp.select_prev_item()
                else
                    fallback()
                end
            end, { "i", "s" }),
        }

        local types = require('cmp.types')

        cmp.setup {
            -- preselect = cmp.PreselectMode.None,
            snippet = {
                expand = function(args) require'luasnip'.lsp_expand(args.body) end,
            },
            view = {
                -- when menu is above, show best result at bottom instead of at top
                -- https://github.com/hrsh7th/nvim-cmp/wiki/Menu-Appearance
                entries = {name='custom', selection_order='near_cursor' },
            },
            mapping = cmp.mapping.preset.insert(mappings),
            sources = cmp.config.sources {
                { name = 'nvim_lsp' },
                { name = 'path', option = {trailing_slash=true} },
                { name = 'nvim_lsp_signature_help' },
                { name = 'luasnip', options = {show_autosnippets=true} },
                -- { name = 'luasnip_choice' },
                { name = 'calc' },
                { name = 'buffer', group_index=2 },
            },
            -- :h cmp-config.formatting.
            formatting = {
                -- https://github.com/onsails/lspkind.nvim
                -- https://github.com/hrsh7th/nvim-cmp/wiki/Menu-Appearance
                format = lspkind.cmp_format {
                    mode='symbol',
                    maxwidth=50,
                    ellipsis_char='…',
                    menu = {
                        buffer        = "",
                        nvim_lsp      = "", -- minimal
                        luasnip       = "", -- "", -- <> is also shown as the type, so it is redudant.
                        nvim_lua      = "",
                        latex_symbols = "",
                        nerdfont      = "󰊪",
                        calc          = "",
                        path          = "/",
                        dictionary    = "",
                        treesitter    = "",
                        zsh           = "󰞷",
                        plugins       = "",
                    }
                },
            },
            -- :h cmp-config.sorting.comparators
            sorting = {
                comparators = {
                    -- default sorting functions except for down prioritize text kind
                    function (entry1, entry2)
                        -- copy of compare.kind to simply push text entries to the bottom
                        -- https://github.com/hrsh7th/nvim-cmp/blob/8a3d2dd7641f75c1b6291311f56454adba79a196/lua/cmp/config/compare.lua#L50
                        local isText1 = entry1:get_kind() == types.lsp.CompletionItemKind.Text
                        local isText2 = entry2:get_kind() == types.lsp.CompletionItemKind.Text
                        if isText1 and not isText2 then return false end
                        if isText2 and not isText1 then return true end
                    end,
                    cmp.config.compare.offset,
                    cmp.config.compare.exact,
                    cmp.config.compare.score,
                    cmp.config.compare.recently_used,
                    cmp.config.compare.locality,
                    cmp.config.compare.kind,
                    cmp.config.compare.sort_text,
                    cmp.config.compare.length,
                    cmp.config.compare.order,
                },
            },
            experimental = {
                ghost_text = { hl_group = 'nontext' },
            }
        }

        cmp.setup.filetype({'markdown', 'tex'}, {
            sources = {
                { name = 'nvim_lsp' },
                { name = 'path', option = {trailing_slash=true} },
                { name = 'nvim_lsp_signature_help' },
                { name = 'luasnip', options = {show_autosnippets=true} },
                { name = 'calc' },
                { name = 'buffer', group_index=2 },
                { name = 'dictionary', keyword_length=3, max_item_count=10, group_index=2 },
            }
        })

        cmp.setup.filetype('lua', {
            sources = {
                { name = 'nvim_lua' },
                { name = 'nvim_lsp' },
                { name = 'plugins' },
                { name = 'path', option = {trailing_slash=true} },
                { name = 'nvim_lsp_signature_help' },
                { name = 'luasnip', options = {show_autosnippets=true} },
                { name = 'nerdfont', },
                { name = 'calc' },
                { name = 'buffer', group_index=2 },
            }
        })

        cmp.setup.filetype('zsh', {
            sources = {
                { name = 'zsh' },
                { name = 'nvim_lsp' },
                { name = 'path', option = {trailing_slash=true} },
                { name = 'nvim_lsp_signature_help' },
                { name = 'luasnip', options = {show_autosnippets=true} },
                { name = 'calc' },
                { name = 'buffer', group_index=2 },
            }
        })

        cmp.setup.filetype('julia', {
            sources = {
                { name = 'plotlyjs' },
                { name = 'nvim_lsp' },
                { name = 'path', option = {trailing_slash=true} },
                { name = 'nvim_lsp_signature_help' },
                { name = 'luasnip', options = {show_autosnippets=true} },
                -- { name = 'luasnip_choice' },
                { name = 'calc' },
                { name = 'buffer', group_index=2 },
            }
        })


        cmp.setup.cmdline(':', {
            sources = {
                { name = 'colorschemes' },
            }
        })

                -- TODO: take inspo from https://github.com/hrsh7th/nvim-cmp/wiki/Menu-Appearance
                -- and link completion menu colors to equivalent things
            end},

            -- alts: hrsh7th/vim-vsnip, SirVer/ultisnips, ...
            {'L3MON4D3/LuaSnip', lazy=true, -- load as cmp dependency
            dependencies="nvim-treesitter/nvim-treesitter", -- depend on treesitter for the ft_func
            -- the make command is optional: https://github.com/L3MON4D3/LuaSnip
            build="make install_jsregexp", },
            {'saadparwaiz1/cmp_luasnip', lazy=true,
            dependencies={'L3MON4D3/LuaSnip', "hrsh7th/nvim-cmp"}, config=function() 
                local luasnip = require "luasnip"
                -- works better to put it here than directly with luasnip, since we need to 
                -- require ft_func and would then replace config anyways
                luasnip.config.set_config {
                    -- https://youtu.be/Dn800rlPIho?t=440
                    -- don't jump back into exited snippet
                    history = true,
                    -- dynamic snippets update as you type
                    updateevents = "TextChanged,TextChangedI",
                    enable_autosnippets = true,
                    store_selection_keys = "<Tab>",
                    -- https://github.com/L3MON4D3/LuaSnip/blob/master/Examples/snippets.lua
                    -- Snippets aren't automatically removed if their text is deleted.
                    -- `delete_check_events` determines on which events (:h events) a check for
                    -- deleted snippets is performed.
                    -- This can be especially useful when `history` is enabled.
                    delete_check_events = "TextChanged",
                    -- treesitter-hl has 100, use something higher (default is 200).
                    ext_base_prio = 300,
                    -- minimal increase in priority.
                    ext_prio_increase = 1,

                    -- get filetype with treesitter instead of default from buffer filetype, to detect filetype in markdown code blocks.
                    -- get error if attempting to set from LuaSnip opts config above.
                    -- https://github.com/L3MON4D3/LuaSnip/blob/master/lua/luasnip/extras/filetype_functions.lua
                    -- "luasnippets" folder no longer auto-included. It also changes from tex->latex.
                    -- ft_func = require("luasnip.extras.filetype_functions").from_cursor_pos
                }

                vim.keymap.set({ "i", "s" }, "<C-k>", function ()
                    -- including expand means ctrl+k will autocomplete the first visible snippet in completion menu
                    if luasnip.expand_or_jumpable() then luasnip.expand_or_jump() end
                end, { silent = true })

                vim.keymap.set({ "i", "s" }, "<C-j>", function ()
                    if luasnip.jumpable(-1) then luasnip.jump(-1) end
                end, { silent = true })

                vim.keymap.set({ "i", "s" }, "<C-l>", function ()
                    if luasnip.choice_active() then
                        luasnip.change_choice(1)
                    else
                        -- fallback to moving right, useful for writing "..."
                        local keys = vim.api.nvim_replace_termcodes('<right>', true,false,true)
                        vim.api.nvim_feedkeys(keys, 'm', false)
                    end
                end)

            end },
            -- custom dicts and spell check that doesn't require spell and spelllang (f3fora/cmp-spell)
            {'uga-rosa/cmp-dictionary', lazy = true, -- load when cmp loads since dependency
            config=function()
                local cmpd = require "cmp_dictionary"
                local rtp = vim.opt.runtimepath:get()[1]

                cmpd.setup {
                    dic = {
                        -- dicts generated with ./spell.sh
                        ["*"] = {
                            rtp .. "/spell/custom.dic",
                            rtp .. "/spell/en.dic",
                        },
                        spelllang = { 
                            da = rtp .. "/spell/da.dic",
                        }
                    },
                    first_case_insensitive = true,
                    async = true, -- from slight but noticeable startup delay to instant.
                }

                function CmpDictUpdate()
                    if vim.bo.iminsert == 1 then
                        cmpd.update()
                    else
                        -- We want spelllang=en,da so we can underline bad 
                        -- spelling for both english and danish words,
                        -- but danish completion should only be shown when 
                        -- iminsert==1.
                        -- Save spelllang to temp var that is probably en,da
                        local spelllang = vim.bo.spelllang
                        vim.bo.spelllang = "en"
                        cmpd.update()
                        vim.bo.spelllang = spelllang
                    end
                end
                vim.defer_fn(CmpDictUpdate, 500)
                -- also trigger it when language is changed. See lua/keymap
                vim.api.nvim_create_autocmd("User", {
                    pattern = "ToggleDansk",
                    callback = CmpDictUpdate
                })
            end},
            -- these default snippets can be replaced with my custom snippets when I have enough
            {"honza/vim-snippets", lazy = true, -- load as cmp dependency
            dependencies={'saadparwaiz1/cmp_luasnip'}, config=function ()
                require("luasnip.loaders.from_snipmate").lazy_load {exclude={"tex", "julia", "all", "_"}}
            end},
            {"rafamadriz/friendly-snippets", lazy = true, -- load as cmp dependency
            dependencies={'saadparwaiz1/cmp_luasnip'}, config=function ()
                require("luasnip.loaders.from_vscode").lazy_load {exclude={"tex", "julia", "license", "global", "all"}}
            end},

            {
                'KadoBOT/cmp-plugins',
                lazy = true, -- loaded when cmp since it is a dependency
                ft = 'lua',
                opts = {files = { "nvim/lua/plugins/" }},
            },

        }
