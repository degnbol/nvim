#!/usr/bin/env lua
return {
    -- completion menu using builtin LSP
    {"hrsh7th/nvim-cmp", dependencies = {
        'hrsh7th/cmp-nvim-lsp', 'hrsh7th/cmp-buffer',
        'hrsh7th/cmp-path',
        'hrsh7th/cmp-nvim-lsp-signature-help',
        'hrsh7th/cmp-nvim-lua', -- neovim Lua API
        'tamago324/cmp-zsh', -- neovim zsh completion
        'onsails/lspkind.nvim', -- pretty pictograms
        'hrsh7th/cmp-calc', -- quick math in completion
    }, config=function()
        -- inspiration from https://vonheikemen.github.io/devlog/tools/setup-nvim-lspconfig-plus-nvim-cmp/
        local cmd = vim.cmd
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
        cmp.setup {
            snippet = {
                expand = function(args) require'luasnip'.lsp_expand(args.body) end,
            },
            window = {
                -- completion = cmp.config.window.bordered(),
                -- documentation = cmp.config.window.bordered(),
            },
            view = {
                -- when menu is above, show best result at bottom instead of at top
                -- https://github.com/hrsh7th/nvim-cmp/wiki/Menu-Appearance
                entries = {name = 'custom', selection_order = 'near_cursor' }
            },
            mapping = cmp.mapping.preset.insert {
                -- the Down and Up calls means we don't move in the list (default) but rather ignore the menu and move the cursor in the file.
                ['<up>'] = cmp.mapping.closeFallback(),
                ['<down>'] = cmp.mapping.closeFallback(),
                ['<C-b>'] = cmp.mapping.scroll_docs(-4),
                ['<C-f>'] = cmp.mapping.scroll_docs(4),
                ['<c-space>'] = cmp.mapping.complete(),
                ['<C-e>'] = cmp.mapping.abort(),
                -- Accept currently selected item. Set `select` to `false` to only confirm explicitly selected items.
                ['<CR>'] = cmp.mapping.confirm({select=false}),
                ["<Tab>"] = cmp.mapping(function(fallback)
                    if cmp.visible() then
                        cmp.select_next_item()
                    elseif require'luasnip'.expand_or_jumpable() then
                        require'luasnip'.expand_or_jump()
                    elseif has_words_before() then
                        cmp.complete()
                    else
                        fallback()
                    end
                end, { "i", "s" }),
                ["<S-Tab>"] = cmp.mapping(function(fallback)
                    if cmp.visible() then
                        cmp.select_prev_item()
                    elseif require'luasnip'.jumpable(-1) then
                        require'luasnip'.jump(-1)
                    else
                        fallback()
                    end
                end, { "i", "s" }),
            },
            sources = cmp.config.sources {
                { name = 'nvim_lsp' },
                { name = 'path', option = {trailing_slash=true} },
                { name = 'nvim_lsp_signature_help' },
                { name = 'luasnip' },
                { name = 'calc' },
                { name = 'buffer' },
            },
            formatting = {
                format = lspkind.cmp_format {
                    mode = 'symbol',
                    maxwidth = 50,
                }
            },
        }

        cmp.setup.filetype({'markdown', 'tex'}, {
            sources = {
                { name = 'nvim_lsp' },
                { name = 'path', option = {trailing_slash=true} },
                { name = 'nvim_lsp_signature_help' },
                { name = 'luasnip' },
                { name = 'calc' },
                { name = 'buffer', group_index=2 },
                { name = 'dictionary', keyword_length=3, max_item_count=10, group_index=2 },
            }
        })

        cmp.setup.filetype('lua', {
            sources = {
                { name = 'nvim_lua' },
                { name = 'nvim_lsp' },
                { name = 'path', option = {trailing_slash=true} },
                { name = 'nvim_lsp_signature_help' },
                { name = 'luasnip' },
                { name = 'calc' },
                { name = 'buffer', group_index=2  },
            }
        })

        cmp.setup.filetype('zsh', {
            sources = {
                { name = 'zsh' },
                { name = 'nvim_lsp' },
                { name = 'path', option = {trailing_slash=true} },
                { name = 'nvim_lsp_signature_help' },
                { name = 'luasnip' },
                { name = 'calc' },
                { name = 'buffer', group_index=2  },
            }
        })
    end},
    {'saadparwaiz1/cmp_luasnip', dependencies={'L3MON4D3/LuaSnip', "hrsh7th/nvim-cmp"}, config=function() 
        local luasnip = require "luasnip"
        -- https://youtu.be/Dn800rlPIho?t=440
        luasnip.config.set_config {
            -- don't jump back into exited snippet
            history = false,
            -- dynamic snippets update as you type
            updateevents = "TextChanged,TextChangedI",
            enable_autosnippets = true,
        }
    end }, -- alts: hrsh7th/vim-vsnip, SirVer/ultisnips, ...
    -- custom dicts and spell check that doesn't require spell and spelllang (f3fora/cmp-spell)
    {'uga-rosa/cmp-dictionary', config=function()
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

        -- We want spelllang=en,da so we can underline bad spelling in both, 
        -- but toggle completion from danish only when iminsert=1.
        function CmpDictUpdate()
            if vim.bo.iminsert == 1 then
                cmpd.update()
            else
                local spelllang = vim.bo.spelllang
                vim.bo.spelllang = "en"
                cmpd.update()
                vim.bo.spelllang = spelllang
            end
        end

        -- was needed with 0 ms even without async=true.
        -- 1000 ms should be a second but it seems it is called instantly after everything else?
        -- Anyways, it works.
        vim.defer_fn(CmpDictUpdate, 1000)

    end},
    {"rafamadriz/friendly-snippets", dependencies={'saadparwaiz1/cmp_luasnip'}, config=function ()
        -- load friendly-snippets with luasnip
        require("luasnip.loaders.from_vscode").lazy_load()
    end}
}
