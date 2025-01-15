#!/usr/bin/env lua
-- works if blink is removed with lazy clean
local using_blink, _ = pcall(require, "blink.cmp")

local cmp_deps
if using_blink then
    cmp_deps = {
        'L3MON4D3/LuaSnip',
        'saadparwaiz1/cmp_luasnip',
    }
else
    cmp_deps = {
        -- 'hrsh7th/cmp-nvim-lsp',
        { "iguanacucumber/mag-nvim-lsp", name = "cmp-nvim-lsp", opts = {} },
        -- 'hrsh7th/cmp-nvim-lua', -- neovim Lua API
        { "iguanacucumber/mag-nvim-lua", name = "cmp-nvim-lua" },
        -- 'degnbol/cmp-buffer',
        { "iguanacucumber/mag-buffer", name = "cmp-buffer" },
        { "iguanacucumber/mag-cmdline", name = "cmp-cmdline" },
        'onsails/lspkind.nvim', -- pretty pictograms
        -- putting completion sources as dependencies so they only load when cmp is loaded.
        'L3MON4D3/LuaSnip',
        -- 'hrsh7th/cmp-path',
        "https://codeberg.org/FelipeLema/cmp-async-path",
        'hrsh7th/cmp-nvim-lsp-signature-help',
        'hrsh7th/cmp-omni',
        -- 'L3MON4D3/cmp-luasnip-choice', -- show choice node choices
        'tamago324/cmp-zsh',         -- neovim zsh completion
        'hrsh7th/cmp-calc',          -- quick math in completion
        'ray-x/cmp-treesitter',      -- treesitter nodes
        'jmbuhr/otter.nvim',         -- TODO: use this for code injected in markdown
        'chrisgrieser/cmp-nerdfont', -- :<search string> to get icons
        'KadoBOT/cmp-plugins',
        'uga-rosa/cmp-dictionary',
        'saadparwaiz1/cmp_luasnip',
        'honza/vim-snippets',
        'rafamadriz/friendly-snippets',
    }
end

return {
    {
        -- hack solution. I want buffer completed words with either capitalization.
        -- It wasn't easy to get cmp to not change capitalization for buffer
        -- completions, so I made this fork of the repo where regular words are
        -- indexed with either capitalization.
        "degnbol/cmp-buffer",
        lazy = true, -- loaded as dependency
        enabled = false,
        -- enabled = not using_blink,
        branch = "patch-1",
    },
    {
        'KadoBOT/cmp-plugins',
        -- enabled = cmp_enabled, -- use compat layer with blink
        lazy = true, -- loaded as dependency
        ft = 'lua',
        opts = { files = { "nvim/lua/plugins/" } },
    },
    -- custom dicts and spell check that doesn't require spell and spelllang (f3fora/cmp-spell)
    {
        'uga-rosa/cmp-dictionary',
        lazy = true, -- loaded as dependency
        enabled = not using_blink,
        -- lazy = true, -- Doesn't work to lazy load.
        config = function()
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
                if vim.bo.iminsert > 0 then
                    -- set dictionary completions with both English and Danish words.
                    cmpd.update()
                else
                    -- We always want spelllang=en,da so we can underline bad
                    -- spelling for both english and danish words,
                    -- but danish completion should only be shown when
                    -- We are writing in Danish, which is indicated by activating iminsert.
                    -- So we save spelllang to temp var (value is probably "en,da")
                    local spelllang = vim.bo.spelllang
                    -- ... then set dictionary completions only in English
                    vim.bo.spelllang = "en"
                    cmpd.update()
                    -- ... but set spelllang back to what it was (probably "en,da")
                    vim.bo.spelllang = spelllang
                end
            end

            vim.defer_fn(CmpDictUpdate, 500)
            -- also trigger it when language is changed. See lua/keymap
            vim.api.nvim_create_autocmd("User", {
                pattern = "ToggleDansk",
                callback = CmpDictUpdate
            })
        end
    },
    {
        -- "hrsh7th/nvim-cmp",
        "iguanacucumber/magazine.nvim",
        name = "nvim-cmp", -- Otherwise highlighting gets messed up
        -- don't disable here so can have separate config for blink.lua
        -- enabled = not using_blink,
        -- event = "InsertEnter" NO, doesn't work, e.g. for query loading luasnip
        dependencies = cmp_deps,
        -- inspiration from https://vonheikemen.github.io/devlog/tools/setup-nvim-lspconfig-plus-nvim-cmp/
        config = function()
            local cmp = require "cmp"

            if using_blink then
                cmp.setup {
                    -- will still complete autosnippets
                    enabled = false,
                    -- preselect = cmp.PreselectMode.None,
                    snippet = {
                        expand = function(args) require 'luasnip'.lsp_expand(args.body) end,
                    },
                    mapping = nil,
                    sources = cmp.config.sources {
                        { name = 'luasnip', options = { show_autosnippets = true } },
                    },
                    formatting = nil,
                    sorting = nil,
                }
                return
            end

            -- menu=show completion menu. menuone=also when only one option. noselect=don't select automatically.
            vim.opt.completeopt = { "menu", "menuone", "noselect" }

            -- for tab completion
            local has_words_before = function()
                local line, col = unpack(vim.api.nvim_win_get_cursor(0))
                return col ~= 0 and
                    vim.api.nvim_buf_get_lines(0, line - 1, line, true)[1]:sub(col, col):match("%s") == nil
            end

            -- add a version of cmp.mapping.close that always does the fallback.
            -- Modified based on https://github.com/hrsh7th/nvim-cmp/blob/main/lua/cmp/config/mapping.lua#L115
            -- NOTE: changing this was highly weird. I had different result in a test file by opening and closing it without making changes to config.
            cmp.mapping.closeFallback = function()
                return function(fallback)
                    require'cmp'.close()
                    fallback()
                end
            end

            -- for tab support, code copied from https://github.com/hrsh7th/nvim-cmp/wiki/Example-mappings#luasnip
            local mappings = {
                -- the Down and Up calls means we don't move in the list (default) but rather ignore the menu and move the cursor in the file.
                ['<up>'] = cmp.mapping.closeFallback(),
                ['<down>'] = cmp.mapping.closeFallback(),
                ['<C-b>'] = cmp.mapping.scroll_docs(-4),
                ['<C-f>'] = cmp.mapping.scroll_docs(4),
                -- complete with c-space instead of ctrl+y or tab. Prefer ctrl-n and -p and ctrl-space.
                ['<C-space>'] = cmp.mapping(function(fallback)
                    if cmp.visible() then
                        cmp.confirm({ select = true })
                    else
                        cmp.complete()
                    end
                end, { "i", "s" }),
                ['<C-e>'] = cmp.mapping.abort(),
                ['<C-c>'] = cmp.mapping.abort(),
                ['<C-n>'] = cmp.mapping(function(fallback)
                    if vim.bo.filetype == "tsv" then
                        fallback()
                    elseif cmp.visible() then
                        -- Hack to avoid luasnip autosnip trigger when browsing through popup menu.
                        -- luasnip is the only autocmd set for the frequent InsertCharPre event,
                        -- where a flag is set to indicate to a later TextChangedI, TextChangedP autocmd that text has been inserted.
                        -- See Luasnip_just_inserted in https://github.com/L3MON4D3/LuaSnip/blob/master/lua/luasnip/config.lua
                        -- InsertCharPre is called after the insertion by cmp.select_next_item.
                        -- We change the eventignore option temporarily to ignore this event.
                        vim.opt.eventignore:append("InsertCharPre")
                        cmp.select_next_item() -- {behavior=cmp.SelectBehavior.Select}
                        -- Resetting has to be scheduled since there is some time between the above call and the trigger.
                        vim.schedule(function()
                            vim.opt.eventignore:remove("InsertCharPre")
                        end)
                    elseif has_words_before() then
                        cmp.complete()
                    else
                        fallback()
                    end
                end, { "i", "s" }),
                ["<C-p>"] = cmp.mapping(function(fallback)
                    if cmp.visible() then
                        -- see above
                        vim.opt.eventignore:append("InsertCharPre")
                        cmp.select_prev_item() -- {behavior=cmp.SelectBehavior.Select}
                        vim.schedule(function()
                            vim.opt.eventignore:remove("InsertCharPre")
                        end)
                    else
                        fallback()
                    end
                end, { "i", "s" }),
                -- I may change tab in future, since I'm using ctrl-N, ctrl-P.
                -- It currently is different in that it expands autocmds if they happen to match the popup selection.
                ["<Tab>"] = cmp.mapping(function(fallback)
                    if vim.bo.filetype == "tsv" then
                        fallback()
                    elseif cmp.visible() then
                        cmp.select_next_item() -- {behavior=cmp.SelectBehavior.Select}
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
                        cmp.select_prev_item() -- {behavior=cmp.SelectBehavior.Select}
                    else
                        fallback()
                    end
                end, { "i", "s" }),
            }

            local types = require 'cmp.types'

            cmp.setup {
                -- preselect = cmp.PreselectMode.None,
                snippet = {
                    expand = function(args) require 'luasnip'.lsp_expand(args.body) end,
                },
                view = {
                    -- when menu is above, show best result at bottom instead of at top
                    -- https://github.com/hrsh7th/nvim-cmp/wiki/Menu-Appearance
                    entries = { name = 'custom', selection_order = 'near_cursor' },
                },
                mapping = cmp.mapping.preset.insert(mappings),
                sources = cmp.config.sources {
                    { name = 'nvim_lsp' },
                    { name = 'path', option = { trailing_slash = true } },
                    { name = 'nvim_lsp_signature_help' },
                    { name = 'luasnip', options = { show_autosnippets = true } },
                    -- { name = 'luasnip_choice' },
                    { name = 'calc' },
                    { name = 'omni', group_index = 2 },
                    { name = 'buffer', group_index = 2 },
                },
                -- :h cmp-config.formatting.
                formatting = {
                    -- https://github.com/onsails/lspkind.nvim
                    -- https://github.com/hrsh7th/nvim-cmp/wiki/Menu-Appearance
                    -- https://github.com/onsails/lspkind.nvim
                    format = require "lspkind".cmp_format {
                        mode = 'symbol',
                        maxwidth = 50,
                        ellipsis_char = '…',
                        menu = {
                            buffer        = "",
                            omni          = "", -- most likely set to syntax keyword completion
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
                        function(entry1, entry2)
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
                -- it's too distracting when inside a line.
                -- I would like it to only be active at end of empty line.
                -- experimental = {
                --     ghost_text = { hl_group = 'nontext' },
                -- }
            }

            cmp.setup.filetype({ 'markdown', 'tex', 'text', 'asciidoc' }, {
                sources = {
                    { name = 'nvim_lsp' },
                    { name = 'path', option = { trailing_slash = true } },
                    { name = 'nvim_lsp_signature_help' },
                    { name = 'luasnip', options = { show_autosnippets = true } },
                    { name = 'calc' },
                    { name = 'buffer', group_index = 2 },
                    { name = 'dictionary', keyword_length = 3, max_item_count = 10, group_index = 2 },
                }
            })

            cmp.setup.filetype('lua', {
                sources = {
                    { name = 'nvim_lua' },
                    { name = 'nvim_lsp' },
                    { name = 'plugins' },
                    { name = 'path', option = { trailing_slash = true } },
                    { name = 'nvim_lsp_signature_help' },
                    { name = 'luasnip', options = { show_autosnippets = true } },
                    { name = 'nerdfont', },
                    { name = 'calc' },
                    { name = 'buffer', group_index = 2 },
                }
            })

            cmp.setup.filetype('zsh', {
                sources = {
                    { name = 'zsh' },
                    { name = 'nvim_lsp' },
                    { name = 'path', option = { trailing_slash = true } },
                    { name = 'nvim_lsp_signature_help' },
                    { name = 'luasnip', options = { show_autosnippets = true } },
                    { name = 'calc' },
                    { name = 'buffer', group_index = 2 },
                }
            })

            cmp.setup.filetype('julia', {
                sources = {
                    { name = 'plotlyjs' },
                    { name = 'nvim_lsp' },
                    { name = 'path', option = { trailing_slash = true } },
                    { name = 'nvim_lsp_signature_help' },
                    { name = 'luasnip', options = { show_autosnippets = true } },
                    -- { name = 'luasnip_choice' },
                    { name = 'calc' },
                    { name = 'buffer', group_index = 2 },
                }
            })

            cmp.setup.filetype('python', {
                sources = {
                    { name = 'pymol_settings' },
                    { name = 'nvim_lsp' },
                    { name = 'path', option = { trailing_slash = true } },
                    { name = 'nvim_lsp_signature_help' },
                    { name = 'luasnip', options = { show_autosnippets = true } },
                    -- { name = 'luasnip_choice' },
                    { name = 'calc' },
                    { name = 'buffer', group_index = 2 },
                }
            })

            cmp.setup.cmdline(':', {
                sources = {
                    { name = 'colorschemes' },
                }
            })

            -- local cmp_case = require "completecase"
            -- cmp.event:on('confirm_done', cmp_case.on_menu_closed())

            -- TODO: take inspo from https://github.com/hrsh7th/nvim-cmp/wiki/Menu-Appearance
            -- and link completion menu colors to equivalent things
        end
    },
    {
        'saadparwaiz1/cmp_luasnip',
        -- enabled = not using_blink, -- still needed
        lazy = true,
        dependencies = {
            'L3MON4D3/LuaSnip',
            -- "hrsh7th/nvim-cmp", -- magazine instead
        },
        config = function()
            local luasnip = require "luasnip"
            local cmp = require "cmp"
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

            -- the rest only if using pure cmp
            if using_blink then return end

            -- like pressing > which looks like forward arrow
            vim.keymap.set({ "i", "s", "n" }, "<C-.>", function()
                -- including expand means this will autocomplete the first visible snippet in completion menu
                -- if luasnip.expand_or_jumpable() then luasnip.expand_or_jump() end
                if luasnip.jumpable(1) then luasnip.jump(1) end
            end, { silent = true })

            -- like pressing < which looks like backwards arrow
            vim.keymap.set({ "i", "s", "n" }, "<C-,>", function()
                if luasnip.jumpable(-1) then luasnip.jump(-1) end
            end, { silent = true })

            -- like pressing ? for choices
            vim.keymap.set({ "i", "s", "n" }, "<C-/>", function()
                if luasnip.choice_active() then
                    luasnip.change_choice(1)
                elseif cmp.visible() then
                    -- This will change an active popup menu listing to only show snippets
                    cmp.complete { config = { sources = { { name = "luasnip" } } } }
                    -- Confirm first selection. Could also use the above by itself for filtering only.
                    cmp.confirm({ select = true })
                end
            end)
            -- with shift to go backwards
            vim.keymap.set({ "i", "s", "n" }, "<C-?>", function()
                if luasnip.choice_active() then
                    luasnip.change_choice(-1)
                end
            end)
        end
    },
}
