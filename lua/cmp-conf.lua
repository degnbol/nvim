#!/usr/bin/env lua
-- inspiration from https://vonheikemen.github.io/devlog/tools/setup-nvim-lspconfig-plus-nvim-cmp/
local cmd = vim.cmd
local cmp = require "cmp"
-- https://github.com/onsails/lspkind.nvim
local lspkind = require "lspkind"
local luasnip = require "luasnip"
local rtp = vim.opt.runtimepath:get()[1]

-- menu=show completion menu. menuone=also when only one option. noselect=don't select automatically.
vim.opt.completeopt = {"menu", "menuone", "noselect"}

-- https://youtu.be/Dn800rlPIho?t=440
luasnip.config.set_config {
    -- don't jump back into exited snippet
    history = false,
    -- dynamic snippets update as you type
    updateevents = "TextChanged,TextChangedI",
    enable_autosnippets = true,
}
-- load friendly-snippets with luasnip
require("luasnip.loaders.from_vscode").lazy_load()

require("cmp_dictionary").setup {
    dic = {
        -- dicts generated with ./spell.sh
        ["*"] = {
            rtp .. "/spell/custom.dic",
            rtp .. "/spell/en.dic",
            rtp .. "/spell/da.dic"},
    },
    first_case_insensitive = true,
}

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
        expand = function(args) require('luasnip').lsp_expand(args.body) end,
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
            elseif luasnip.expand_or_jumpable() then
                luasnip.expand_or_jump()
            elseif has_words_before() then
                cmp.complete()
            else
                fallback()
            end
        end, { "i", "s" }),
        ["<S-Tab>"] = cmp.mapping(function(fallback)
            if cmp.visible() then
                cmp.select_prev_item()
            elseif luasnip.jumpable(-1) then
                luasnip.jump(-1)
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
    sources = cmp.config.sources {
        { name = 'nvim_lsp' },
        { name = 'path' },
        { name = 'nvim_lsp_signature_help' },
        { name = 'luasnip' },
        { name = 'calc' },
        { name = 'buffer' },
        { name = 'dictionary', keyword_length = 3 },
    }
})

cmp.setup.filetype('lua', {
    sources = cmp.config.sources {
        { name = 'nvim_lua' },
        { name = 'nvim_lsp' },
        { name = 'path' },
        { name = 'nvim_lsp_signature_help' },
        { name = 'luasnip' },
        { name = 'calc' },
        { name = 'buffer' },
    }
})
