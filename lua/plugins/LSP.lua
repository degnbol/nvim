#!/usr/bin/env lua
return {
    -- LSP. For a given file, either complete with cmp (builtin recommended, but e.g. jedi language servers is slow), coc (old, not using builtin LSP), or coq (hacks builtin LSP)
    {"neovim/nvim-lspconfig", config=function() require'lsp' end},
    "folke/neodev.nvim", -- signature on nvim lua calls which helps in messing with plugins etc
    -- add :LspInstall <language> and :Mason for conveniently installing LSP language specific servers
    {
        "williamboman/mason-lspconfig.nvim",
        dependencies={"neovim/nvim-lspconfig", "williamboman/mason.nvim"},
        config=function() 
            -- https://github.com/williamboman/mason-lspconfig.nvim
            -- naming: https://github.com/williamboman/mason-lspconfig.nvim/blob/main/doc/server-mapping.md
            -- config help: https://github.com/neovim/nvim-lspconfig/blob/master/doc/server_configurations.md
            -- also see lsp.lua
            require("mason").setup()
            require("mason-lspconfig").setup {
                ensure_installed = {
                    "bashls",
                    -- https://github.com/neovim/nvim-lspconfig/blob/master/doc/server_configurations.md#omnisharp
                    "csharp_ls",
                    "jedi_language_server",
                    "pyright",
                    "pylsp",
                    "julials",
                    "ltex",
                    "texlab",
                    "r_language_server",
                    "lua_ls",
                    "vimls",
                    "rust_analyzer",
                    "tsserver", -- javascript. MasonInstall typescript-language-server
                }
            }
        end
    },
    -- completion menu using builtin LSP
    {"hrsh7th/nvim-cmp", dependencies = {
        'hrsh7th/cmp-nvim-lsp', 'hrsh7th/cmp-buffer',
        'hrsh7th/cmp-path',
        'hrsh7th/cmp-nvim-lsp-signature-help',
        'hrsh7th/cmp-nvim-lua', -- neovim Lua API
        'tamago324/cmp-zsh', -- neovim zsh completion
        'onsails/lspkind.nvim', -- pretty pictograms
        'hrsh7th/cmp-calc', -- quick math in completion
    }, config=function() require'cmp-conf' end},
    {'saadparwaiz1/cmp_luasnip', dependencies={'L3MON4D3/LuaSnip', "hrsh7th/nvim-cmp"}, config=function() require"luasnip-conf" end }, -- alts: hrsh7th/vim-vsnip, SirVer/ultisnips, ...
    -- custom dicts and spell check that doesn't require spell and spelllang (f3fora/cmp-spell)
    {'uga-rosa/cmp-dictionary', config=function() require"cmp_dict" end},
    {"j-hui/fidget.nvim", config=function() require"fidget-conf" end}, -- corner print what LSP is running
    "rafamadriz/friendly-snippets",
}

