#!/usr/bin/env lua
-- https://github.com/williamboman/mason-lspconfig.nvim
-- naming: https://github.com/williamboman/mason-lspconfig.nvim/blob/main/doc/server-mapping.md
-- config help: https://github.com/neovim/nvim-lspconfig/blob/master/doc/server_configurations.md
-- also see lsp.lua
require("mason").setup()
require("mason-lspconfig").setup {
    ensure_installed = {
        "bashls",
        "jedi_language_server",
        "pyright",
        "pylsp",
        "julials",
        "ltex",
        "texlab",
        "r_language_server",
        "sumneko_lua",
        "vimls",
    }
}
