
return {
    settings = {
        Lua = {
            diagnostics = {
                -- https://www.reddit.com/r/neovim/comments/11mdxex/lspconfig_cannot_access_configuration_for_lua_ls/
                -- Get the language server to recognize the `vim` global
                globals = { 'vim' },
            },
            workspace = {
                -- https://www.reddit.com/r/neovim/comments/11mdxex/lspconfig_cannot_access_configuration_for_lua_ls/
                -- Make the server aware of Neovim runtime files
                library = vim.api.nvim_get_runtime_file("", true),
                -- stopped a prompt at startup asking to configure project as luv.
                checkThirdParty = false
            },
            completion = { callSnippet = "Replace" },
        }
    }
}
