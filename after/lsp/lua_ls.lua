return {
    settings = {
        Lua = {
            diagnostics = {
                globals = { 'vim', 'lfs' },
                severity = {
                    ["inject-field"] = "Warning",
                },
            },
            completion = {
                autoRequire = false,
                callSnippet = "Replace",
            },
            workspace = {
                checkThirdParty = false,
            },
            runtime = {
                version = "LuaJIT",
            },
        }
    }
}
