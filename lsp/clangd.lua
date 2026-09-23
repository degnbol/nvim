return {
    -- `compile_db` is required at call time: vim.lsp.enable loads this file on every startup.
    root_dir = function(bufnr, on_dir) require("compile_db").root_dir("clangd", bufnr, on_dir) end,
    before_init = function(params, config)
        local result = config.root_dir and require("compile_db").result(config.root_dir)
        if not result then return end
        params.initializationOptions = vim.tbl_extend("force", params.initializationOptions or {}, {
            compilationDatabasePath = result.dir,
            -- Used only when the database has no entries, e.g. for a header-only project:
            -- clangd infers the command of a file without an entry from the nearest entry.
            fallbackFlags = result.flags,
        })
    end,
}
