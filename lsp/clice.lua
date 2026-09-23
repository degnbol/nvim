--- Whether a directory has a clice config file, at either path clice reads.
--- @param root string workspace root directory
--- @return boolean
local function has_clice_config(root)
    return vim.uv.fs_stat(vim.fs.joinpath(root, "clice.toml")) ~= nil
        or vim.uv.fs_stat(vim.fs.joinpath(root, ".clice", "config.toml")) ~= nil
end

return {
    -- `compile_db` is required at call time: vim.lsp.enable loads this file on every startup.
    root_dir = function(bufnr, on_dir) require("compile_db").root_dir("clice", bufnr, on_dir) end,
    -- Skipped when the project has its own config: initializationOptions `rules`
    -- replace the rules in clice.toml instead of merging with them.
    before_init = function(params, config)
        if not config.root_dir or has_clice_config(config.root_dir) then return end
        local options = {
            project = { cache_dir = vim.fs.joinpath(vim.fn.stdpath("cache"), "clice", vim.fn.sha256(config.root_dir)) },
        }
        local result = require("compile_db").result(config.root_dir)
        if result then
            -- `default_command` applies only to files without a compilation database entry.
            options.rules = { {
                compile_commands = { result.dir },
                default_command = vim.list_extend({ "clang" }, result.flags),
            } }
        end
        params.initializationOptions = vim.tbl_extend("force", params.initializationOptions or {}, options)
    end,
}
