--- Whether a directory has a clice config file, at either path clice reads.
--- @param root string workspace root directory
--- @return boolean
local function has_clice_config(root)
    return vim.uv.fs_stat(vim.fs.joinpath(root, "clice.toml")) ~= nil
        or vim.uv.fs_stat(vim.fs.joinpath(root, ".clice", "config.toml")) ~= nil
end

return {
    -- `default_command` applies only to files without a compilation database entry.
    -- Skipped when the project has its own config: initializationOptions `rules`
    -- replace the rules in clice.toml instead of merging with them.
    before_init = function(params, config)
        if config.root_dir and has_clice_config(config.root_dir) then return end
        params.initializationOptions = vim.tbl_extend("force", params.initializationOptions or {}, {
            rules = { { default_command = vim.list_extend({ "clang" }, require "c_fallback_flags") } },
        })
    end,
}
