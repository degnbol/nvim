local M = {}

--- Add every local plugin under <root>/modules to the runtimepath.
--- @param root string
--- @return string[] bases Plugin directories, in rtp order.
function M.add(root)
    --- @type string[]
    local bases = {}
    for mod, kind in vim.fs.dir(root .. "/modules") do
        if kind == "directory" then
            local base = root .. "/modules/" .. mod
            vim.opt.runtimepath:append(base)
            local after_dir = base .. "/after"
            if vim.uv.fs_stat(after_dir) then vim.opt.runtimepath:append(after_dir) end
            table.insert(bases, base)
        end
    end
    return bases
end

return M
