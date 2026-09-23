local M = {}

--- Compile flags (include paths, defines) that pkg-config reports for packages.
--- Packages pkg-config fails on are left out, with one warning that gives pkg-config's
--- error for each. If pkg-config is not installed, returns no flags, with a warning.
--- @param packages string[] pkg-config package names, e.g. "glib-2.0"
--- @return string[] flags in package order, each listed once
function M.cflags(packages)
    if #packages > 0 and vim.fn.executable("pkg-config") == 0 then
        vim.notify("pkg-config not found; no flags for " .. table.concat(packages, ", "), vim.log.levels.WARN)
        return {}
    end
    local flags, failures = {}, {}
    for _, name in ipairs(packages) do
        local result = vim.system({ "pkg-config", "--cflags", name }, { text = true }):wait()
        if result.code == 0 then
            vim.list_extend(flags, vim.split(vim.trim(result.stdout), "%s+", { trimempty = true }))
        else
            table.insert(failures, name .. ": " .. vim.trim(result.stderr))
        end
    end
    if #failures > 0 then
        vim.notify("pkg-config failed:\n" .. table.concat(failures, "\n"), vim.log.levels.WARN)
    end
    return vim.list.unique(flags)
end

return M
