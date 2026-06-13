local M = {}

--- Resolve a plugin's web URL from package-spec source strings.
--- Matches the spec entry whose source path ends in `/<repo>`. A full URL is
--- returned verbatim (honouring non-github remotes), a bare `account/repo`
--- slug is treated as a github shorthand (what `gh()` expands to).
--- @param repo string bare repository name, without account or host
--- @param lines string[] lines of the package spec file
--- @return string|nil url
function M.resolve(repo, lines)
    local pat = [["([^"]*/]] .. vim.pesc(repo) .. [[)"]]
    for _, line in ipairs(lines) do
        local src = line:match(pat)
        if src then
            if src:match("^https?://") then
                return src
            end
            return "https://github.com/" .. src
        end
    end
end

return M
