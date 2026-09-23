local async = require "utils/async"

local M = {}

--- Header names of the `#include` directives in C/C++ source lines.
--- A directive inside a block comment is reported too.
--- @param lines string[]
--- @return { angle: string[], quoted: string[] } headers of `#include <…>` and `#include "…"`, in order, each once
function M.includes(lines)
    local angle, quoted = {}, {}
    for _, line in ipairs(lines) do
        local directive = line:match("^%s*#%s*include%s*(.*)")
        if directive then
            table.insert(angle, directive:match("^<([^>]+)>"))
            table.insert(quoted, directive:match('^"([^"]+)"'))
        end
    end
    return { angle = vim.list.unique(angle), quoted = vim.list.unique(quoted) }
end

--- The `#include <…>` search dirs listed in a compiler's `-v` output, without framework dirs.
--- @param verbose_output string stderr of e.g. `clang -E -v -`
--- @return string[] dirs in search order
function M.parse_search_dirs(verbose_output)
    local dirs = {}
    local inside = false
    for line in vim.gsplit(verbose_output, "\n") do
        if vim.startswith(line, "#include <...> search starts here:") then
            inside = true
        elseif vim.startswith(line, "End of search list.") then
            break
        elseif inside and not line:find("%(framework directory%)$") then
            table.insert(dirs, vim.trim(line))
        end
    end
    return dirs
end

local warned_missing = {}

--- A compiler's default `#include <…>` search dirs. If the compiler is not installed,
--- warns once per compiler and reports no dirs.
--- @async
--- @param compiler string executable, e.g. "clang"
--- @param language "c"|"c++"
--- @return string[] dirs in search order
function M.search_dirs(compiler, language)
    if vim.fn.executable(compiler) == 0 then
        if not warned_missing[compiler] then
            warned_missing[compiler] = true
            vim.notify(compiler .. " not found; no default include dirs", vim.log.levels.WARN)
        end
        return {}
    end
    local result = async.system({ compiler, "-x" .. language, "-E", "-v", "-" }, { stdin = "", text = true })
    return M.parse_search_dirs(result.stderr)
end

--- Where a header resolves in a list of include dirs.
--- @param header string name as written in the directive, e.g. "glib/gtypes.h"
--- @param dirs string[] include dirs, in search order
--- @return string|nil path of the header in the first dir that holds it
function M.find_header(header, dirs)
    for _, dir in ipairs(dirs) do
        local path = vim.fs.joinpath(dir, header)
        if vim.uv.fs_stat(path) then return path end
    end
end

return M
