local async = require "utils/async"
local util = require "utils/init"
local vim_async = require "vim._async"

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

local max_processes = 16

--- The whitespace-separated words of a command's stdout.
--- @async
--- @param cmd string[]
--- @return string[]|nil words nil on a non-zero exit
--- @return string stderr trimmed
local function words(cmd)
    local result = async.system(cmd, { text = true })
    if result.code ~= 0 then return nil, vim.trim(result.stderr) end
    return vim.split(vim.trim(result.stdout), "%s+", { trimempty = true }), ""
end

--- Cache key of the installed `.pc` files: the search path variables and the mtimes of
--- the search dirs. A `.pc` edited in place does not change it.
--- @async
--- @return string key
local function index_key()
    local path, libdir = vim.env.PKG_CONFIG_PATH, vim.env.PKG_CONFIG_LIBDIR
    local pc_path = libdir or (words({ "pkg-config", "--variable", "pc_path", "pkg-config" }) or {})[1] or ""
    local parts = { path or "", libdir or "" }
    for _, dir in ipairs(vim.split((path or "") .. ":" .. pc_path, ":", { trimempty = true })) do
        local stat = vim.uv.fs_stat(dir)
        table.insert(parts, dir .. "=" .. (stat and stat.mtime.sec .. "." .. stat.mtime.nsec or "-"))
    end
    return table.concat(parts, "\n")
end

--- Each installed package's own include dirs, from at most `max_processes` pkg-config
--- runs at a time.
--- @async
--- @return table<string, string[]> index package name → include dirs
--- @return string[] failures one message per package pkg-config failed on; those are left out
local function build_index()
    local packages, stderr = words({ "pkg-config", "--list-package-names" })
    if not packages then return {}, { stderr } end
    local index, failures = {}, {}
    local jobs = vim.tbl_map(function(package)
        return function()
            -- Stays set if the job raises: vim._async.join drops errors of its jobs.
            failures[package] = package .. ": did not run"
            local flags, err = words({ "pkg-config", "--maximum-traverse-depth=1", "--cflags-only-I", package })
            if flags then
                failures[package] = nil
                index[package] = vim.tbl_map(function(flag) return (flag:gsub("^%-I", "")) end, flags)
            else
                failures[package] = package .. ": " .. err
            end
        end
    end, packages)
    vim_async.join(max_processes, jobs)
    return index, vim.tbl_values(failures)
end

--- @return string path
local function index_cache_path()
    return vim.fs.joinpath(vim.fn.stdpath("cache"), "pkg_config_include_dirs.json")
end

--- The index for a cache key: the cached one if its key matches, else a new one, cached
--- when pkg-config failed on no package.
--- @async
--- @param key string see `index_key`
--- @return table<string, string[]> index
local function load_index(key)
    local text = util.readtext(index_cache_path())
    if text then
        local ok, cached = pcall(vim.json.decode, text)
        if ok and type(cached) == "table" and cached.key == key then return cached.index end
    end
    local index, failures = build_index()
    if #failures > 0 then
        vim.notify("pkg-config failed, index not cached:\n" .. table.concat(failures, "\n"), vim.log.levels.WARN)
        return index
    end
    vim.fn.mkdir(vim.fn.stdpath("cache"), "p")
    local out, err = io.open(index_cache_path(), "w")
    if not out then
        vim.notify("pkg-config index not cached: " .. err, vim.log.levels.WARN)
        return index
    end
    out:write(vim.json.encode({ key = key, index = index }))
    out:close()
    return index
end

local load_index_shared = async.shared(load_index)

--- Each installed package's *own* include dirs: its `-I` flags, without those of the
--- packages it requires. Cached on disk until the `.pc` search path or its dirs change.
--- Concurrent calls share one build. If pkg-config is not installed, warns once and
--- reports no packages.
--- @async
--- @return table<string, string[]> index package name → include dirs
function M.include_dir_index()
    if vim.fn.executable("pkg-config") == 0 then
        vim.notify_once("pkg-config not found; no include dirs from packages", vim.log.levels.WARN)
        return {}
    end
    return load_index_shared(index_key())
end

--- Whether a package is named after a header: its name, without a version suffix
--- and ignoring case, is the header's first path component or its stem, with or
--- without a "lib" prefix. E.g. "sdl2" for "SDL2/SDL.h", "glib-2.0" for "glib.h".
--- @param package string package name
--- @param header string header name as written in `#include`
--- @return boolean
local function is_named_after(package, header)
    local name = package:lower():gsub("%-[%d.]+$", "")
    local first = (header:lower():match("^[^/]+") or ""):gsub("%.[^.]*$", "")
    local stem = vim.fs.basename(header):lower():gsub("%.[^.]*$", "")
    return vim.list_contains({ first, stem, "lib" .. first, "lib" .. stem }, name)
end

--- The packages that provide a set of headers. Per header, this is the package with the
--- longest include dir that holds the header; on a tie, a package named after the header
--- (see `is_named_after`), then the first by name. Headers that no package holds are skipped.
--- @param headers string[] header names as written in `#include`, e.g. "glib.h"
--- @param index table<string, string[]> package name → include dirs of the package itself
--- @return string[] packages sorted, each once
function M.providers(headers, index)
    local packages = vim.tbl_keys(index)
    table.sort(packages)
    local providers = {}
    for _, header in ipairs(headers) do
        local best -- { package, dir, named }
        for _, package in ipairs(packages) do
            for _, dir in ipairs(index[package]) do
                if vim.uv.fs_stat(vim.fs.joinpath(dir, header)) then
                    local named = is_named_after(package, header)
                    if not best or #dir > #best.dir or (#dir == #best.dir and named and not best.named) then
                        best = { package = package, dir = dir, named = named }
                    end
                end
            end
        end
        if best then table.insert(providers, best.package) end
    end
    table.sort(providers)
    return vim.list.unique(providers)
end

return M
