-- A generated compile_commands.json for C/C++ projects without one, in the cache dir.
-- Server behaviour: .claude/skills/lsp/references/c-cpp.md
local async = require "utils/async"
local c_includes = require "utils/c_includes"
local pkg_config = require "utils/pkg_config"
local vim_async = require "vim._async"

local M = {}

local source_exts = { c = true, cc = true, cpp = true, cxx = true, ["c++"] = true }
local header_exts = { h = true, hh = true, hpp = true, hxx = true }
local max_entries = 20000
local files_per_chunk = 200

--- The compilation database a C/C++ file already has: `compile_commands.json`,
--- `build/compile_commands.json` or `compile_flags.txt` in the file's dir or an ancestor,
--- else `compile_commands.json` in an immediate subdir of `root`.
--- @param path string file
--- @param root string|nil workspace root, nil to skip the subdir search
--- @return string|nil path of the database, nil if there is none
function M.find_own(path, root)
    local found = vim.fs.find({ "compile_commands.json", "build/compile_commands.json", "compile_flags.txt" },
        { upward = true, path = vim.fs.dirname(path) })[1]
    if found or not root then return found end
    for name, type in vim.fs.dir(root) do
        local candidate = vim.fs.joinpath(root, name, "compile_commands.json")
        if type == "directory" and vim.uv.fs_stat(candidate) then return candidate end
    end
end

--- The generated database dir for a project.
--- @param root string project root
--- @return string dir
function M.cache_dir(root)
    return vim.fs.joinpath(vim.fn.stdpath("cache"), "compile_db", vim.fn.sha256(root))
end

--- @param name string dir name
--- @return boolean
local function is_build_or_hidden_dir(name)
    return vim.startswith(name, ".") or name == "build" or name:find("^build[-_]") ~= nil
end

--- C/C++ files of a project, outside dot-dirs and `build`, `build-*`, `build_*` dirs.
--- Stops after `max_entries` walked entries, with a warning.
--- @param root string project root
--- @return { sources: string[], headers: string[] } files absolute paths, in walk order
function M.collect_files(root)
    local files = { sources = {}, headers = {} }
    local n_entries = 0
    local walk = vim.fs.dir(root, {
        depth = math.huge,
        skip = function(dir) return not is_build_or_hidden_dir(vim.fs.basename(dir)) end,
    })
    for name, type in walk do
        n_entries = n_entries + 1
        if n_entries > max_entries then
            vim.notify(("compile_db: stopped after %d entries in %s"):format(max_entries, root), vim.log.levels.WARN)
            break
        end
        local ext = type == "file" and name:match("%.([^./]+)$")
        if source_exts[ext] then
            table.insert(files.sources, vim.fs.joinpath(root, name))
        elseif header_exts[ext] then
            table.insert(files.headers, vim.fs.joinpath(root, name))
        end
    end
    return files
end

--- The `#include`d headers of files. Reads `files_per_chunk` files per main loop
--- iteration, so the editor stays responsive. Files that cannot be read are left out,
--- with one warning.
--- @async
--- @param paths string[] files
--- @return table<string, { angle: string[], quoted: string[] }> includes path → headers, see `c_includes.includes`
local function read_includes(paths)
    local includes, unreadable = {}, {}
    for i, path in ipairs(paths) do
        local file, err = io.open(path)
        local content = file and file:read("*a")
        if file then file:close() end
        if content then
            includes[path] = c_includes.includes(vim.split(content, "\n"))
        else
            table.insert(unreadable, err or path)
        end
        if i % files_per_chunk == 0 then vim_async.await(1, vim.schedule) end
    end
    if #unreadable > 0 then
        vim.notify("compile_db: cannot read\n" .. table.concat(unreadable, "\n"), vim.log.levels.WARN)
    end
    return includes
end

--- Dirs from which a header name resolves to one of a set of files.
--- @param header string name as written in `#include`, e.g. "sub/x.h"
--- @param paths string[] candidate files
--- @return string[] dirs
local function dirs_holding(header, paths)
    local suffix = "/" .. header
    return vim.iter(paths)
        :filter(function(path) return vim.endswith(path, suffix) end)
        :map(function(path) return path:sub(1, -#suffix - 1) end)
        :totable()
end

--- Compile flags for a C/C++ project: `-I` for each project dir that holds an included
--- header (except a quoted one next to its includer), then the pkg-config flags of the
--- packages that provide the headers found neither in the project nor in clang's
--- default search dirs.
--- @async
--- @param files { sources: string[], headers: string[] } absolute paths of the project's files
--- @return string[] flags
function M.infer_flags(files)
    local paths = vim.list_extend(vim.list_extend({}, files.sources), files.headers)
    local project_dirs, outside = {}, {}
    for path, headers in pairs(read_includes(paths)) do
        -- Only `#include "…"` searches the includer's own dir.
        local unresolved = vim.tbl_filter(function(header)
            return not c_includes.find_header(header, { vim.fs.dirname(path) })
        end, headers.quoted)
        for _, header in ipairs(vim.list_extend(unresolved, headers.angle)) do
            local dirs = dirs_holding(header, files.headers)
            if #dirs > 0 then
                vim.list_extend(project_dirs, dirs)
            else
                table.insert(outside, header)
            end
        end
    end
    project_dirs = vim.list.unique(project_dirs)
    table.sort(project_dirs)

    local system_dirs = vim.list_extend(c_includes.search_dirs("clang", "c"), c_includes.search_dirs("clang", "c++"))
    local missing = vim.tbl_filter(function(header)
        return not c_includes.find_header(header, system_dirs)
    end, vim.list.unique(outside))
    local flags = vim.tbl_map(function(dir) return "-I" .. dir end, project_dirs)
    return vim.list_extend(flags, pkg_config.cflags(pkg_config.providers(missing, pkg_config.include_dir_index())))
end

--- Write a compile_commands.json with one `clang` command per source.
--- Raises an error if the file cannot be opened for writing.
--- @param dir string database dir, created if missing
--- @param root string working directory of each command
--- @param sources string[] absolute paths
--- @param flags string[] compile flags shared by all sources
function M.write(dir, root, sources, flags)
    vim.fn.mkdir(dir, "p")
    local entries = vim.tbl_map(function(file)
        local arguments = vim.list_extend(vim.list_extend({ "clang" }, flags), { "-c", file })
        return { directory = root, file = file, arguments = arguments }
    end, sources)
    local out = assert(io.open(vim.fs.joinpath(dir, "compile_commands.json"), "w"))
    out:write(vim.json.encode(entries))
    out:close()
end

--- A generated database: its dir, and the compile flags all its entries share.
--- @alias compile_db.Result { dir: string, flags: string[] }

--- Root → last generated database.
--- @type table<string, compile_db.Result>
local results = {}
--- Names of the LSP configs whose root is found by `root_dir`.
--- @type table<string, true>
local users = {}

--- Generate a project's database in `cache_dir(root)`, and keep the result for
--- `result(root)`, which a failure clears. Concurrent calls for a root share one run.
--- @async
--- @param root string project root
--- @return compile_db.Result result
M.generate = async.shared(function(root)
    results[root] = nil
    local files = M.collect_files(root)
    local flags = M.infer_flags(files)
    local dir = M.cache_dir(root)
    M.write(dir, root, files.sources, flags)
    results[root] = { dir = dir, flags = flags }
    return results[root]
end)

--- Run `generate`, and report a failure with a notification.
--- @param root string project root
--- @param callback fun(result: compile_db.Result|nil) nil on failure; called on the main loop
local function generate_then(root, callback)
    vim_async.run(function() return M.generate(root) end, function(err, result)
        vim.schedule(function()
            if err then vim.notify("compile_db: " .. err, vim.log.levels.ERROR) end
            callback(result)
        end)
    end)
end

--- An LSP `root_dir` function body. It finds the root from the config's `root_markers`.
--- If the buffer's file has no database of its own (`find_own`) and none was generated
--- for the root yet, it generates one before it reports the root, also after a failure.
--- If the file has its own database, it forgets the kept `result(root)`.
--- The root is not reported if the buffer is deleted during generation.
--- @param name string LSP config name
--- @param bufnr integer
--- @param on_dir fun(root: string|nil)
function M.root_dir(name, bufnr, on_dir)
    users[name] = true
    local root = vim.fs.root(bufnr, vim.lsp.config[name].root_markers)
    if not root then return on_dir(nil) end
    if M.find_own(vim.api.nvim_buf_get_name(bufnr), root) then
        results[root] = nil
        return on_dir(root)
    end
    if results[root] then return on_dir(root) end
    generate_then(root, function()
        if vim.api.nvim_buf_is_valid(bufnr) then on_dir(root) end
    end)
end

--- Generate the database again for the roots of a buffer's clients whose root is found
--- by `root_dir`, then restart those clients, so they read new `#include`s. A root with
--- its own database is skipped, with a warning.
--- @param bufnr integer
function M.regenerate(bufnr)
    local path = vim.api.nvim_buf_get_name(bufnr)
    local clients_by_root = {}
    for _, client in ipairs(vim.lsp.get_clients({ bufnr = bufnr })) do
        if users[client.name] and client.root_dir then
            clients_by_root[client.root_dir] = clients_by_root[client.root_dir] or {}
            table.insert(clients_by_root[client.root_dir], client)
        end
    end
    if next(clients_by_root) == nil then
        return vim.notify("compile_db: no client in this buffer uses a generated database", vim.log.levels.WARN)
    end
    for root, clients in pairs(clients_by_root) do
        local own = M.find_own(path, root)
        if own then
            vim.notify("compile_db: " .. root .. " has its own database " .. own, vim.log.levels.WARN)
        else
            generate_then(root, function(result)
                if not result then return end
                -- The function behind `:lsp restart`; there is no public Lua one.
                for _, client in ipairs(clients) do client:_restart() end
            end)
        end
    end
end

--- The last database generated for a root.
--- @param root string project root
--- @return compile_db.Result|nil result nil if none, or if the last generation failed
function M.result(root)
    return results[root]
end

return M
