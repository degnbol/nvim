---@diagnostic disable: undefined-global
local compile_db = require("compile_db")

--- Make a project tree under a new temp dir.
--- @param files table<string, string[]> path relative to the root → file lines
--- @return string root realpath of the project root
local function project(files)
    local root = vim.fn.tempname()
    for path, lines in pairs(files) do
        local full = vim.fs.joinpath(root, path)
        vim.fn.mkdir(vim.fs.dirname(full), "p")
        vim.fn.writefile(lines, full)
    end
    vim.fn.mkdir(root, "p")
    return assert(vim.uv.fs_realpath(root))
end

--- Every file under a dir, relative to it.
--- @param root string
--- @return string[] paths sorted
local function tree(root)
    local paths = {}
    for name, type in vim.fs.dir(root, { depth = math.huge }) do
        if type == "file" then table.insert(paths, name) end
    end
    table.sort(paths)
    return paths
end

describe("compile_db.find_own", function()
    it("finds a database upward from the file", function()
        local root = project({ ["compile_commands.json"] = { "[]" }, ["src/a/x.c"] = {} })
        assert.are.equal(vim.fs.joinpath(root, "compile_commands.json"),
            compile_db.find_own(vim.fs.joinpath(root, "src/a/x.c"), root))
    end)

    it("finds build/compile_commands.json in an ancestor", function()
        local root = project({ ["build/compile_commands.json"] = { "[]" }, ["src/x.c"] = {} })
        assert.are.equal(vim.fs.joinpath(root, "build/compile_commands.json"),
            compile_db.find_own(vim.fs.joinpath(root, "src/x.c"), root))
    end)

    it("finds a database one level below the root", function()
        local root = project({ ["out/compile_commands.json"] = { "[]" }, ["src/x.c"] = {} })
        assert.are.equal(vim.fs.joinpath(root, "out/compile_commands.json"),
            compile_db.find_own(vim.fs.joinpath(root, "src/x.c"), root))
    end)

    it("finds compile_flags.txt", function()
        local root = project({ ["compile_flags.txt"] = { "-Iinc" }, ["x.c"] = {} })
        assert.are.equal(vim.fs.joinpath(root, "compile_flags.txt"),
            compile_db.find_own(vim.fs.joinpath(root, "x.c"), root))
    end)

    it("finds nothing in a project without one", function()
        local root = project({ ["src/x.c"] = {}, ["a/b/compile_commands.json"] = { "[]" } })
        assert.is_nil(compile_db.find_own(vim.fs.joinpath(root, "src/x.c"), root))
    end)
end)

describe("compile_db.collect_files", function()
    it("keeps C/C++ sources and headers outside dot-dirs and build dirs", function()
        local root = project({
            ["a.c"] = {}, ["b.cpp"] = {}, ["sub/c.cc"] = {}, ["inc/d.h"] = {}, ["inc/e.hpp"] = {},
            ["notes.txt"] = {}, ["x.m"] = {},
            [".git/f.c"] = {}, ["build/g.c"] = {}, ["build-debug/h.c"] = {}, ["build_rel/i.h"] = {},
        })
        local files = compile_db.collect_files(root)
        table.sort(files.sources)
        table.sort(files.headers)
        assert.are.same({
            vim.fs.joinpath(root, "a.c"), vim.fs.joinpath(root, "b.cpp"), vim.fs.joinpath(root, "sub/c.cc"),
        }, files.sources)
        assert.are.same({ vim.fs.joinpath(root, "inc/d.h"), vim.fs.joinpath(root, "inc/e.hpp") }, files.headers)
    end)

    it("stops at the entry cap with a warning", function()
        local root = vim.fn.tempname()
        vim.fn.mkdir(root, "p")
        for i = 1, 20010 do
            local file = io.open(vim.fs.joinpath(root, i .. ".c"), "w")
            assert(file):close()
        end
        local messages = {}
        local notify = vim.notify
        vim.notify = function(msg) table.insert(messages, msg) end
        local ok, files = pcall(compile_db.collect_files, root)
        vim.notify = notify
        assert(ok, files)
        assert.are.equal(20000, #files.sources)
        assert.are.equal(1, #messages)
    end)
end)

--- Run an async function and wait for its result.
--- @param fn async fun(): any
--- @return any result
local function run(fn)
    return require("vim._async").run(fn):wait()
end

describe("compile_db.infer_flags", function()
    local has_glib = vim.system({ "pkg-config", "--exists", "glib-2.0" }):wait().code == 0

    --- Flags inferred for a project.
    --- @param root string
    --- @return string[] flags
    local function flags_of(root)
        return run(function() return compile_db.infer_flags(compile_db.collect_files(root)) end)
    end

    it("puts project include dirs before pkg-config flags", function()
        if not has_glib then return pending("glib not installed") end
        local root = project({
            ["src/main.c"] = { "#include <config.h>", "#include <glib.h>", '#include "local.h"', "#include <stdio.h>" },
            ["src/local.h"] = {},
            ["gen/config.h"] = {},
        })
        local flags = flags_of(root)
        assert.are.equal("-I" .. vim.fs.joinpath(root, "gen"), flags[1])
        assert.are.same(require("utils/pkg_config").cflags({ "glib-2.0" }), vim.list_slice(flags, 2))
    end)

    it("gives an angle include next to its includer a -I, and a quoted one none", function()
        local root = project({
            ["src/main.c"] = { "#include <config.h>", '#include "local.h"' },
            ["src/config.h"] = {},
            ["src/local.h"] = {},
        })
        assert.are.same({ "-I" .. vim.fs.joinpath(root, "src") }, flags_of(root))
    end)
end)

describe("compile_db.generate", function()
    it("writes one entry per source outside the project, once for concurrent calls", function()
        local root = project({ ["a.c"] = { '#include "a.h"' }, ["a.h"] = {}, ["sub/b.c"] = {} })
        local tasks = {}
        for i = 1, 2 do tasks[i] = require("vim._async").run(function() return compile_db.generate(root) end) end
        local results = vim.tbl_map(function(task) return task:wait() end, tasks)
        assert.are.equal(results[1], results[2])
        assert.are.same({ "a.c", "a.h", "sub/b.c" }, tree(root))

        local dir = results[1].dir
        assert.are.equal(compile_db.cache_dir(root), dir)
        local entries = vim.json.decode(table.concat(vim.fn.readfile(vim.fs.joinpath(dir, "compile_commands.json")), "\n"))
        local files = vim.iter(entries):map(function(entry) return entry.file end):totable()
        table.sort(files)
        assert.are.same({ vim.fs.joinpath(root, "a.c"), vim.fs.joinpath(root, "sub/b.c") }, files)
        assert.are.equal(root, entries[1].directory)
        assert.are.equal("clang", entries[1].arguments[1])
        assert.are.same(results[1], compile_db.result(root))
    end)
end)

describe("compile_db.root_dir", function()
    vim.lsp.config("compile_db_test", { root_markers = { ".git" } })

    --- The root that `root_dir` passes to `on_dir` for a file, once any generation is done.
    --- @param path string
    --- @return string|nil root
    local function root_of(path)
        local buf = vim.fn.bufadd(path)
        local root, done = nil, false
        compile_db.root_dir("compile_db_test", buf, function(dir) root, done = dir, true end)
        vim.wait(60000, function() return done end)
        return root
    end

    it("generates a database, and forgets it once the project has its own", function()
        local root = project({ [".git/HEAD"] = {}, ["a.c"] = {} })
        local path = vim.fs.joinpath(root, "a.c")
        assert.are.equal(root, root_of(path))
        assert.truthy(compile_db.result(root))

        vim.fn.writefile({ "[]" }, vim.fs.joinpath(root, "compile_commands.json"))
        assert.are.equal(root, root_of(path))
        assert.is_nil(compile_db.result(root))
    end)

    it("reports the root after a failed generation, with an error", function()
        local root = project({ [".git/HEAD"] = {}, ["a.c"] = {} })
        local write, notify = compile_db.write, vim.notify
        local messages = {}
        compile_db.write = function() error("disk full", 0) end
        vim.notify = function(msg) table.insert(messages, msg) end
        local ok, reported = pcall(root_of, vim.fs.joinpath(root, "a.c"))
        compile_db.write, vim.notify = write, notify
        assert(ok, reported)
        assert.are.equal(root, reported)
        assert.is_nil(compile_db.result(root))
        assert.are.same({ "compile_db: disk full" }, messages)
    end)
end)
