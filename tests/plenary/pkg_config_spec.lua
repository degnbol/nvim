---@diagnostic disable: undefined-global
local cflags = require("utils/pkg_config").cflags

--- Run fn with vim.notify replaced by a recorder.
--- @param fn function
--- @return any result return value of fn
--- @return string[] messages each message passed to vim.notify
local function capture_notify(fn)
    local messages = {}
    local notify = vim.notify
    vim.notify = function(msg) table.insert(messages, msg) end
    local ok, result = pcall(fn)
    vim.notify = notify
    assert(ok, result)
    return result, messages
end

describe("pkg_config.cflags", function()
    local has_glib = vim.system({ "pkg-config", "--exists", "glib-2.0", "gobject-2.0" }):wait().code == 0

    it("returns the package's include directory as separate flags", function()
        if not has_glib then return pending("glib not installed") end
        local flags = cflags({ "glib-2.0" })
        assert.truthy(vim.iter(flags):any(function(flag) return flag:match("^%-I.*glib%-2%.0$") end))
        assert.is_nil(vim.iter(flags):find(function(flag) return flag:find("%s") end))
    end)

    it("lists a flag shared by two packages once", function()
        if not has_glib then return pending("glib not installed") end
        local flags = cflags({ "glib-2.0", "gobject-2.0" })
        assert.are.same(vim.list.unique(vim.deepcopy(flags)), flags)
    end)

    it("skips an unknown package with a warning and keeps the others", function()
        if not has_glib then return pending("glib not installed") end
        local flags, messages = capture_notify(function() return cflags({ "no-such-package-xyz", "glib-2.0" }) end)
        assert.are.same(cflags({ "glib-2.0" }), flags)
        assert.are.equal(1, #messages)
        assert.truthy(messages[1]:find("no-such-package-xyz", 1, true))
    end)

    it("returns no flags, with a warning, when pkg-config is not installed", function()
        local path = vim.env.PATH
        vim.env.PATH = ""
        local ok, flags, messages = pcall(capture_notify, function() return cflags({ "glib-2.0" }) end)
        vim.env.PATH = path
        assert(ok, flags)
        assert.are.same({}, flags)
        assert.are.equal(1, #messages)
        assert.truthy(messages[1]:find("pkg-config not found", 1, true))
    end)

    it("returns no flags for no packages", function()
        assert.are.same({}, cflags({}))
    end)
end)

describe("pkg_config.providers", function()
    local providers = require("utils/pkg_config").providers
    local tmp = vim.fn.tempname()
    for _, dir in ipairs({ "inc", "inc/foo", "other", "shared/SDL2" }) do vim.fn.mkdir(vim.fs.joinpath(tmp, dir), "p") end
    vim.fn.writefile({}, vim.fs.joinpath(tmp, "inc", "foo", "foo.h"))
    vim.fn.writefile({}, vim.fs.joinpath(tmp, "other", "bar.h"))
    vim.fn.writefile({}, vim.fs.joinpath(tmp, "shared", "SDL2", "SDL.h"))
    vim.fn.writefile({}, vim.fs.joinpath(tmp, "shared", "ssl.h"))
    local index = {
        broad = { vim.fs.joinpath(tmp, "inc") },
        foo = { vim.fs.joinpath(tmp, "inc", "foo") },
        zbar = { vim.fs.joinpath(tmp, "other") },
        abar = { vim.fs.joinpath(tmp, "other") },
        Qt6Platform = { vim.fs.joinpath(tmp, "shared") },
        sdl2 = { vim.fs.joinpath(tmp, "shared") },
        libssl = { vim.fs.joinpath(tmp, "shared") },
    }

    it("maps each header to the package whose own dir holds it, skipping unowned ones", function()
        assert.are.same({ "abar", "foo" }, providers({ "foo.h", "bar.h", "none.h" }, index))
    end)

    it("prefers the longest dir, then a package named after the header, then the first name", function()
        assert.are.same({ "broad" }, providers({ "foo/foo.h" }, index))
        assert.are.same({ "foo" }, providers({ "foo.h" }, index))
        assert.are.same({ "abar" }, providers({ "bar.h" }, index))
        assert.are.same({ "sdl2" }, providers({ "SDL2/SDL.h" }, index))
        assert.are.same({ "libssl" }, providers({ "ssl.h" }, index))
    end)
end)

describe("pkg_config.include_dir_index", function()
    local include_dir_index = require("utils/pkg_config").include_dir_index
    local has_glib = vim.system({ "pkg-config", "--exists", "glib-2.0" }):wait().code == 0

    --- Run fn with vim.system replaced by a recorder.
    --- @param fn function
    --- @return any result return value of fn
    --- @return string[][] cmds each command passed to vim.system
    local function capture_system(fn)
        local cmds = {}
        local system = vim.system
        vim.system = function(cmd, ...) table.insert(cmds, cmd); return system(cmd, ...) end
        local ok, result = pcall(fn)
        vim.system = system
        assert(ok, result)
        return result, cmds
    end

    --- Start include_dir_index calls at once and wait for all of them.
    --- @param n_calls integer
    --- @return table<string, string[]>[] indexes one per call
    local function indexes_now(n_calls)
        local tasks = {}
        for i = 1, n_calls do tasks[i] = require("vim._async").run(include_dir_index) end
        return vim.tbl_map(function(task) return task:wait() end, tasks)
    end

    it("builds once for concurrent calls, and gives glib its own dirs only", function()
        if not has_glib then return pending("glib not installed") end
        local indexes, cmds = capture_system(function() return indexes_now(2) end)
        local n_lists = #vim.tbl_filter(function(cmd) return cmd[2] == "--list-package-names" end, cmds)
        assert.are.equal(1, n_lists)
        assert.are.same(indexes[1], indexes[2])
        local glib_dirs = indexes[1]["glib-2.0"]
        assert.truthy(vim.iter(glib_dirs):any(function(dir) return dir:match("include/glib%-2%.0$") end))
        assert.is_nil(vim.iter(glib_dirs):find(function(dir) return dir:find("pcre2", 1, true) end))
    end)

    it("serves a later call from the cache", function()
        if not has_glib then return pending("glib not installed") end
        local first = indexes_now(1)[1]
        local indexes, cmds = capture_system(function() return indexes_now(1) end)
        assert.are.same(first, indexes[1])
        assert.are.equal(1, #cmds) -- the pc_path lookup for the cache key
    end)
end)
