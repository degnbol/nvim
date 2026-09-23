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
