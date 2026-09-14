---@diagnostic disable: undefined-global
-- Tests for lua/utils/paths.lua: resolve_path (pure single-candidate
-- resolution) and resolve_location_under_cursor (the gf / K resolver, with its
-- :lnum[:col] suffix parsing), plus the gf mapping in plugin/paths.lua that
-- consumes it.

local paths = require("utils.paths")

describe("utils.paths.resolve_path", function()
    it("returns an absolute token normalized", function()
        assert.are.equal("/etc/hosts", paths.resolve_path("/etc/hosts", 0))
    end)

    it("expands ~ to the home directory", function()
        assert.are.equal(vim.env.HOME .. "/x", paths.resolve_path("~/x", 0))
    end)

    it("expands $ENV", function()
        assert.are.equal(vim.env.HOME .. "/x", paths.resolve_path("$HOME/x", 0))
    end)

    it("anchors a relative token to the buffer's own directory", function()
        local buf = vim.api.nvim_create_buf(false, true)
        vim.api.nvim_buf_set_name(buf, "/tmp/somedir/note.md")
        assert.are.equal("/tmp/somedir/sib.lua", paths.resolve_path("./sib.lua", buf))
    end)

    it("resolves .. against the buffer's directory", function()
        local buf = vim.api.nvim_create_buf(false, true)
        vim.api.nvim_buf_set_name(buf, "/tmp/a/b/note.md")
        assert.are.equal("/tmp/a/sib.lua", paths.resolve_path("../sib.lua", buf))
    end)

    it("anchors a relative token to cwd for a nameless buffer", function()
        local buf = vim.api.nvim_create_buf(false, true)
        assert.are.equal(vim.uv.cwd() .. "/sib.lua", paths.resolve_path("./sib.lua", buf))
    end)

    it("returns nil for an empty token", function()
        assert.is_nil(paths.resolve_path("", 0))
    end)
end)

describe("utils.paths.resolve_location_under_cursor", function()
    local tmp

    --- Resolve with `lines` in the window's buffer and the cursor on the first
    --- byte of the plain substring `at`.
    local function locate(lines, at)
        local buf = vim.api.nvim_create_buf(false, true)
        vim.api.nvim_buf_set_lines(buf, 0, -1, false, lines)
        vim.api.nvim_win_set_buf(0, buf)
        for i, line in ipairs(lines) do
            local s = line:find(at, 1, true)
            if s then
                vim.api.nvim_win_set_cursor(0, { i, s - 1 })
                return paths.resolve_location_under_cursor(0)
            end
        end
        error("cursor target " .. at .. " not in the buffer")
    end

    --- /tmp is /private/tmp on macOS, so compare through resolve().
    local function assert_path(expected, actual)
        assert.are.equal(vim.fn.resolve(expected), vim.fn.resolve(actual or ""))
    end

    before_each(function()
        tmp = vim.fn.tempname()
        vim.fn.mkdir(tmp .. "/src", "p")
        for _, name in ipairs({ "x.lua", "x.md", "a,b.py", "src/x.py" }) do
            assert.are.equal(0, vim.fn.writefile({ "1", "2", "3" }, tmp .. "/" .. name))
        end
    end)

    after_each(function()
        vim.fs.rm(tmp, { recursive = true, force = true })
    end)

    it("returns the line number with the cursor on the path", function()
        local path, lnum, col = locate({ tmp .. "/x.lua:170" }, tmp)
        assert_path(tmp .. "/x.lua", path)
        assert.are.equal(170, lnum)
        assert.is_nil(col)
    end)

    it("returns the line number with the cursor on the digits", function()
        local path, lnum = locate({ tmp .. "/x.lua:170" }, "170")
        assert_path(tmp .. "/x.lua", path)
        assert.are.equal(170, lnum)
    end)

    it("returns the column of a :lnum:col suffix", function()
        local path, lnum, col = locate({ tmp .. "/x.lua:170:5" }, tmp)
        assert_path(tmp .. "/x.lua", path)
        assert.are.equal(170, lnum)
        assert.are.equal(5, col)
    end)

    it("returns no line number for a bare path", function()
        local path, lnum = locate({ tmp .. "/x.lua" }, tmp)
        assert_path(tmp .. "/x.lua", path)
        assert.is_nil(lnum)
    end)

    it("strips a #anchor", function()
        local path, lnum = locate({ tmp .. "/x.md#install" }, tmp)
        assert_path(tmp .. "/x.md", path)
        assert.is_nil(lnum)
    end)

    it("rejects a :0 suffix", function()
        local path, lnum = locate({ tmp .. "/x.lua:0" }, tmp)
        assert_path(tmp .. "/x.lua", path)
        assert.is_nil(lnum)
    end)

    it("falls back to <cfile> when the scan's token is a shorter tail", function()
        -- ',' is in 'isfname' but not in PATH_RUN, so the scan yields "b.py:12".
        local path, lnum = locate({ tmp .. "/a,b.py:12" }, "b.py")
        assert_path(tmp .. "/a,b.py", path)
        assert.is_nil(lnum)
    end)

    it("ignores a digit run too long to be a line number", function()
        assert.is_nil(locate({ "abc:12345678901234567890" }, "12345"))
    end)

    it("keeps the line number when only the column is too long", function()
        local path, lnum, col = locate({ tmp .. "/x.lua:170:12345678901234567890" }, tmp)
        assert_path(tmp .. "/x.lua", path)
        assert.are.equal(170, lnum)
        assert.is_nil(col)
    end)

    it("resolves a bare name against the buffer's own directory", function()
        local buf = vim.api.nvim_create_buf(false, true)
        vim.api.nvim_buf_set_name(buf, tmp .. "/note.md")
        vim.api.nvim_buf_set_lines(buf, 0, -1, false, { "x.lua:2" })
        vim.api.nvim_win_set_buf(0, buf)
        vim.api.nvim_win_set_cursor(0, { 1, 0 })
        local path, lnum = paths.resolve_location_under_cursor(0)
        assert_path(tmp .. "/x.lua", path)
        assert.are.equal(2, lnum)
    end)

    it("returns nil for a path-shaped token that does not exist", function()
        assert.is_nil(locate({ tmp .. "/nope.lua:12" }, tmp))
    end)

    it("ignores a time of day", function()
        assert.is_nil(locate({ "meeting at 12:30" }, "12:30"))
    end)

    it("returns the line number for a buffer-var path", function()
        local path, lnum = locate({ "PEEKROOT=" .. tmp, "$PEEKROOT/src/x.py:170" }, "$PEEKROOT")
        assert_path(tmp .. "/src/x.py", path)
        assert.are.equal(170, lnum)
    end)
end)

describe("gf mapping", function()
    -- The runner starts nvim with --noplugin (plenary test_harness.lua:90), so
    -- the mapping is not installed until this file is sourced. Sourcing also
    -- prepends the repo and all of $PATH to 'path', which findfile() reads:
    -- restore it, or the resolver tests resolve against a different 'path'
    -- than they were written for.
    local saved_path = vim.o.path
    vim.cmd.runtime("plugin/paths.lua")
    vim.o.path = saved_path

    local dir
    local FIXTURE_LINES = 5

    --- Numbered lines, `first` in place of line 1 where given.
    local function fixture_lines(first)
        local lines = { first or "line 1" }
        for i = 2, FIXTURE_LINES do
            lines[i] = "line " .. i
        end
        return lines
    end

    --- Write a fixture file `name`, and leave the window on a buffer holding
    --- that path plus `suffix`, cursor on its first byte. A fixture per case:
    --- `:edit` restores the position it remembers for an already-loaded file,
    --- so a reused one starts wherever the previous case left it.
    local function gf_from(suffix, name)
        local path = dir .. "/" .. name
        assert.are.equal(0, vim.fn.writefile(fixture_lines(), path))
        vim.cmd("enew!") -- bang: a failed case leaves a modified buffer behind
        vim.api.nvim_buf_set_lines(0, 0, -1, false, { path .. suffix })
        vim.api.nvim_win_set_cursor(0, { 1, 0 })
        return path
    end

    --- Buffer name (symlinks resolved, for /tmp vs /private/tmp) and cursor
    --- after the mapping. `normal` without the bang, or it is bypassed.
    local function gf()
        vim.cmd("normal gf")
        local r, c = unpack(vim.api.nvim_win_get_cursor(0))
        return vim.fn.resolve(vim.api.nvim_buf_get_name(0)), r, c
    end

    before_each(function()
        dir = vim.fn.tempname()
        vim.fn.mkdir(dir, "p")
    end)

    after_each(function()
        vim.fs.rm(dir, { recursive = true, force = true })
    end)

    it("lands on the :lnum line", function()
        local path = gf_from(":3", "a.lua")
        local name, r, c = gf()
        assert.are.equal(vim.fn.resolve(path), name)
        assert.are.same({ 3, 0 }, { r, c })
    end)

    it("lands on the :lnum:col column", function()
        gf_from(":3:5", "b.lua")
        local _, r, c = gf()
        assert.are.same({ 3, 4 }, { r, c })
    end)

    it("clamps a line number past the end of the file", function()
        gf_from(":9999", "c.lua")
        local _, r = gf()
        assert.are.equal(FIXTURE_LINES, r)
    end)

    it("opens a name holding # and %", function()
        -- Both are Ex-command metacharacters (alternate and current file), so
        -- the name has to reach the buffer without going through one. The
        -- :lnum keeps this on the resolver branch; a bare name would fall
        -- through to built-in gf and prove nothing about ours.
        local path = gf_from(":3", "a#b%c.lua")
        local name, r = gf()
        assert.are.equal(vim.fn.resolve(path), name)
        assert.are.equal(3, r)
    end)

    it("opens at line 1 when the suffix is line 0", function()
        local path = gf_from(":0", "d.lua")
        local name, r = gf()
        assert.are.equal(vim.fn.resolve(path), name)
        assert.are.equal(1, r)
    end)

    it("opens a bare path at line 1", function()
        local path = gf_from("", "e.lua")
        local name, r = gf()
        assert.are.equal(vim.fn.resolve(path), name)
        assert.are.equal(1, r)
    end)

    it("jumps within the current buffer without reloading it", function()
        -- The fixture points at itself, and `dir` is under a symlinked /tmp or
        -- /var, so the path never matches the buffer name literally. Opening
        -- the file a second time would lose the edit below, and <C-o> only
        -- comes back from the same buffer if the departure was recorded.
        local path = dir .. "/self.lua"
        assert.are.equal(0, vim.fn.writefile(fixture_lines(path .. ":4"), path))
        vim.cmd.edit(path)
        vim.api.nvim_buf_set_lines(0, 1, 2, false, { "edited" })
        vim.api.nvim_win_set_cursor(0, { 1, 0 })

        local name, r = gf()
        assert.are.equal(vim.fn.resolve(path), name)
        assert.are.equal(4, r)
        assert.is_true(vim.bo.modified) -- not reloaded from disk

        vim.cmd("normal! \15") -- <C-o>
        assert.are.equal(vim.fn.resolve(path), vim.fn.resolve(vim.api.nvim_buf_get_name(0)))
        assert.are.equal(1, vim.api.nvim_win_get_cursor(0)[1])
    end)
end)
