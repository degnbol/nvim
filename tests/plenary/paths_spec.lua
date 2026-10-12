---@diagnostic disable: undefined-global
-- Tests for lua/utils/paths.lua: resolve_path (pure single-candidate
-- resolution), resolve_location_under_cursor (the location resolver, with its
-- :lnum[:col] suffix parsing) and the 'includeexpr' name transforms, plus the
-- gf-family mappings and FileType wrapper in plugin/paths.lua.

local paths = require("utils.paths")

-- The runner starts nvim with --noplugin (plenary test_harness.lua:90), so the
-- mappings and autocmds are not installed until this file is sourced. Sourcing
-- also adds the repo and all of $PATH to 'path', which findfile() reads:
-- restore it, or the resolver tests resolve against a different 'path' than
-- they were written for.
local saved_path = vim.o.path
vim.cmd.runtime("plugin/paths.lua")
vim.o.path = saved_path

--- /tmp is /private/tmp on macOS, so compare through resolve().
local function assert_path(expected, actual)
    assert.are.equal(vim.fn.resolve(expected), vim.fn.resolve(actual or ""))
end

--- A temp dir holding .git/, src/, and five-line x.lua and x.md.
local function make_repo()
    local dir = vim.fn.tempname()
    vim.fn.mkdir(dir .. "/.git", "p")
    vim.fn.mkdir(dir .. "/src", "p")
    for _, name in ipairs({ "x.lua", "x.md" }) do
        assert.are.equal(0, vim.fn.writefile({ "line 1", "line 2", "line 3", "line 4", "line 5" },
            dir .. "/" .. name))
    end
    return dir
end

--- A fresh buffer in the current window, named `name`, holding `lines`.
local function show(name, lines)
    vim.cmd("enew!")
    vim.api.nvim_buf_set_name(0, name)
    vim.api.nvim_buf_set_lines(0, 0, -1, false, lines)
    return vim.api.nvim_get_current_buf()
end

describe("utils.paths.expand_dollars", function()
    local repo

    before_each(function() repo = make_repo() end)
    after_each(function() vim.fs.rm(repo, { recursive = true, force = true }) end)

    it("expands $(git root)", function()
        local buf = show(repo .. "/note.md", {})
        assert_path(repo .. "/x", paths.expand_dollars("$(git root)/x", buf))
    end)

    it("expands an environment variable", function()
        vim.env.PATHS_SPEC_VAR = "/a"
        assert.are.equal("/a/x", paths.expand_dollars("$PATHS_SPEC_VAR/x", 0))
        vim.env.PATHS_SPEC_VAR = nil
    end)

    it("expands a buffer NAME=value assignment", function()
        local buf = show(repo .. "/note.md", { "PEEK='/p'" })
        assert.are.equal("/p/x", paths.expand_dollars("$PEEK/x", buf))
    end)

    it("keeps unresolvable parts literal", function()
        local buf = show(repo .. "/note.md", {})
        assert.are.equal("$(ls)/$NOPE_X/y", paths.expand_dollars("$(ls)/$NOPE_X/y", buf))
    end)
end)

describe("utils.paths.transform_name", function()
    local repo

    before_each(function() repo = make_repo() end)
    after_each(function() vim.fs.rm(repo, { recursive = true, force = true }) end)

    it("removes a trailing #anchor", function()
        local buf = show(repo .. "/note.md", {})
        assert.are.equal(repo .. "/x.md", paths.transform_name(repo .. "/x.md#install", buf))
    end)

    it("expands a buffer $NAME", function()
        local buf = show(repo .. "/note.md", { "R=" .. repo })
        assert.are.equal(repo .. "/x.md", paths.transform_name("$R/x.md", buf))
    end)

    it("returns a directory", function()
        local buf = show(repo .. "/note.md", { "R=" .. repo })
        assert.are.equal(repo .. "/src", paths.transform_name("$R/src", buf))
    end)

    it("returns nil when no candidate exists", function()
        local buf = show(repo .. "/note.md", {})
        assert.is_nil(paths.transform_name(repo .. "/nope.md#a", buf))
    end)
end)

describe("utils.paths.includeexpr", function()
    local repo

    before_each(function() repo = make_repo() end)
    after_each(function() vim.fs.rm(repo, { recursive = true, force = true }) end)

    --- A lua buffer, so the runtime ftplugin sets its own 'includeexpr'.
    local function lua_buffer()
        local buf = show(repo .. "/init.lua", {})
        vim.bo.filetype = "lua"
        return buf
    end

    --- The name of the buffer gf opens from the first byte of `line` in a lua
    --- buffer. Through gf, as the filetype's expression reads v:fname, which
    --- is set only while 'includeexpr' is evaluated.
    local function gf_in_lua(line)
        lua_buffer()
        vim.api.nvim_buf_set_lines(0, 0, -1, false, { line })
        vim.api.nvim_win_set_cursor(0, { 1, 0 })
        vim.cmd("normal gf")
        return vim.api.nvim_buf_get_name(0)
    end

    it("moves the filetype's expression to b:ft_includeexpr", function()
        lua_buffer()
        assert.are.equal(paths.INCLUDEEXPR, vim.bo.includeexpr)
        assert.is_not_nil(vim.b.ft_includeexpr)
    end)

    it("resolves through the filetype's expression", function()
        assert.is_truthy(gf_in_lua("utils.paths"):match("/lua/utils/paths%.lua$"))
    end)

    it("resolves through the global transforms after the filetype's expression", function()
        -- Absolute: the runtime lua ftplugin removes "." from 'path'.
        assert_path(repo .. "/x.md", gf_in_lua(repo .. "/x.md#a"))
    end)

    it("keeps the filetype's expression when FileType fires again", function()
        lua_buffer()
        local ft_expr = vim.b.ft_includeexpr
        vim.api.nvim_exec_autocmds("FileType", { pattern = "lua", modeline = false })
        assert.are.equal(ft_expr, vim.b.ft_includeexpr)
        assert.are.equal(paths.INCLUDEEXPR, vim.bo.includeexpr)
    end)

    it("drops the filetype's expression when the filetype changes", function()
        lua_buffer()
        vim.cmd("setlocal filetype=text")
        assert.is_nil(vim.b.ft_includeexpr)
        assert.are.equal(paths.INCLUDEEXPR, vim.bo.includeexpr)
    end)

    it("returns the name unchanged when nothing resolves", function()
        show(repo .. "/note.md", {})
        assert.are.equal("nope.md#a", paths.includeexpr("nope.md#a"))
    end)
end)

describe("gf family mappings", function()
    local repo

    before_each(function()
        repo = make_repo()
        vim.cmd("silent! tabonly | silent! only")
    end)
    after_each(function() vim.fs.rm(repo, { recursive = true, force = true }) end)

    --- Leave a buffer named note.md in the repo holding `line`, cursor on the
    --- first byte of the plain substring `at`.
    local function from(line, at)
        show(repo .. "/note.md", { line })
        vim.api.nvim_win_set_cursor(0, { 1, line:find(at, 1, true) - 1 })
    end

    local function cursor() return vim.api.nvim_win_get_cursor(0) end

    it("<C-w>f opens $(git root)/x:lnum:col in a split", function()
        from("$(git root)/x.lua:3:2", "x.lua")
        vim.cmd("normal \23f")
        assert.are.equal(2, #vim.api.nvim_tabpage_list_wins(0))
        assert_path(repo .. "/x.lua", vim.api.nvim_buf_get_name(0))
        assert.are.same({ 3, 1 }, cursor())
    end)

    it("<C-w>gf opens $(git root)/x:lnum:col in a tab", function()
        from("$(git root)/x.lua:3:2", "x.lua")
        vim.cmd("normal \23gf")
        assert.are.equal(2, #vim.api.nvim_list_tabpages())
        assert_path(repo .. "/x.lua", vim.api.nvim_buf_get_name(0))
        assert.are.same({ 3, 1 }, cursor())
    end)

    it("gf on `git` in $(git root)/x opens the git-root file", function()
        -- vim.go: the current buffer may still be the previous case's x.lua,
        -- whose ftplugin sets a local 'path' that vim.o would read and write.
        local saved = vim.go.path
        vim.go.path = saved .. ",/usr/bin"
        from("$(git root)/x.lua", "git")
        vim.cmd("normal gf")
        vim.go.path = saved
        assert_path(repo .. "/x.lua", vim.api.nvim_buf_get_name(0))
    end)

    it("gf prefers the buffer's directory over an earlier 'path' entry", function()
        vim.fn.mkdir(repo .. "/sub", "p")
        assert.are.equal(0, vim.fn.writefile({ "sub" }, repo .. "/sub/x.md"))
        local saved = vim.go.path
        vim.go.path = repo .. ",.,,"
        for i, name in ipairs({ "x.md", "x.md#a" }) do
            show(repo .. "/sub/note" .. i .. ".md", { name })
            vim.api.nvim_win_set_cursor(0, { 1, 0 })
            vim.cmd("normal gf")
            assert_path(repo .. "/sub/x.md", vim.api.nvim_buf_get_name(0))
        end
        vim.go.path = saved
    end)

    it("never runs a backtick command from a buffer NAME=value", function()
        local marker = repo .. "/ran"
        from("X=`touch " .. marker .. "`", "X")
        vim.api.nvim_buf_set_lines(0, -1, -1, false, { "$X/foo" })
        vim.api.nvim_win_set_cursor(0, { 2, 0 })
        paths.resolve_location_under_cursor()
        pcall(vim.api.nvim_command, "normal gf")
        assert.is_nil(vim.uv.fs_stat(marker))
    end)

    it("<C-w>f on a miss reports E447 as the native key does", function()
        from("nope.md", "nope")
        -- An error message inside :normal arrives here as an exception.
        local ok, err = pcall(vim.api.nvim_command, "normal \23f")
        assert.is_false(ok)
        assert.is_truthy(tostring(err):match("E447:"))
        assert.is_nil(tostring(err):match("E5108")) -- no Lua traceback
        assert.are.equal(1, #vim.api.nvim_tabpage_list_wins(0))
    end)

    it("gF on x.md#a:4 opens line 4 through includeexpr", function()
        from("x.md#a:4", "x")
        vim.cmd("normal gF")
        assert_path(repo .. "/x.md", vim.api.nvim_buf_get_name(0))
        assert.are.equal(4, cursor()[1])
    end)
end)

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

    it("anchors a relative token to cwd for a URL-named buffer", function()
        local buf = vim.api.nvim_create_buf(false, true)
        vim.api.nvim_buf_set_name(buf, "agentic://1/input")
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
                return paths.resolve_location_under_cursor()
            end
        end
        error("cursor target " .. at .. " not in the buffer")
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
        local path, lnum = paths.resolve_location_under_cursor()
        assert_path(tmp .. "/x.lua", path)
        assert.are.equal(2, lnum)
    end)

    it("follows a symlink before .. in a path relative to the buffer", function()
        vim.fn.mkdir(tmp .. "/real/sub", "p")
        vim.fn.mkdir(tmp .. "/base", "p")
        assert.are.equal(0, vim.fn.writefile({ "1" }, tmp .. "/real/g.txt"))
        assert(vim.uv.fs_symlink(tmp .. "/real/sub", tmp .. "/base/lnk"))
        local buf = vim.api.nvim_create_buf(false, true)
        vim.api.nvim_buf_set_name(buf, tmp .. "/base/note.md")
        vim.api.nvim_buf_set_lines(buf, 0, -1, false, { "lnk/../g.txt" })
        vim.api.nvim_win_set_buf(0, buf)
        vim.api.nvim_win_set_cursor(0, { 1, 0 })
        assert_path(tmp .. "/real/g.txt", paths.resolve_location_under_cursor())
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
