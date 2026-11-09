local tags = require "symbol-ls.tags"

local M = {}

--- Call `fn` with a callback and wait for it. Fails after 3 s.
--- @param fn fun(done: fun(...))
--- @return any ... the callback's arguments
function M.await(fn)
    local args
    fn(function(...) args = { ... } end)
    assert(vim.wait(3000, function() return args ~= nil end), "callback not called")
    return unpack(args)
end

local tags_dir = tags.dir

--- Point `tags.dir` at a new directory on 'runtimepath' holding
--- `tests/fixtures/lua_tags.scm` as the lua tags query, current for `revision`.
--- Undo with `reset_tags_dir`.
--- @param revision string
--- @return string dir
function M.fixture_tags_dir(revision)
    local dir = vim.fn.tempname()
    vim.fn.mkdir(dir .. "/queries/lua", "p")
    vim.fn.writefile(vim.fn.readfile("tests/fixtures/lua_tags.scm"), dir .. "/queries/lua/tags.scm")
    vim.fn.writefile({ vim.json.encode { lua = { revision = revision, tags = true } } }, dir .. "/revisions.json")
    tags.dir = function() return dir end
    vim.opt.runtimepath:append(dir)
    return dir
end

--- Restore `tags.dir` and take the fixture directory off 'runtimepath'.
function M.reset_tags_dir()
    vim.opt.runtimepath:remove(tags.dir())
    tags.dir = tags_dir
end

--- New git repository with `files` (path relative to the repo → content)
--- tracked and `untracked` written but not added.
--- @param files table<string, string>
--- @param untracked table<string, string>|nil
--- @return string root resolved, as getcwd() gives it after a chdir there
function M.git_repo(files, untracked)
    local root = vim.fn.tempname()
    vim.fn.mkdir(root, "p")
    root = assert(vim.uv.fs_realpath(root))
    for path, content in pairs(vim.tbl_extend("force", files, untracked or {})) do
        vim.fn.mkdir(vim.fs.dirname(root .. "/" .. path), "p")
        vim.fn.writefile(vim.split(content, "\n"), root .. "/" .. path)
    end
    assert(vim.system({ "git", "init", "-q" }, { cwd = root }):wait().code == 0)
    local add = vim.list_extend({ "git", "add" }, vim.tbl_keys(files))
    assert(vim.system(add, { cwd = root }):wait().code == 0)
    return root
end

return M
