---@diagnostic disable: undefined-global
-- Tests for lua/completion/path.lua: blink's path dirname with $ forms expanded.

local paths = require("utils.paths")

local blink_dir = vim.fn.stdpath("data") .. "/site/pack/core/opt/blink.cmp"
if not paths.is_dir(blink_dir) then
    pending("blink.cmp is not installed at " .. blink_dir)
    return
end
vim.opt.runtimepath:append(blink_dir)

local cp = require("completion.path")

--- /tmp is /private/tmp on macOS, so compare through resolve().
local function assert_path(expected, actual)
    assert.are.equal(vim.fn.resolve(expected), actual)
end

--- A fresh buffer named `name`, holding `lines` and then `line`, with a blink
--- context for the cursor at the end of `line`.
local function context_at_end(name, lines, line)
    vim.cmd("enew!")
    vim.api.nvim_buf_set_name(0, name)
    vim.api.nvim_buf_set_lines(0, 0, -1, false, vim.list_extend(vim.list_extend({}, lines), { line }))
    local start_col = line:find("[%w_%-]*$")
    return {
        bufnr = vim.api.nvim_get_current_buf(),
        line = line,
        bounds = { start_col = start_col, length = #line - start_col + 1 },
    }
end

describe("completion.path.expanded_dirname", function()
    local repo, saved_root

    before_each(function()
        repo = vim.fn.tempname()
        vim.fn.mkdir(repo .. "/.git", "p")
        vim.fn.mkdir(repo .. "/src", "p")
        saved_root = vim.env.ROOT
        vim.env.ROOT = repo
    end)
    after_each(function()
        vim.env.ROOT = saved_root
        vim.env.PCS_DIR = nil
        vim.fs.rm(repo, { recursive = true, force = true })
    end)

    it("expands $(git root) with the cursor right after /", function()
        local ctx = context_at_end(repo .. "/x.sh", {}, "ls $(git root)/src/")
        assert_path(repo .. "/src", cp.expanded_dirname(cp.opts, ctx))
    end)

    it("expands an environment variable", function()
        local ctx = context_at_end(vim.fn.tempname() .. ".sh", {}, "ls $ROOT/src/")
        assert_path(repo .. "/src", cp.expanded_dirname(cp.opts, ctx))
    end)

    it("prefers a buffer NAME=value line over the environment", function()
        vim.env.PCS_DIR = "/nonexistent"
        local ctx = context_at_end(vim.fn.tempname() .. ".sh", { "PCS_DIR=" .. repo }, "ls $PCS_DIR/src/")
        assert_path(repo .. "/src", cp.expanded_dirname(cp.opts, ctx))
    end)

    it("collapses the // from a value that ends in /", function()
        vim.env.PCS_DIR = repo .. "/"
        local ctx = context_at_end(vim.fn.tempname() .. ".sh", {}, "ls $PCS_DIR/src/")
        assert_path(repo .. "/src", cp.expanded_dirname(cp.opts, ctx))
    end)

    it("shifts the keyword start with the expansion", function()
        local ctx = context_at_end(vim.fn.tempname() .. ".sh", {}, "ls $ROOT/src/al")
        assert_path(repo .. "/src", cp.expanded_dirname(cp.opts, ctx))
    end)

    it("returns nil for a directory that does not exist", function()
        local ctx = context_at_end(vim.fn.tempname() .. ".sh", {}, "ls $ROOT/nope/")
        assert.is_nil(cp.expanded_dirname(cp.opts, ctx))
    end)
end)
