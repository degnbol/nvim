---@diagnostic disable: undefined-global
-- Tests for lua/completion/cmp_zsh_capture.zsh: the words it lists for a line.

local shared = vim.fn.stdpath("data") .. "/site/pack/core/opt/cmp-zsh/bin/cmp_capture_shared.zsh"
if vim.fn.filereadable(shared) == 0 then
    pending("cmp-zsh is not installed at " .. shared)
    return
end
local script = vim.fn.getcwd() .. "/lua/completion/cmp_zsh_capture.zsh"

--- The distinct words the capture lists for `line`, run in `dir`, sorted.
local function capture(line, dir)
    local result = vim.system({ "zsh", script, line }, { cwd = dir, env = { CMP_ZSH_SHARED = shared } }):wait()
    assert.are.equal(0, result.code, result.stderr)
    local words = vim.tbl_map(function(l) return (l:gsub("\r$", ""):gsub(" %-%- .*", "")) end,
        vim.split(vim.trim(result.stdout), "\n"))
    table.sort(vim.list.unique(words))
    return words
end

describe("cmp_zsh_capture.zsh", function()
    local dir

    before_each(function()
        dir = vim.fn.tempname()
        vim.fn.mkdir(dir .. "/src", "p")
        for _, name in ipairs({ "x.lua", "y.l", "src/align.py" }) do
            assert.are.equal(0, vim.fn.writefile({}, dir .. "/" .. name))
        end
    end)
    after_each(function() vim.fs.rm(dir, { recursive = true, force = true }) end)

    it("lists the matches of a glob", function()
        assert.are.same({ "x.lua", "y.l" }, capture("ls *.l", dir))
    end)

    it("lists the contents of a $ directory, not its expansion", function()
        assert.are.same({ "$PWD/src/align.py" }, capture("ls $PWD/src/", dir))
    end)
end)
