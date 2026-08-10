---@diagnostic disable: undefined-global
-- Tests for lua/path_highlight.lua (anchored-path detection + prose gate +
-- existence resolution) and utils.paths.resolve_path (pure single-candidate
-- resolution). The decoration provider's extmark placement is not exercised
-- here — headless has no redraw to drive it — but every piece it composes is.

local paths = require("utils.paths")
local ph = require("path_highlight")

--- Fresh scratch buffer with the given lines, filetype markdown, parser primed
--- (injections parsed so the prose gate can see code_span/link_destination).
local function md_buf(src)
    local buf = vim.api.nvim_create_buf(false, true)
    vim.api.nvim_buf_set_lines(buf, 0, -1, false, vim.split(src, "\n", { plain = true }))
    vim.bo[buf].filetype = "markdown"
    vim.treesitter.get_parser(buf, "markdown"):parse(true)
    return buf
end

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

describe("path_highlight.scan", function()
    local function tokens(line)
        local out = {}
        for _, m in ipairs(ph.scan(line)) do
            out[#out + 1] = m.token
        end
        return out
    end

    it("detects each anchor form", function()
        assert.are.same({ "./foo/bar.lua" }, tokens("see ./foo/bar.lua here"))
        assert.are.same({ "../up.md" }, tokens("go ../up.md"))
        assert.are.same({ "/abs/path" }, tokens("at /abs/path end"))
        assert.are.same({ "~/y.txt" }, tokens("open ~/y.txt now"))
        assert.are.same({ "$HOME/z" }, tokens("under $HOME/z ok"))
    end)

    it("reports the 0-indexed byte column", function()
        local m = ph.scan("see ./x")
        assert.are.equal(1, #m)
        assert.are.equal(4, m[1].col)
    end)

    it("rejects an unanchored slash token", function()
        assert.are.same({}, tokens("this and/or that"))
        assert.are.same({}, tokens("a bare foo/bar token"))
    end)

    it("rejects a lone slash or bare anchor", function()
        assert.are.same({}, tokens("a / b"))
        assert.are.same({}, tokens("trailing ./ only"))
        assert.are.same({}, tokens("home ~/ alone"))
    end)

    it("finds multiple candidates on one line", function()
        assert.are.same({ "./a", "/b/c" }, tokens("./a and /b/c"))
    end)
end)

describe("path_highlight.in_prose", function()
    it("accepts a path in plain prose", function()
        local buf = md_buf("see ./foo.lua now")
        assert.is_true(ph.in_prose(buf, 0, 4))
    end)

    it("rejects a path inside an inline code span", function()
        local buf = md_buf("see `./foo.lua` now")
        assert.is_false(ph.in_prose(buf, 0, 5))
    end)

    it("rejects a path inside a link destination", function()
        local buf = md_buf("see [x](./foo.lua) now")
        assert.is_false(ph.in_prose(buf, 0, 8))
    end)

    it("rejects a path inside a known-language fenced code block", function()
        local buf = md_buf("```lua\nlocal p = ./foo.lua\n```")
        assert.is_false(ph.in_prose(buf, 1, 10))
    end)

    it("rejects a path inside a bare / unknown-language fenced block", function()
        -- No injected language, so language_for_range stays markdown; the
        -- ancestor climb must catch fenced_code_block.
        assert.is_false(ph.in_prose(md_buf("```\n./foo.lua\n```"), 1, 0))
        assert.is_false(ph.in_prose(md_buf("```text\n./foo.lua\n```"), 1, 0))
    end)

    it("rejects a path inside an indented code block", function()
        assert.is_false(ph.in_prose(md_buf("    ./foo.lua"), 0, 4))
    end)
end)

describe("path_highlight.existing_len", function()
    local tmp, existing

    before_each(function()
        tmp = vim.fn.tempname()
        vim.fn.mkdir(tmp, "p")
        existing = tmp .. "/real.md"
        assert.are.equal(0, vim.fn.writefile({ "x" }, existing))
    end)

    after_each(function()
        vim.fs.rm(tmp, { recursive = true, force = true })
    end)

    it("returns the token length for an existing absolute path", function()
        assert.are.equal(#existing, ph.existing_len(existing, 0))
    end)

    it("returns nil for a non-existent path", function()
        assert.is_nil(ph.existing_len(tmp .. "/nope.md", 0))
    end)

    it("strips trailing sentence punctuation to resolve", function()
        -- "…/real.md." — raw fails, stripped resolves; length excludes the dot.
        assert.are.equal(#existing, ph.existing_len(existing .. ".", 0))
    end)

    it("prefers the raw token when it exists", function()
        local dir = tmp .. "/.." -- ends in dots but is a real directory
        assert.are.equal(#dir, ph.existing_len(dir, 0))
    end)
end)
