---@diagnostic disable: undefined-global
-- Tests for lua/textobjects.lua through mini.ai, and the operator-pending `r`
-- map in plugin/keymaps.lua.

vim.cmd.packadd "mini.nvim"
vim.cmd.packadd "nvim-treesitter-textobjects"
-- Its plugin/filetypes.lua gives filetype tex the latex parser.
vim.cmd.packadd "nvim-treesitter"
-- ftplugin/tex.lua calls vimtex functions. vimtex reads its settings, which
-- keep it off `ic`/`ac`, when it loads.
local vimtex_spec = vim.iter(require "plugins/tex"):find(function(spec) return spec[1] == "vimtex" end)
assert(vimtex_spec, "no vimtex spec in lua/plugins/tex.lua")
vimtex_spec.before()
vim.cmd.packadd "vimtex"
-- Language options come from vim.filetype.get_option, which sources ftplugins
-- only with `:filetype plugin on`. It caches the result per filetype, so the
-- setting must not change between tests.
vim.cmd "filetype plugin on"
require("mini.ai").setup {
    custom_textobjects = vim.tbl_map(function(o) return o.spec end, require "textobjects"),
}

--- Open a scratch buffer of `lines` with the cursor on the first `target`.
--- @param lines string[]
--- @param filetype string "" for none. Its parser is started when one exists.
--- @param commentstring string
--- @param target string text to place the cursor on
local function open(lines, filetype, commentstring, target)
    local buf = vim.api.nvim_create_buf(false, true)
    vim.api.nvim_set_current_buf(buf)
    vim.api.nvim_buf_set_lines(buf, 0, -1, false, lines)
    vim.bo[buf].filetype = filetype
    vim.bo[buf].commentstring = commentstring
    local parser = vim.treesitter.get_parser(buf, nil, { error = false })
    if parser then parser:parse(true) end
    for row, line in ipairs(lines) do
        local col = line:find(target, 1, true)
        if col then return vim.api.nvim_win_set_cursor(0, { row, col - 1 }) end
    end
    error("no " .. target)
end

--- Yank `motion` in the current buffer, then wipe the buffer.
--- @param motion string e.g. "ic"
--- @return string yanked
local function yank(motion)
    vim.fn.setreg('"', "")
    vim.cmd.normal("y" .. motion)
    vim.api.nvim_buf_delete(0, { force = true })
    return vim.fn.getreg('"')
end

describe("comment textobject", function()
    it("python: inner drops '# ', outer keeps it", function()
        open({ "x = 1", "# hello there", "y = 2" }, "python", "# %s", "there")
        assert.are.equal("hello there", yank "ic")
        open({ "x = 1", "# hello there", "y = 2" }, "python", "# %s", "there")
        assert.are.equal("# hello there", yank "ac")
    end)
    it("python: cursor on the leader", function()
        open({ "# hello" }, "python", "# %s", "#")
        assert.are.equal("hello", yank "ic")
    end)
    it("trailing comment after code", function()
        open({ "x = 1  # note" }, "python", "# %s", "note")
        assert.are.equal("note", yank "ic")
    end)
    it("julia line comment", function()
        open({ "x = 1", "# hello" }, "julia", "# %s", "hello")
        assert.are.equal("hello", yank "ic")
    end)
    it("latex falls back to full-line comments", function()
        open({ "text", "% a remark", "more" }, "tex", "% %s", "remark")
        assert.are.equal("a remark", yank "ic")
        open({ "text", "% a remark", "more" }, "tex", "% %s", "remark")
        assert.are.equal("% a remark", yank "ac")
    end)
    it("C /* */ drops the 'comments' pair that 'commentstring' lacks", function()
        open({ "int x; /* note */" }, "c", "// %s", "note")
        assert.are.equal("note", yank "ic")
        open({ "/* one", "   two */", "int x;" }, "c", "// %s", "two")
        assert.are.equal("one\n   two", yank "ic")
    end)
    it("drops the longest leader in 'comments'", function()
        open({ "--- doc", "local x" }, "lua", "-- %s", "doc")
        assert.are.equal("doc", yank "ic")
        open({ "/// doc", "fn f() {}" }, "rust", "// %s", "doc")
        assert.are.equal("doc", yank "ic")
        open({ "//! doc", "fn f() {}" }, "rust", "// %s", "doc")
        assert.are.equal("doc", yank "ic")
    end)
    it("lua block comments keep only the query's @comment.inner", function()
        open({ "--[[ x ]]", "local y" }, "lua", "-- %s", "x")
        assert.are.equal("x", yank "ic")
        open({ "--[==[", "  y", "]==]" }, "lua", "-- %s", "y")
        assert.are.equal("y", yank "ic")
    end)
    it("inner of the comment under the cursor when outside its inner part", function()
        open({ "--[[ x ]]", "--[[ z ]]" }, "lua", "-- %s", "]]")
        assert.are.equal("x", yank "ic")
        open({ "# a  ", "# b" }, "python", "# %s", "  ")
        assert.are.equal("a", yank "ic")
    end)
    it("inner of an empty comment is empty after its delimiters", function()
        open({ "--[[]]" }, "lua", "-- %s", "[[")
        vim.cmd.normal "cicX"
        assert.are.same({ "--[[X]]" }, vim.api.nvim_buf_get_lines(0, 0, -1, true))
        vim.api.nvim_buf_delete(0, { force = true })
    end)
    it("selects nothing without a comment syntax", function()
        open({ "# not a comment" }, "", "", "not")
        assert.are.equal("", yank "ic")
    end)
    it("python fence in markdown uses python's comment", function()
        open({ "text", "```python", "x = 1  # inside", "```" }, "markdown", "<!-- %s -->", "inside")
        assert.are.equal("inside", yank "ic")
    end)
end)

describe("class textobject", function()
    it("python class", function()
        open({ "class A:", "    x = 1", "y = 2" }, "python", "# %s", "x")
        assert.are.equal("x = 1", yank "iC")
        open({ "class A:", "    x = 1", "y = 2" }, "python", "# %s", "x")
        assert.are.equal("class A:\n    x = 1", yank "aC")
    end)
end)

describe("tex", function()
    it("f is a command, its body the argument under the cursor", function()
        open({ "a \\textbf{bold} b" }, "tex", "% %s", "bold")
        assert.are.equal("bold", yank "if")
        open({ "a \\textbf{bold} b" }, "tex", "% %s", "bold")
        assert.are.equal("\\textbf{bold}", yank "af")
        open({ "a \\section[short]{Long} b" }, "tex", "% %s", "Long")
        assert.are.equal("Long", yank "if")
        open({ "a \\section[short]{Long} b" }, "tex", "% %s", "Long")
        assert.are.equal("\\section[short]{Long}", yank "af")
    end)
    it("F is a command name", function()
        open({ "a $\\alpha_{x}$ b" }, "tex", "% %s", "alpha")
        assert.are.equal("alpha", yank "iF")
        open({ "a \\section*{S} b" }, "tex", "% %s", "section")
        assert.are.equal("\\section*", yank "aF")
        open({ "\\makeatletter \\@ifstar\\foo" }, "tex", "% %s", "ifstar")
        assert.are.equal("@ifstar", yank "iF")
    end)
    it("F is the name of the command around the cursor", function()
        open({ "x \\textbf{\\emph{y}} z" }, "tex", "% %s", "y")
        assert.are.equal("\\emph", yank "aF")
    end)
end)

describe("rest of paragraph", function()
    vim.cmd.runtime "plugin/keymaps.lua"
    local para = { "a", "b", "c", "", "d", "e" }

    it("from the middle of a paragraph", function()
        open(para, "", "", "b")
        assert.are.equal("b\nc\n", yank "r")
    end)
    it("from a blank line takes it and the next paragraph", function()
        open(para, "", "", "")
        vim.api.nvim_win_set_cursor(0, { 4, 0 })
        assert.are.equal("\nd\ne\n", yank "r")
    end)
    it("on a blank last line", function()
        open({ "a", "" }, "", "", "a")
        vim.api.nvim_win_set_cursor(0, { 2, 0 })
        assert.are.equal("\n", yank "r")
    end)
    it("last paragraph without trailing blank line", function()
        open(para, "", "", "d")
        assert.are.equal("d\ne\n", yank "r")
    end)
end)
