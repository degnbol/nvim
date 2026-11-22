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

--- Run normal-mode `keys` in the current buffer, then wipe the buffer.
--- @param keys string e.g. "yic"
--- @return string register `"` after `keys`, "" when they yank nothing
local function register_after(keys)
    vim.fn.setreg('"', "")
    vim.cmd.normal(keys)
    vim.api.nvim_buf_delete(0, { force = true })
    return vim.fn.getreg('"')
end

--- Yank `motion` in the current buffer, then wipe the buffer.
--- @param motion string e.g. "ic"
--- @return string yanked
local function yank(motion) return register_after("y" .. motion) end

--- Run normal-mode `keys` in the current buffer, then wipe the buffer.
--- @param keys string e.g. "dak"
--- @return string[] lines of the buffer after `keys`
local function lines_after(keys)
    vim.cmd.normal(keys)
    local lines = vim.api.nvim_buf_get_lines(0, 0, -1, true)
    vim.api.nvim_buf_delete(0, { force = true })
    return lines
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

describe("code fence textobject", function()
    local markdown = { "a", "", "```py", "x = 1", "y = 2", "```", "", "b" }
    it("markdown: inner is the content, outer the fenced block", function()
        open(markdown, "markdown", "<!-- %s -->", "x =")
        assert.are.equal("x = 1\ny = 2\n", yank "iF")
        open(markdown, "markdown", "<!-- %s -->", "x =")
        assert.are.equal("```py\nx = 1\ny = 2\n```\n", yank "aF")
    end)
    it("julia docstring", function()
        local julia = { '"""', "Doc.", "", "```julia", "f(1)", "```", '"""', "f(x) = x" }
        open(julia, "julia", "# %s", "f(1)")
        assert.are.equal("f(1)\n", yank "iF")
        open(julia, "julia", "# %s", "f(1)")
        assert.are.equal("```julia\nf(1)\n```\n", yank "aF")
    end)
end)

describe("key/value textobject", function()
    -- Keys have more than one character: in Visual mode the reference region
    -- of a one-character key equals the key.
    local nested = { "f(aa=g(bb=11), cc=22)" }

    it("nested pairs: the pair around the cursor", function()
        open(nested, "julia", "# %s", "bb")
        assert.are.equal("11", yank "iv")
        open(nested, "julia", "# %s", "11")
        assert.are.equal("bb", yank "ik")
    end)
    it("a count selects an enclosing pair", function()
        open(nested, "julia", "# %s", "bb")
        assert.are.equal("g(bb=11)", yank "2iv")
    end)
    it("repeating in Visual mode steps out", function()
        open(nested, "julia", "# %s", "bb")
        assert.are.equal("aa", register_after "vikiky")
        open({ "f(aa=g(b=1))" }, "julia", "# %s", "b=")
        assert.are.equal("aa", register_after "viky")
        open({ "f(aa=g(bé=11))" }, "julia", "# %s", "bé")
        assert.are.equal("aa", register_after "vikiky")
        -- The selection of `ak` ends on the newline.
        open({ '{"x": {"aa":', "1}}" }, "json", "", "aa")
        assert.are.equal('"x": ', register_after "vakaky")
    end)
    it("python keyword arguments, separator and its whitespace", function()
        local lines = { "f(aa=1, bb = 22)" }
        open(lines, "python", "# %s", "22")
        assert.are.equal("bb", yank "ik")
        open(lines, "python", "# %s", "= ")
        assert.are.equal("22", yank "iv")
        open(lines, "python", "# %s", " =")
        assert.are.equal("22", yank "iv")
        open(lines, "python", "# %s", "bb")
        assert.are.same({ "f(aa=1, 22)" }, lines_after "dak")
        open(lines, "python", "# %s", "bb")
        assert.are.same({ "f(aa=1, bb)" }, lines_after "dav")
    end)
    it("no pair around the cursor: mini.ai's search", function()
        open({ "f(aa=1)" }, "python", "# %s", "f")
        assert.are.equal("1", yank "iv")
        open({ "f(aa=1)" }, "python", "# %s", "f")
        vim.b.miniai_config = { search_method = "cover" }
        assert.are.equal("", yank "iv")
        open({ "f(aa=1)" }, "python", "# %s", ")")
        vim.b.miniai_config = { search_method = "cover_or_prev" }
        assert.are.equal("1", yank "iv")
        open({ "f()", "", "", "x = 1" }, "python", "# %s", "f")
        vim.b.miniai_config = { n_lines = 1 }
        assert.are.equal("", yank "iv")
    end)
    it("value on the next line", function()
        open({ "aa:", "  bb: 1" }, "yaml", "# %s", "aa")
        assert.are.equal("aa:\n  ", yank "ak")
    end)

    -- One case per node that the after/queries add.
    it("python", function()
        open({ 'd = {"aa": 11}' }, "python", "# %s", "aa")
        assert.are.equal("11", yank "iv")
        open({ "def g(xx=11): pass" }, "python", "# %s", "xx")
        assert.are.equal("11", yank "iv")
        open({ "def g(yy: int = 22): pass" }, "python", "# %s", "yy")
        assert.are.equal("22", yank "iv")
    end)
    it("r: named arguments and parameter defaults, not positional arguments", function()
        open({ "f(aa = 11, 22)" }, "r", "# %s", "aa")
        assert.are.equal("11", yank "iv")
        open({ "f(aa = 11, 22)" }, "r", "# %s", "22")
        assert.are.equal("", yank "iv")
        open({ "function(xx = 11) xx" }, "r", "# %s", "xx")
        assert.are.equal("11", yank "iv")
    end)
    it("lua table fields, brackets in the key", function()
        local lines = { 't = { aa = 11, ["b"] = 2, [x] = 3, 4 }' }
        open(lines, "lua", "-- %s", "aa")
        assert.are.equal("11", yank "iv")
        open(lines, "lua", "-- %s", '"b"')
        assert.are.equal('["b"]', yank "ik")
        open(lines, "lua", "-- %s", '"b"')
        assert.are.same({ "t = { aa = 11, 2, [x] = 3, 4 }" }, lines_after "dak")
        open(lines, "lua", "-- %s", '"b"')
        assert.are.same({ 't = { aa = 11, ["b"], [x] = 3, 4 }' }, lines_after "dav")
        -- One pair for `[x] = 3`: the next enclosing one is `t = …`.
        open(lines, "lua", "-- %s", "x]")
        assert.are.equal("[x]", yank "ik")
        open(lines, "lua", "-- %s", "x]")
        assert.are.equal("t", yank "2ik")
    end)
    it("julia =>", function()
        open({ "Dict(:aa => 11)" }, "julia", "# %s", "aa")
        assert.are.equal("11", yank "iv")
    end)
    it("yaml flow mapping", function()
        open({ "a: {cc: 22}" }, "yaml", "# %s", "cc")
        assert.are.equal("22", yank "iv")
    end)
    it("toml: a trailing comment is no part of the pair", function()
        open({ "aa = 11 # cc" }, "toml", "# %s", "aa")
        assert.are.equal("11", yank "iv")
        open({ "aa = 11 # cc" }, "toml", "# %s", "cc")
        assert.are.equal("", yank "iv")
        open({ "bb = { xx = 11 }" }, "toml", "# %s", "xx")
        assert.are.equal("11", yank "iv")
    end)
    it("json", function()
        open({ '{"aa": 11}' }, "json", "", "aa")
        assert.are.equal("11", yank "iv")
    end)
    it("ini: the value without the whitespace after '='", function()
        open({ "[s]", "aa = 11 ; c", "bb=22" }, "dosini", "; %s", "aa")
        assert.are.equal("11 ; c", yank "iv")
    end)
    it("zsh assignments before a command", function()
        -- after/plugin/treesitter.lua does this in the config.
        vim.treesitter.language.register("zsh", "sh.zsh")
        open({ "aa=11 bb=22 cmd" }, "zsh", "# %s", "bb")
        assert.are.equal("22", yank "iv")
    end)
    it("typst named arguments and let", function()
        open({ "#f(aa: /* c */ 11)" }, "typst", "// %s", "aa")
        assert.are.equal("11", yank "iv")
        open({ "#let xx = 22" }, "typst", "// %s", "xx")
        assert.are.equal("22", yank "iv")
    end)
    it("make", function()
        open({ "AA := 11" }, "make", "# %s", "AA")
        assert.are.equal("11", yank "iv")
    end)
    it("latex key=value options", function()
        open({ "\\usepackage[aa=11]{p}" }, "tex", "% %s", "aa")
        assert.are.equal("11", yank "iv")
    end)
    it("bibtex", function()
        open({ "@article{k, title = {T}}" }, "bib", "% %s", "title")
        assert.are.equal("{T}", yank "iv")
    end)
    it("vim dictionary", function()
        open({ 'let d = {"aa": 11}' }, "vim", '" %s', "aa")
        assert.are.equal("11", yank "iv")
    end)
    it("python fence in markdown", function()
        open({ "text", "```python", "f(aa=11)", "```" }, "markdown", "<!-- %s -->", "aa")
        assert.are.equal("11", yank "iv")
    end)

    describe("without a tree-sitter pair query", function()
        it("the argument around the cursor", function()
            open({ "f(a: 1, b: 2)" }, "", "", "b")
            assert.are.equal("2", yank "iv")
        end)
        it("the line", function()
            open({ "export XX=11" }, "", "", "XX")
            assert.are.equal("11", yank "iv")
            open({ "aa = f( )" }, "", "", " )")
            assert.are.equal("f( )", yank "iv")
            open({ "- key: v" }, "", "", "v")
            assert.are.equal("key", yank "ik")
        end)
    end)
end)

describe("split_pair", function()
    local split_pair = require("utils/init").split_pair
    --- @param text string
    --- @return string[]|nil key and value
    local function parts(text)
        local pair = split_pair(text)
        return pair and { text:sub(pair.key[1], pair.key[2]), text:sub(pair.value[1], pair.value[2]) }
    end

    it("splits at the first separator", function()
        assert.are.same({ "A", "1" }, parts "A := 1")
        assert.are.same({ "A", "1" }, parts "A?=1")
        assert.are.same({ "url", "http://x:8080" }, parts "url: http://x:8080")
        assert.are.same({ "X", "1" }, parts "export X=1")
        assert.are.same({ "key", "v" }, parts "- key: v")
    end)
    it("skips comparisons and separators in quotes and brackets", function()
        assert.is_nil(parts "a == b")
        assert.are.same({ '"a=b"', "c" }, parts '"a=b": c')
        assert.is_nil(parts "f(a: 1, b: 2)")
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
