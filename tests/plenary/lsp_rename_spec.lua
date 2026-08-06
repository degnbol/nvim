---@diagnostic disable: undefined-global
local argparse = require "lsp_rename.argparse"
local lsp = require "utils.lsp"

describe("utils.lsp.workspace_edit_files", function()
    local edit = { range = {}, newText = "x" }

    it("normalizes the `changes` shape", function()
        local files = lsp.workspace_edit_files({ changes = { ["file:///a.py"] = { edit } } })
        assert.are.equal(1, #files)
        assert.are.equal("file:///a.py", files[1].uri)
        assert.are.same({ edit }, files[1].edits)
    end)

    it("normalizes the `documentChanges` shape", function()
        local files = lsp.workspace_edit_files({
            documentChanges = {
                { textDocument = { uri = "file:///a.py", version = 1 }, edits = { edit } },
            },
        })
        assert.are.equal("file:///a.py", files[1].uri)
        assert.are.same({ edit }, files[1].edits)
    end)

    it("skips resource operations without edits", function()
        local files = lsp.workspace_edit_files({
            documentChanges = {
                { kind = "create", uri = "file:///new.py" },
                { textDocument = { uri = "file:///a.py", version = 1 }, edits = { edit } },
            },
        })
        assert.are.equal(1, #files)
        assert.are.equal("file:///a.py", files[1].uri)
    end)

    it("returns the live edits array so callers can append", function()
        local we = { changes = { ["file:///a.py"] = {} } }
        local files = lsp.workspace_edit_files(we)
        table.insert(files[1].edits, edit)
        assert.are.same({ edit }, we.changes["file:///a.py"])
    end)
end)

describe("lsp_rename.argparse dest derivation", function()
    -- Shared truth with the argparse diagnostic: call snippet → derived dest.
    local cases = {
        ['add_argument("infiles")'] = "infiles",
        ['add_argument("-o", "--out")'] = "out",
        ['add_argument("--out-file")'] = "out_file",
        ['add_argument("-o")'] = "o",
        ['add_argument("-o", "--out", "--alias")'] = "out",
        ['add_argument("--out", dest="target")'] = "target",
        ['add_argument("infiles", dest="target")'] = "target",
    }
    for src, dest in pairs(cases) do
        it(src .. " → " .. dest, function()
            assert.are.equal(dest, argparse.dest_of(src))
        end)
    end
end)

describe("lsp_rename.argparse edits", function()
    it("rewrites a positional name string", function()
        local edits = argparse.edits('p.add_argument("infiles", nargs="+")', "infiles", "inputs")
        assert.are.equal(1, #edits)
        assert.are.same({ srow = 0, scol = 16, erow = 0, ecol = 23, newText = "inputs" }, edits[1])
    end)

    it("rewrites the long flag a flag-derived dest came from", function()
        local edits = argparse.edits('p.add_argument("-o", "--out")', "out", "output")
        assert.are.equal(1, #edits)
        assert.are.equal("--output", edits[1].newText)
        -- targets the `--out` content, leaving `-o` alone
        assert.are.equal(22, edits[1].scol)
    end)

    it("maps `_` to `-` when rewriting a flag", function()
        local edits = argparse.edits('p.add_argument("--out")', "out", "out_file")
        assert.are.equal("--out-file", edits[1].newText)
    end)

    it("rewrites only dest= when present, never the option string", function()
        local edits = argparse.edits('p.add_argument("--out", dest="out")', "out", "output")
        assert.are.equal(1, #edits)
        assert.are.equal("output", edits[1].newText)
        assert.are.equal(30, edits[1].scol) -- the "out" inside dest=, not --out
    end)

    it("ignores calls whose dest does not match the old name", function()
        local edits = argparse.edits('p.add_argument("--other")', "out", "output")
        assert.are.same({}, edits)
    end)

    it("never touches the option string when dest= is non-literal", function()
        local edits = argparse.edits('p.add_argument("--out", dest=cfg.name)', "out", "output")
        assert.are.same({}, edits)
    end)

    it("ignores f-string option names it cannot resolve", function()
        local edits = argparse.edits('p.add_argument(f"--{x}")', "out", "output")
        assert.are.same({}, edits)
    end)
end)

describe("lsp_rename.argparse dest_at", function()
    -- Dest configured by the option string the cursor sits in. `mark` is a
    -- unique substring of the source; the cursor is placed at its first byte.
    local function dest_at(src, mark)
        local byte = assert(src:find(mark, 1, true), "mark not found")
        local root = vim.treesitter.get_string_parser(src, "python"):parse()[1]:root()
        local node = assert(root:named_descendant_for_range(0, byte - 1, 0, byte - 1))
        return argparse.dest_at(node, src)
    end

    it("resolves from the long flag string", function()
        assert.are.equal("cool_value", dest_at('p.add_argument("-v", "--cool-value")', "cool-value"))
    end)

    it("resolves from a different option string of the same call", function()
        -- cursor on `-v`, dest comes from `--cool-value` — proves containment
        -- checks the full option set, not just the token that derived the dest.
        assert.are.equal("cool_value", dest_at('p.add_argument("-v", "--cool-value")', "-v"))
    end)

    it("resolves from a positional string", function()
        assert.are.equal("infiles", dest_at('p.add_argument("infiles")', "infiles"))
    end)

    it("resolves the dest from either the option string or the dest= value", function()
        assert.are.equal("target", dest_at('p.add_argument("--out", dest="target")', "--out"))
        assert.are.equal("target", dest_at('p.add_argument("--out", dest="target")', "target"))
    end)

    it("returns nil for a non-option value string (nargs, default)", function()
        assert.is_nil(dest_at('p.add_argument("--out", nargs="+")', "+"))
        assert.is_nil(dest_at('p.add_argument("--out", default="foo")', "foo"))
    end)

    it("returns nil for a string in a non-add_argument call", function()
        assert.is_nil(dest_at('self.helper("cool_value")', "cool_value"))
    end)

    it("returns nil when the cursor is not on a string", function()
        assert.is_nil(dest_at('p.add_argument("--out")', "add_argument"))
    end)

    it("returns nil for a non-literal dest= option string", function()
        assert.is_nil(dest_at('p.add_argument("--out", dest=cfg.name)', "--out"))
    end)
end)

describe("lsp_rename.anchor.position", function()
    local anchor = require "lsp_rename.anchor"

    -- Load `lines` into a python buffer and return the anchor position for the
    -- cursor at the first byte of `mark` (a unique substring). Returns a
    -- `{ row, col }` table or nil, for stable `assert.are.same`.
    local function position(lines, mark)
        local buf = vim.api.nvim_create_buf(false, true)
        vim.api.nvim_buf_set_lines(buf, 0, -1, false, lines)
        vim.bo[buf].filetype = "python"
        vim.treesitter.get_parser(buf, "python"):parse()

        local row, col
        for i, line in ipairs(lines) do
            local byte = line:find(mark, 1, true)
            if byte then row, col = i - 1, byte - 1 break end
        end
        assert(row, "mark not found")

        local arow, acol = anchor.position(buf, row, col)
        vim.api.nvim_buf_delete(buf, { force = true })
        return arow and { arow, acol } or nil
    end

    local tap = {
        "class Args(Tap):",
        '    cool_value: str = "foo"',
        "    def configure(self):",
        '        self.add_argument("-v", "--cool-value")',
    }

    it("resolves an option string to its Tap attribute declaration", function()
        assert.are.same({ 1, 4 }, position(tap, "cool-value")) -- `cool_value` on line index 1, col 4
    end)

    it("resolves an annotation-only attribute (no default value)", function()
        local annotated = vim.deepcopy(tap)
        annotated[2] = "    cool_value: str"
        assert.are.same({ 1, 4 }, position(annotated, "cool-value"))
    end)

    it("returns nil when the Tap class has no matching attribute", function()
        local no_attr = vim.deepcopy(tap)
        no_attr[2] = '    other: str = "foo"'
        assert.is_nil(position(no_attr, "cool-value"))
    end)

    it("returns nil when the enclosing class does not subclass Tap", function()
        local not_tap = vim.deepcopy(tap)
        not_tap[1] = "class Args(object):"
        assert.is_nil(position(not_tap, "cool-value"))
    end)

    it("returns nil for a module-level add_argument (no enclosing class)", function()
        assert.is_nil(position({ 'parser.add_argument("--out")' }, "--out"))
    end)

    it("returns nil when the cursor is not on an option string", function()
        assert.is_nil(position(tap, "configure"))
    end)
end)
