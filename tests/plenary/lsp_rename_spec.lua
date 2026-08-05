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
