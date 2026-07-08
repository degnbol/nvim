---@diagnostic disable: undefined-global
local gloss = require "typst_glossary"

describe("typst_glossary", function()
    it("ref_key strips @ and modifier suffix", function()
        assert.are.equal("ts", gloss.ref_key("@ts"))
        assert.are.equal("ts", gloss.ref_key("@ts:pl"))
        assert.are.equal("ts", gloss.ref_key("@ts:short"))
        assert.are.equal("fig", gloss.ref_key("@fig:f"))
        assert.are.equal("π-bond", gloss.ref_key("@π-bond"))
    end)

    it("entries finds dict pairs by {row, col}, dequotes keys, skips string fields", function()
        local content = table.concat({
            "#let g = (",
            '  ts: ( short: "TS" ),',
            '  "1-2-shift": ( short: "shift" ),',
            "  plain: 3,",
            ")",
        }, "\n")
        local positions = gloss.entries(content)
        assert.are.same({ 1, 2 }, positions["ts"])
        assert.are.same({ 2, 2 }, positions["1-2-shift"])
        assert.is_nil(positions["short"])
        assert.is_nil(positions["plain"])
    end)

    it("imports keeps relative file paths, drops package imports", function()
        local content = table.concat({
            '#import "preamble.typ": *',
            '#import "./glossary.typ": glossary-local',
            '#import "@preview/glossy:0.9.0": glossary',
            '#import "@local/degnlib:0.1.0": style',
        }, "\n")
        assert.are.same(
            { "/doc/preamble.typ", "/doc/glossary.typ" },
            gloss.imports(content, "/doc")
        )
    end)
end)
