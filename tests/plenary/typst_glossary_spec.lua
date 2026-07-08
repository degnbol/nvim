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

    it("fields extracts an entry's string/content values, stripping brackets", function()
        local content = table.concat({
            "#let g = (",
            '  mep: (',
            '    short: "MEP",',
            "    long-fmt: [2-#emph[C]-methyl],",
            "    description: [The @ts pathway.],",
            "  ),",
            "  ts: ( short: \"TS\", long: \"terpene synthase\" ),",
            ")",
        }, "\n")
        local mep = gloss.fields(content, "mep")
        assert.are.same({ text = "MEP", kind = "string" }, mep["short"])
        assert.are.same({ text = "2-#emph[C]-methyl", kind = "content" }, mep["long-fmt"])
        assert.are.same({ text = "The @ts pathway.", kind = "content" }, mep["description"])
        assert.are.same({ text = "terpene synthase", kind = "string" }, gloss.fields(content, "ts")["long"])
        assert.is_nil(gloss.fields(content, "absent"))
    end)

    it("render prefers -fmt forms, keeps short even if it differs from the key", function()
        local md = gloss.render({
            short = { text = "kcat", kind = "string" },
            ["short-fmt"] = { text = '$k_"cat"$', kind = "content" },
            long = { text = "turnover number", kind = "string" },
        })
        assert.are.equal('```typst\n$k_"cat"$ — turnover number\n```', md)
    end)

    it("render fences content descriptions, omits the dash when no long form", function()
        assert.are.equal(
            "```typst\nMEP\n```\n\n```typst\nThe @ts pathway.\n```",
            gloss.render({
                short = { text = "MEP", kind = "string" },
                description = { text = "The @ts pathway.", kind = "content" },
            })
        )
        assert.are.equal(
            "```typst\nTS\n```\n\nA plain string description.",
            gloss.render({
                short = { text = "TS", kind = "string" },
                description = { text = "A plain string description.", kind = "string" },
            })
        )
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
