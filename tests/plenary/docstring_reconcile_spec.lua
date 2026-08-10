---@diagnostic disable: undefined-global
local reconcile = require("docstring.reconcile").reconcile
local model = require "docstring.model"

--- A model whose params are `entries`, summary "S.".
local function with_params(entries)
    local m = model.empty()
    m.summary = "S."
    m.params = entries
    return m
end

--- SigParam list from `{ name, type }` pairs (type optional).
local function sig(...)
    local params = {}
    for _, spec in ipairs({ ... }) do
        params[#params + 1] = { name = spec[1], type = spec[2] }
    end
    return params
end

describe("docstring reconcile params", function()
    it("appends a fresh typeless entry for a param added at the end", function()
        local m = with_params({
            { name = "a", type = nil, description = "First." },
            { name = "b", type = nil, description = "Second." },
        })
        local out = reconcile(sig({ "a" }, { "b" }, { "c" }), m)
        assert.are.same({
            { name = "a", type = nil, description = "First." },
            { name = "b", type = nil, description = "Second." },
            { name = "c", type = nil, description = "" },
        }, out.params)
    end)

    it("inserts a fresh entry in the middle at signature position", function()
        local m = with_params({
            { name = "a", type = nil, description = "First." },
            { name = "b", type = nil, description = "Second." },
        })
        local out = reconcile(sig({ "a" }, { "c" }, { "b" }), m)
        assert.are.same({ "a", "c", "b" }, vim.tbl_map(function(p) return p.name end, out.params))
        assert.are.equal("", out.params[2].description)
    end)

    it("reorders doc entries to signature order, keeping descriptions", function()
        local m = with_params({
            { name = "a", type = nil, description = "First." },
            { name = "b", type = nil, description = "Second." },
        })
        local out = reconcile(sig({ "b" }, { "a" }), m)
        assert.are.same({
            { name = "b", type = nil, description = "Second." },
            { name = "a", type = nil, description = "First." },
        }, out.params)
    end)

    it("surfaces a param absent from the signature as an orphan", function()
        local m = with_params({
            { name = "a", type = nil, description = "Kept." },
            { name = "old", type = nil, description = "Removed from signature." },
        })
        local out = reconcile(sig({ "a" }), m)
        assert.are.same({ { name = "a", type = nil, description = "Kept." } }, out.params)
        assert.are.same({ { name = "old", type = nil, description = "Removed from signature." } }, out.unmatched)
    end)

    it("reclaims an orphan's description when its param returns to the signature", function()
        local m = with_params({ { name = "a", type = nil, description = "Kept." } })
        m.unmatched = { { name = "old", type = nil, description = "Was orphaned." } }
        local out = reconcile(sig({ "a" }, { "old" }), m)
        assert.are.same({
            { name = "a", type = nil, description = "Kept." },
            { name = "old", type = nil, description = "Was orphaned." },
        }, out.params)
        assert.are.same({}, out.unmatched)
    end)

    it("preserves a documented type when the signature param is unannotated", function()
        local m = with_params({ { name = "x", type = "int", description = "A count." } })
        local out = reconcile(sig({ "x" }), m)
        assert.are.equal("int", out.params[1].type)
    end)

    it("does not inject a signature type into a typeless entry", function()
        local m = with_params({ { name = "x", type = nil, description = "A count." } })
        local out = reconcile(sig({ "x", "int" }), m)
        assert.is_nil(out.params[1].type)
    end)

    it("refreshes a documented type from the signature annotation", function()
        local m = with_params({ { name = "x", type = "float", description = "A count." } })
        local out = reconcile(sig({ "x", "int" }), m)
        assert.are.equal("int", out.params[1].type)
    end)

    it("generates fresh typeless entries from a signature with no docstring", function()
        local out = reconcile(sig({ "a" }, { "b", "int" }), model.empty())
        assert.are.equal("", out.summary)
        assert.are.same({
            { name = "a", type = nil, description = "" },
            { name = "b", type = nil, description = "" },
        }, out.params)
    end)

    it("passes summary, returns, raises and other through unchanged", function()
        local m = with_params({ { name = "a", type = nil, description = "A." } })
        m.returns = { type = "bool", description = "Ok." }
        m.raises = { { type = "ValueError", description = "Bad." } }
        m.other = { { header = "Example", body = ">>> f()" } }
        local out = reconcile(sig({ "a" }), m)
        assert.are.equal("S.", out.summary)
        assert.are.same({ type = "bool", description = "Ok." }, out.returns)
        assert.are.same({ { type = "ValueError", description = "Bad." } }, out.raises)
        assert.are.same({ { header = "Example", body = ">>> f()" } }, out.other)
    end)

    it("does not mutate the input model's entries", function()
        local entry = { name = "x", type = "float", description = "A count." }
        local m = with_params({ entry })
        reconcile(sig({ "x", "int" }), m)
        assert.are.equal("float", entry.type)
    end)

    it("is idempotent: a second reconcile is a no-op", function()
        local m = with_params({
            { name = "b", type = "float", description = "Second." },
            { name = "old", type = nil, description = "Gone." },
        })
        m.returns = { type = "int", description = "Count." }
        local params = sig({ "a" }, { "b", "int" })
        local once = reconcile(params, m)
        local twice = reconcile(params, once)
        assert.are.same(once, twice)
    end)
end)
