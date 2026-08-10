---@diagnostic disable: undefined-global
local google = require "docstring.style.google"
local model = require "docstring.model"

--- A model with one entry of every kind, exercising typed + typeless params, a
--- multiline description, a return, a raise, an unmodeled section, and an orphan.
local function full_model()
    local m = model.empty()
    m.summary = "Do a thing.\n\nA longer explanation\nover two lines."
    m.params = {
        { name = "x", type = nil, description = "The x value." },
        { name = "y", type = "int", description = "The y value,\ncontinued on a second line." },
        { name = "*args", type = nil, description = "Extra positionals." },
    }
    m.returns = { type = "bool", description = "Whether it worked." }
    m.raises = { { type = "ValueError", description = "When x is negative." } }
    m.other = { { header = "Example", body = ">>> f(1, 2)\nTrue" } }
    m.unmatched = { { name = "old", type = nil, description = "Removed from the signature." } }
    return m
end

describe("docstring google style", function()
    it("round-trips a full model: parse(render(m)) == m", function()
        local m = full_model()
        assert.are.same(m, google.parse(google.render(m)))
    end)

    it("round-trips an empty model", function()
        local m = model.empty()
        assert.are.same(m, google.parse(google.render(m)))
    end)

    it("round-trips a summary-only model", function()
        local m = model.empty()
        m.summary = "Just a summary."
        assert.are.same(m, google.parse(google.render(m)))
    end)

    it("parses a typeless Google docstring", function()
        local text = table.concat({
            "Add two numbers.",
            "",
            "Args:",
            "    a: First addend.",
            "    b: Second addend.",
            "",
            "Returns:",
            "    The sum.",
        }, "\n")
        local m = google.parse(text)
        assert.are.equal("Add two numbers.", m.summary)
        assert.are.same({ name = "a", type = nil, description = "First addend." }, m.params[1])
        assert.are.same({ name = "b", type = nil, description = "Second addend." }, m.params[2])
        assert.are.same({ type = nil, description = "The sum." }, m.returns)
    end)

    it("parses a typed Args entry with a parenthesised type", function()
        local text = "S.\n\nArgs:\n    n (int): A count."
        local m = google.parse(text)
        assert.are.same({ name = "n", type = "int", description = "A count." }, m.params[1])
    end)

    it("parses a multiline description by indent continuation", function()
        local text = "S.\n\nArgs:\n    x: First line.\n        Second line."
        local m = google.parse(text)
        assert.are.equal("First line.\nSecond line.", m.params[1].description)
    end)

    it("routes an orphan section back to unmatched, not other", function()
        local text = "S.\n\n" .. google.ORPHAN_HEADER .. ":\n    gone: No longer a param."
        local m = google.parse(text)
        assert.are.equal(0, #m.other)
        assert.are.same({ name = "gone", type = nil, description = "No longer a param." }, m.unmatched[1])
    end)

    it("keeps an unmodeled section verbatim in other", function()
        local text = "S.\n\nNote:\n    Be careful."
        local m = google.parse(text)
        assert.are.same({ header = "Note", body = "Be careful." }, m.other[1])
    end)

    it("does not mistake a colon in a return description for a type", function()
        local m = model.empty()
        m.summary = "S."
        m.returns = { type = nil, description = "The parsed config: see below." }
        assert.are.same(m, google.parse(google.render(m)))
    end)

    it("round-trips a subscripted return type", function()
        local m = model.empty()
        m.summary = "S."
        m.returns = { type = "dict[str, int]", description = "The counts." }
        assert.are.same(m, google.parse(google.render(m)))
    end)

    it("parses a simple return type but leaves prose alone", function()
        assert.are.same(
            { type = "bool", description = "Whether it worked." },
            google.parse("S.\n\nReturns:\n    bool: Whether it worked.").returns
        )
        assert.are.same(
            { type = nil, description = "See the note: it varies." },
            google.parse("S.\n\nReturns:\n    See the note: it varies.").returns
        )
    end)

    it("renders a blank line inside a description without trailing whitespace", function()
        local m = model.empty()
        m.summary = "S."
        m.params = { { name = "x", type = nil, description = "First.\n\nAfter a gap." } }
        local text = google.render(m)
        for _, line in ipairs(vim.split(text, "\n", { plain = true })) do
            assert.is_falsy(line:match("^%s+$"), "blank line has trailing whitespace: " .. vim.inspect(line))
        end
        assert.are.same(m, google.parse(text))
    end)

    it("round-trips a dotted exception type in Raises", function()
        local m = model.empty()
        m.summary = "S."
        m.raises = { { type = "socket.error", description = "On connection failure." } }
        assert.are.same(m, google.parse(google.render(m)))
    end)
end)

describe("docstring locate_doc_entries", function()
    it("returns name-token ranges for Args and orphan entries", function()
        local text = table.concat({
            "S.",
            "",
            "Args:",
            "    alpha: A.",
            "    beta: B.",
        }, "\n")
        local located = google.locate_doc_entries(text)
        assert.are.equal(2, #located)
        assert.are.same({ name = "alpha", srow = 3, scol = 4, erow = 3, ecol = 9 }, located[1])
        assert.are.same({ name = "beta", srow = 4, scol = 4, erow = 4, ecol = 8 }, located[2])
    end)

    it("locates the name only, not the type or description", function()
        local text = "S.\n\nArgs:\n    n (int): A count."
        local located = google.locate_doc_entries(text)
        assert.are.same({ name = "n", srow = 3, scol = 4, erow = 3, ecol = 5 }, located[1])
    end)

    it("skips Raises entries", function()
        local text = "S.\n\nRaises:\n    ValueError: Bad."
        assert.are.equal(0, #google.locate_doc_entries(text))
    end)
end)
