---@diagnostic disable: undefined-global
-- Tests for lua/utils/iabbrev.lua — vim-abolish :Abolish replacement.

local ab = require("utils.iabbrev")

describe("iabbrev.mixedcase", function()
    it("uppercases first char of all-lowercase", function()
        assert.are.equal("Hello", ab.mixedcase("hello"))
    end)
    it("uppercases first char and preserves mixed case", function()
        -- camelcase lowercases first → "tHere", mixedcase re-uppers → "THere"
        assert.are.equal("THere", ab.mixedcase("THere"))
        assert.are.equal("Th", ab.mixedcase("Th"))
    end)
    it("turns ALL-UPPER into Title", function()
        assert.are.equal("Abc", ab.mixedcase("ABC"))
    end)
    it("converts snake_case to MixedCase", function()
        assert.are.equal("FooBar", ab.mixedcase("foo_bar"))
        assert.are.equal("FooBar", ab.mixedcase("FOO_BAR"))
    end)
    it("treats dash like underscore", function()
        assert.are.equal("FooBar", ab.mixedcase("foo-bar"))
    end)
end)

describe("iabbrev.expand_braces", function()
    it("returns the input unchanged when there are no braces", function()
        assert.are.same({ { "foo", "bar" } }, ab.expand_braces("foo", "bar"))
    end)

    it("expands single brace in lhs with literal rhs", function()
        -- No brace in rhs → all expansions map to literal rhs.
        local pairs = ab.expand_braces("liek", "like")
        assert.are.same({ { "liek", "like" } }, pairs)
    end)

    it("pairs lhs and rhs braces by index", function()
        -- paren → parenthesis, parens → parentheses
        local pairs = ab.expand_braces("paren{,s}", "parenthes{i,e}s")
        assert.are.same({
            { "paren", "parenthesis" },
            { "parens", "parentheses" },
        }, pairs)
    end)

    it("treats empty {} in rhs as 'copy lhs target verbatim'", function()
        local pairs = ab.expand_braces("TH{ere,en,e,is}", "Th{}")
        assert.are.same({
            { "THere", "There" },
            { "THen", "Then" },
            { "THe", "The" },
            { "THis", "This" },
        }, pairs)
    end)

    it("expands multiple braces recursively", function()
        -- {a,b}{c,d} pairs index-wise: ac→ec, bd→fd
        local pairs = ab.expand_braces("{a,b}{c,d}", "{e,f}{}")
        assert.are.same({
            { "ac", "ec" },
            { "ad", "ed" },
            { "bc", "fc" },
            { "bd", "fd" },
        }, pairs)
    end)

    it("handles empty target in middle (leading comma)", function()
        local pairs = ab.expand_braces("aa{,s}", "amino acid{,s}")
        assert.are.same({
            { "aa", "amino acid" },
            { "aas", "amino acids" },
        }, pairs)
    end)

    it("wraps rhs replacements modulo their length", function()
        -- 3 targets, 1 replacement → all map to that replacement
        local pairs = ab.expand_braces("{a,b,c}", "{x}")
        assert.are.same({
            { "a", "x" },
            { "b", "x" },
            { "c", "x" },
        }, pairs)
    end)
end)

describe("iabbrev.case_variants", function()
    it("emits 3 distinct keys for all-lowercase input", function()
        local out = ab.case_variants("liek", "like")
        assert.are.equal("like", out.liek)
        assert.are.equal("Like", out.Liek)
        assert.are.equal("LIKE", out.LIEK)
    end)

    it("includes the literal lhs alongside variants", function()
        local out = ab.case_variants("THere", "There")
        -- mixedcase("THere") = "THere" (idempotent), so literal and mixed overlap.
        assert.are.equal("There", out.THere)
        assert.are.equal("there", out.there)
        assert.are.equal("THERE", out.THERE)
    end)
end)

describe("iabbrev.iabbrev (side effects)", function()
    -- Use a unique prefix per test so we don't collide with real abbreviations
    -- or pollute global state for later tests.
    local function clear(lhs)
        pcall(function() vim.cmd("iunabbrev " .. lhs) end)
    end

    it("creates an iabbrev visible to :iabbrev <lhs>", function()
        clear("zzliek"); clear("zzLiek"); clear("zzLIEK")
        ab.iabbrev("zzliek", "like")
        local out = vim.api.nvim_exec2("iabbrev zzliek", { output = true }).output
        assert.is_truthy(out:find("like"))
    end)

    it("with cases=true defines the UPPER variant too", function()
        clear("zzfoo"); clear("Zzfoo"); clear("ZZFOO")
        ab.iabbrev("zzfoo", "bar", true)
        local out = vim.api.nvim_exec2("iabbrev ZZFOO", { output = true }).output
        assert.is_truthy(out:find("BAR"))
    end)

    it("with cases=false defines only the literal", function()
        clear("zzIm"); clear("zzim"); clear("zzIM")
        ab.iabbrev("zzIm", "I'm", false)
        local lower = vim.api.nvim_exec2("iabbrev zzim", { output = true }).output
        assert.is_truthy(lower:find("No abbreviation") or lower == "")
    end)
end)
