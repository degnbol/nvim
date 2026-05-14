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

    it("with predicate, dispatch returns rhs when predicate passes", function()
        clear("zzpred")
        ab.iabbrev("zzpred", "passed", false, true, function() return true end)
        assert.are.equal("passed", ab._dispatch_lookup("zzpred"))
    end)

    it("with predicate, dispatch returns lhs when predicate fails", function()
        clear("zzpred2")
        ab.iabbrev("zzpred2", "passed", false, true, function() return false end)
        assert.are.equal("zzpred2", ab._dispatch_lookup("zzpred2"))
    end)
end)

describe("tex prose predicate (via plugin/abbreviations.lua)", function()
    -- vimtex isn't loaded under minimal_init, so we mock the three functions
    -- the predicate calls. Each test installs its own mock state.
    local plugin_loaded = false
    local function ensure_loaded()
        if plugin_loaded then return end
        dofile(vim.fn.expand("~/dotfiles/config/nvim/plugin/abbreviations.lua"))
        plugin_loaded = true
    end

    --- @param state { math?: boolean, in_doc?: boolean, env?: string, cmd?: string }
    local function mock(state)
        vim.fn["vimtex#syntax#in_mathzone"] = function() return state.math and 1 or 0 end
        vim.fn["vimtex#env#is_inside"] = function(_)
            return state.in_doc and { 1, 1 } or { 0, 0 }
        end
        vim.fn["vimtex#env#get_inner"] = function()
            return state.env and { name = state.env } or {}
        end
        vim.fn["vimtex#cmd#get_current"] = function()
            return state.cmd and { name = state.cmd } or {}
        end
    end

    -- After mocking, re-fire the tex FileType autocmd to register iabbrevs.
    -- We dispatch on a known lhs (`algo`) through _dispatch_lookup, which
    -- triggers the predicate.
    local function dispatch(state)
        ensure_loaded()
        mock(state)
        local buf = vim.api.nvim_create_buf(false, true)
        vim.api.nvim_set_current_buf(buf)
        pcall(function() vim.bo[buf].filetype = "tex" end)
        return ab._dispatch_lookup("algo")
    end

    it("body prose expands (in document, no env, no cmd)", function()
        assert.are_not.equal("algo", dispatch({ in_doc = true }))
    end)
    it("preamble skips (not in document)", function()
        assert.are.equal("algo", dispatch({ in_doc = false }))
    end)
    it("math skips", function()
        assert.are.equal("algo", dispatch({ in_doc = true, math = true }))
    end)
    it("verbatim env skips (not in prose env list)", function()
        assert.are.equal("algo", dispatch({ in_doc = true, env = "verbatim" }))
    end)
    it("itemize env expands", function()
        assert.are_not.equal("algo", dispatch({ in_doc = true, env = "itemize" }))
    end)
    it("figure env + caption cmd expands", function()
        assert.are_not.equal("algo",
            dispatch({ in_doc = true, env = "figure", cmd = "caption" }))
    end)
    it("figure env + includegraphics cmd skips", function()
        assert.are.equal("algo",
            dispatch({ in_doc = true, env = "figure", cmd = "includegraphics" }))
    end)
    it("label cmd skips", function()
        assert.are.equal("algo", dispatch({ in_doc = true, cmd = "label" }))
    end)
    it("section cmd expands", function()
        assert.are_not.equal("algo", dispatch({ in_doc = true, cmd = "section" }))
    end)
    it("starred section* cmd expands", function()
        assert.are_not.equal("algo", dispatch({ in_doc = true, cmd = "section*" }))
    end)
    it("starred figure* env + caption expands", function()
        assert.are_not.equal("algo",
            dispatch({ in_doc = true, env = "figure*", cmd = "caption" }))
    end)
end)

describe("typst prose predicate (via plugin/abbreviations.lua)", function()
    -- These tests load the abbreviations plugin so the FileType autocmd fires.
    -- The typst predicate uses treesitter to detect markup vs code mode.
    local plugin_loaded = false
    local function ensure_loaded()
        if plugin_loaded then return end
        dofile(vim.fn.expand("~/dotfiles/config/nvim/plugin/abbreviations.lua"))
        plugin_loaded = true
    end

    local function dispatch_at(lines, lhs, row, col)
        ensure_loaded()
        local buf = vim.api.nvim_create_buf(false, true)
        vim.api.nvim_buf_set_lines(buf, 0, -1, false, lines)
        vim.api.nvim_set_current_buf(buf)
        pcall(function() vim.bo[buf].filetype = "typst" end)
        vim.treesitter.get_parser(buf, "typst"):parse(true)
        vim.api.nvim_win_set_cursor(0, { row + 1, col })
        return ab._dispatch_lookup(lhs)
    end

    -- Find the cursor column just after `lhs` in the first matching line.
    local function find_after(lines, lhs)
        for i, line in ipairs(lines) do
            local idx = line:find(lhs, 1, true)
            if idx then return i - 1, idx - 1 + #lhs end
        end
    end

    local function expand_check(label, lines, lhs)
        it(label .. " — expands", function()
            local row, col = find_after(lines, lhs)
            local got = dispatch_at(lines, lhs, row, col)
            assert.are_not.equal(lhs, got)
        end)
    end

    local function skip_check(label, lines, lhs)
        it(label .. " — skips", function()
            local row, col = find_after(lines, lhs)
            assert.are.equal(lhs, dispatch_at(lines, lhs, row, col))
        end)
    end

    expand_check("plain markup",       { "Prose with algo here" }, "algo")
    expand_check("content arg [...]",  { "#text(size: 12pt)[Prose algo here]" }, "algo")
    expand_check("nested content arg", { "#text[Body algo more]" }, "algo")

    skip_check("inline math $...$",    { "Math: $algo$ done" }, "algo")
    skip_check("raw block ``` ```",    { "```", "algo = 1", "```" }, "algo")
    skip_check("raw span `...`",       { "Inline `algo` raw" }, "algo")
    skip_check("#let binding",         { "#let algo = 1" }, "algo")
    skip_check("#import path",         { '#import "algo.typ"' }, "algo")
    skip_check("#show expression",     { "#show heading: algo" }, "algo")
    skip_check("#set arg value",       { "#set page(margin: algo)" }, "algo")
    skip_check("function call name",   { "#algo(x, y)" }, "algo")
    skip_check("code arg in call",     { "#text(weight: algo)[Body]" }, "algo")
end)
