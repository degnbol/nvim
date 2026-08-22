---@diagnostic disable: undefined-global
-- Tests for after/queries/perl/injections.scm (regex nested in every pattern)
-- and the interpolation repaint in after/queries/perl/highlights.scm. The
-- zsh -> perl step that feeds them lives in zsh_injections_spec.lua.

-- minimal_init.lua only adds the config dir to rtp; its after/ subdir (where
-- both queries live) is not searched unless appended explicitly.
vim.opt.runtimepath:append(vim.fn.getcwd() .. "/after")

--- Text of every regex region injected into `src`, parsed as perl.
--- @param src string
--- @return string[]
local function regex_regions(src)
    local parser = vim.treesitter.get_string_parser(src, "perl")
    parser:parse(true)
    local regions = {}
    for lang, child in pairs(parser:children()) do
        if lang == "regex" then
            for _, tree in ipairs(child:trees()) do
                regions[#regions + 1] =
                    vim.treesitter.get_node_text(tree:root(), src)
            end
        end
    end
    return regions
end

--- Assert `src` injects exactly one regex region, holding `expected`.
local function assert_regex(src, expected)
    local regions = regex_regions(src)
    assert.are.equal(1, #regions,
        "expected 1 regex region, got " .. vim.inspect(regions))
    assert.are.equal(expected, regions[1])
end

--- Highlight captures landing on a node whose text equals `text`, parsing `src`
--- as perl. Each entry is { group = "...", priority = <number|nil> }.
--- @param src string
--- @param text string
--- @return { group: string, priority: number|nil }[]
local function captures_on(src, text)
    local root = vim.treesitter.get_string_parser(src, "perl")
        :parse(true)[1]:root()
    local query = vim.treesitter.query.get("perl", "highlights")
    local captures = {}
    for id, node, meta in query:iter_captures(root, src) do
        if vim.treesitter.get_node_text(node, src) == text then
            captures[#captures + 1] = {
                group = query.captures[id],
                priority = tonumber(meta.priority),
            }
        end
    end
    return captures
end

--- True if `captures` holds a capture of `group` at `priority`.
local function has_capture(captures, group, priority)
    for _, capture in ipairs(captures) do
        if capture.group == group and capture.priority == priority then
            return true
        end
    end
    return false
end

describe("perl regex injection", function()
    it("injects a substitution pattern", function()
        assert_regex([[s/\d+/x/;]], [[\d+]])
    end)

    it("injects a match pattern", function()
        assert_regex([[if ($x =~ /^a[bc]+$/) { }]], [[^a[bc]+$]])
    end)

    it("injects a qr// pattern", function()
        assert_regex([[my $re = qr/\d{2,}/i;]], [[\d{2,}]])
    end)

    it("injects a pattern in braces: s{…}{…}", function()
        assert_regex([[s{\s+}{ }g;]], [[\s+]])
    end)

    it("excludes the delimiters and the replacement", function()
        assert_regex([[s/a/b/g;]], "a")
    end)

    it("keeps the pattern contiguous across escapes", function()
        assert_regex(
            [[s/(`[\w\/.-]+\.(?:lua|md)):[\d,-]+`/$1`/g;]],
            [[(`[\w\/.-]+\.(?:lua|md)):[\d,-]+`]]
        )
    end)

    it("does not inject a transliteration character map", function()
        assert.are.equal(0, #regex_regions([[tr/a-z/A-Z/;]]))
    end)
end)

describe("perl interpolation repaint", function()
    it("repaints an interpolated scalar at priority 101", function()
        assert.is_true(
            has_capture(captures_on([[s/$foo\d+/bar/;]], "$foo"), "variable", 101),
            "expected @variable at priority 101 on $foo")
    end)

    it("repaints element access at priority 101", function()
        assert.is_true(
            has_capture(captures_on([[s/$h{k}/x/;]], "$h{k}"), "variable", 101),
            "expected @variable at priority 101 on $h{k}")
    end)

    it("leaves a pattern with no interpolation to the regex layer", function()
        for _, capture in ipairs(captures_on([[s/\d+/x/;]], [[\d+]])) do
            assert.are_not.equal(101, capture.priority)
        end
    end)
end)
