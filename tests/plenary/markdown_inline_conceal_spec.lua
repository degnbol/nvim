---@diagnostic disable: undefined-global
-- Tests for after/queries/markdown_inline/highlights.scm:
--   1. Strikethrough conceal — only `~~double~~` (nested-strikethrough parse),
--      not stray single tildes.
--   2. Autolink bracket dimming — `<` and `>` of uri_autolink/email_autolink
--      get @punctuation.bracket via #head!/#tail! directives.

-- Source plugin/treesitter.lua to register #head!/#tail! directives.
-- plugin/* isn't reliably loaded before PlenaryBustedFile runs.
vim.cmd.source(vim.fn.getcwd() .. "/plugin/treesitter.lua")

--- Parse `src` as markdown, run the markdown_inline highlights query, and
--- return one entry per capture: { capture, sr, sc, er, ec, text, conceal,
--- priority }. Ranges reflect any `metadata.range` overrides set by
--- directives.
local function captures(src)
    local buf = vim.api.nvim_create_buf(false, true)
    vim.api.nvim_set_current_buf(buf)
    local lines = vim.split(src, "\n", { plain = true })
    vim.api.nvim_buf_set_lines(buf, 0, -1, false, lines)
    vim.bo[buf].filetype = "markdown"
    local md_parser = vim.treesitter.get_parser(buf, "markdown")
    md_parser:parse(true)
    local query = vim.treesitter.query.get("markdown_inline", "highlights")
    local out = {}
    md_parser:for_each_tree(function(tree, ltree)
        if ltree:lang() ~= "markdown_inline" then return end
        for id, node, metadata in query:iter_captures(tree:root(), buf) do
            local sr, sc, er, ec = node:range()
            local m = metadata[id]
            if m and m.range then
                sr, sc, _, er, ec, _ = unpack(m.range)
            end
            local ok, txt = pcall(vim.api.nvim_buf_get_text, buf, sr, sc, er, ec, {})
            local prio = (m and m.priority) or metadata.priority
            out[#out + 1] = {
                capture = query.captures[id],
                sr = sr, sc = sc, er = er, ec = ec,
                text = ok and table.concat(txt, "\n") or nil,
                conceal = metadata.conceal,
                priority = prio and tonumber(prio) or nil,
            }
        end
    end)
    pcall(vim.api.nvim_buf_delete, buf, { force = true })
    return out
end

describe("markdown_inline strikethrough conceal", function()
    it("conceals all four ~ of ~~double~~ at priority 250", function()
        local caps = captures("| | ~~strike~~ |")
        -- The `~~strike~~` text sits at cols 5–15. Find every conceal
        -- directive whose target is a `~` byte.
        local tilde_p250 = {}
        for _, c in ipairs(caps) do
            if c.text == "~" and c.conceal == "" and c.priority == 250 then
                tilde_p250[#tilde_p250 + 1] = c
            end
        end
        assert.are.equal(4, #tilde_p250,
            "expected 4 priority-250 conceal=\"\" captures on the four ~ chars, got "
            .. #tilde_p250)
    end)

    it("highlights inner strikethrough as @markup.strikethrough.double", function()
        local caps = captures("~~strike~~")
        local hits = {}
        for _, c in ipairs(caps) do
            if c.capture == "markup.strikethrough.double" then
                hits[#hits + 1] = c
            end
        end
        assert.are.equal(1, #hits)
        -- Range is the inner strikethrough (excludes the outermost ~ on each side).
        assert.are.equal("~strike~", hits[1].text)
    end)

    it("does not re-conceal single tildes", function()
        local caps = captures("~14 vs ~7")
        -- Discriminator: nested-only conceal rules must NOT fire on the
        -- single-tilde parser false-positive. The user's existing
        -- priority-200 unconceal is what keeps them visible; we just need
        -- to verify our new priority-250 rules don't match here.
        for _, c in ipairs(caps) do
            if c.text == "~" then
                assert.is_not.equal(250, c.priority,
                    "single-tilde ~ should not get priority-250 conceal")
            end
        end
    end)
end)

describe("markdown_inline autolink bracket dim", function()
    it("captures < and > of uri_autolink as @punctuation.bracket", function()
        local caps = captures("<https://example.com>")
        local punct = {}
        for _, c in ipairs(caps) do
            if c.capture == "punctuation.bracket" then
                punct[#punct + 1] = c
            end
        end
        assert.are.equal(2, #punct)
        -- Sort by start col for deterministic comparison.
        table.sort(punct, function(a, b) return a.sc < b.sc end)
        assert.are.equal("<", punct[1].text)
        assert.are.equal(">", punct[2].text)
    end)

    it("captures < and > of email_autolink as @punctuation.bracket", function()
        local caps = captures("<user@example.com>")
        local punct = {}
        for _, c in ipairs(caps) do
            if c.capture == "punctuation.bracket" then
                punct[#punct + 1] = c
            end
        end
        assert.are.equal(2, #punct)
        table.sort(punct, function(a, b) return a.sc < b.sc end)
        assert.are.equal("<", punct[1].text)
        assert.are.equal(">", punct[2].text)
    end)

    it("captures inline-link punctuation as @punctuation.bracket", function()
        local caps = captures("[label](https://example.com)")
        local punct = {}
        for _, c in ipairs(caps) do
            if c.capture == "punctuation.bracket" then
                punct[#punct + 1] = c
            end
        end
        table.sort(punct, function(a, b) return a.sc < b.sc end)
        local texts = {}
        for _, c in ipairs(punct) do texts[#texts + 1] = c.text end
        assert.are.same({ "[", "]", "(", ")" }, texts)
    end)

    it("captures image punctuation as @punctuation.bracket", function()
        local caps = captures("![alt](https://example.com/img.png)")
        local punct = {}
        for _, c in ipairs(caps) do
            if c.capture == "punctuation.bracket" then
                punct[#punct + 1] = c
            end
        end
        table.sort(punct, function(a, b) return a.sc < b.sc end)
        local texts = {}
        for _, c in ipairs(punct) do texts[#texts + 1] = c.text end
        assert.are.same({ "!", "[", "]", "(", ")" }, texts)
    end)
end)
