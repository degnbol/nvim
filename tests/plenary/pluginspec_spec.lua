---@diagnostic disable: undefined-global
-- Tests for lua/utils/pluginspec.lua

local pluginspec = require("utils.pluginspec")

-- Controlled fixture covering the spec forms in pack_specs.lua: gh() shorthand
-- slugs, full non-github URLs, and noise (version strings, comments).
local FIXTURE = {
    [[    gh("tpope/vim-repeat"),]],
    [[    gh("echasnovski/mini.nvim"),]],
    [[    { src = gh("L3MON4D3/LuaSnip"), version = vim.version.range("2") },]],
    [[    { src = gh("gert7/srt.nvim"), version = "main" },]],
    [[    { src = "https://gitlab.com/foo/bar" },]],
    [[    { src = "https://git.sr.ht/~someone/cool.nvim" },]],
    [[    -- a comment mentioning unrelated/text]],
}

describe("pluginspec.resolve", function()
    it("expands a github shorthand slug", function()
        assert.are.equal(
            "https://github.com/tpope/vim-repeat",
            pluginspec.resolve("vim-repeat", FIXTURE)
        )
    end)

    it("handles dotted repo names without misparsing the pattern", function()
        assert.are.equal(
            "https://github.com/echasnovski/mini.nvim",
            pluginspec.resolve("mini.nvim", FIXTURE)
        )
    end)

    it("resolves a slug inside a src = gh(...) table entry", function()
        assert.are.equal(
            "https://github.com/L3MON4D3/LuaSnip",
            pluginspec.resolve("LuaSnip", FIXTURE)
        )
    end)

    it("opens a full non-github URL verbatim", function()
        assert.are.equal(
            "https://gitlab.com/foo/bar",
            pluginspec.resolve("bar", FIXTURE)
        )
    end)

    it("honours a non-github host with a dotted repo", function()
        assert.are.equal(
            "https://git.sr.ht/~someone/cool.nvim",
            pluginspec.resolve("cool.nvim", FIXTURE)
        )
    end)

    it("matches only on the final path segment, not a prefix", function()
        -- "vim" is a prefix of "vim-repeat" but never a trailing segment here.
        assert.is_nil(pluginspec.resolve("vim", FIXTURE))
    end)

    it("returns nil for an unknown repo (caller falls through)", function()
        assert.is_nil(pluginspec.resolve("nonexistent", FIXTURE))
    end)
end)
