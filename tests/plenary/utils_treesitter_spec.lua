---@diagnostic disable: undefined-global
local ts = require "utils/treesitter"

--- Scratch buffer of `lines` with `filetype`, its parser parsed when it has one.
--- @param lines string[]
--- @param filetype string
--- @return integer buf
local function open(lines, filetype)
    local buf = vim.api.nvim_create_buf(false, true)
    vim.api.nvim_buf_set_lines(buf, 0, -1, false, lines)
    vim.bo[buf].filetype = filetype
    local parser = vim.treesitter.get_parser(buf, nil, { error = false })
    if parser then parser:parse(true) end
    return buf
end

-- Language options come from vim.filetype.get_option, which sources ftplugins
-- only with `:filetype plugin on`. It caches the result per filetype, so the
-- setting must not change between tests.
vim.cmd "filetype plugin on"

describe("option_at", function()
    it("is the value of the language tree's filetype", function()
        local buf = open({ "local x = 1" }, "lua")
        assert.are.equal(":---,:--", ts.option_at(buf, { 0, 0 }, "comments"))
    end)

    it("skips an injected language that is not a filetype", function()
        local injections = ts.read_query("lua", "injections")
        vim.treesitter.query.set("lua", "injections",
            '((comment_content) @injection.content (#set! injection.language "markdown_inline"))')
        local buf = open({ "--- doc text" }, "lua")
        -- What a language that is not a filetype gets from get_option.
        local global_commentstring = vim.go.commentstring
        vim.go.commentstring = "#%s"
        local commentstring = ts.option_at(buf, { 0, 6 }, "commentstring")
        vim.go.commentstring = global_commentstring
        vim.treesitter.query.set("lua", "injections", injections)
        assert.are.equal("markdown_inline", ts.lang_tree_at(buf, { 0, 6 }):lang())
        assert.are.equal("-- %s", commentstring)
    end)

    it("is the buffer's value without a parser", function()
        local buf = open({ "x" }, "")
        vim.bo[buf].comments = ":;"
        assert.are.equal(":;", ts.option_at(buf, { 0, 0 }, "comments"))
    end)
end)

describe("commentstring_at", function()
    it("is the buffer option without a parser", function()
        local buf = open({ "x" }, "")
        vim.bo[buf].commentstring = "; %s"
        assert.are.equal("; %s", ts.commentstring_at(buf, { 0, 0 }))
    end)

    it("is the commentstring of the language tree's filetype", function()
        local buf = open({ "local x = 1" }, "lua")
        vim.bo[buf].commentstring = "; %s"
        assert.are.equal("-- %s", ts.commentstring_at(buf, { 0, 0 }))
    end)

    -- Last: replaces lua's highlights query for the rest of this nvim.
    it("prefers a capture's bo.commentstring metadata", function()
        vim.treesitter.query.set("lua", "highlights",
            '((identifier) @variable (#set! bo.commentstring "## %s"))')
        local buf = open({ "local x = 1" }, "lua")
        vim.treesitter.start(buf, "lua")
        assert.are.equal("## %s", ts.commentstring_at(buf, { 0, 6 }))
    end)
end)
