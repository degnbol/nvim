---@diagnostic disable: undefined-global
local tags = require "symbol-ls.tags"
local helpers = require "tests.helpers"

describe("symbol-ls.tags.definitions", function()
    before_each(function() helpers.fixture_tags_dir("r1") end)
    after_each(helpers.reset_tags_dir)
    local src = table.concat({
        "local function foo() end",
        "function M.foo(foo) end",
        "t.foo = 1",
        "M.foo()",
        "local function foobar() end",
    }, "\n")

    it("returns the function definitions only, at their @name", function()
        local hits = tags.definitions(src, "lua", "foo", "/repo/a.lua")
        local SymbolKind = vim.lsp.protocol.SymbolKind
        assert.are.same({ 0, 1 }, vim.tbl_map(function(h) return h.info.location.range.start.line end, hits))
        for _, hit in ipairs(hits) do
            assert.are.equal("foo", hit.info.name)
            assert.are.equal(SymbolKind.Function, hit.info.kind)
            assert.are.equal(vim.uri_from_fname("/repo/a.lua"), hit.info.location.uri)
            assert.is_nil(hit.client)
        end
        assert.are.same({ start = { line = 1, character = 11 }, ["end"] = { line = 1, character = 14 } },
            hits[2].info.location.range)
        assert.are.same({ 1, 0, 1, 23 }, hits[2].definition.range)
    end)

    it("starts a definition at the comments directly above it", function()
        local hits = tags.definitions("-- far\n\n-- near\n-- nearer\nfunction foo() end", "lua", "foo", "/a.lua")
        assert.are.equal(2, hits[1].definition.first_row)
    end)
end)

describe("symbol-ls.tags directives", function()
    it("strip! and select-adjacent! give the adjacent comments as @doc", function()
        local query = vim.treesitter.query.parse("lua", [[
            ((comment)* @doc
             .
             (function_declaration) @definition.function
             (#strip! @doc "^[\\s\\-]+|^\\s$")
             (#select-adjacent! @doc @definition.function))
        ]])
        local src = "-- far\n\n-- near\n--  nearer\nfunction foo() end"
        local root = vim.treesitter.get_string_parser(src, "lua"):parse()[1]:root()
        local docs = {}
        for _, match, metadata in query:iter_matches(root, src) do
            local id = assert(vim.iter(pairs(query.captures)):find(function(_, name) return name == "doc" end))
            assert.are.equal(2, #match[id])
            docs[#docs + 1] = metadata[id].text
        end
        assert.are.same({ "near\nnearer" }, docs)
    end)
end)

describe("symbol-ls.tags.ensure", function()
    local net_request = vim.net.request
    local notify = vim.notify
    local urls, notes, dir
    --- @type table<string, string|false> path suffix → body, false for 404
    local remote

    before_each(function()
        dir = helpers.fixture_tags_dir("r1")
        urls, notes, remote = {}, {}, {}
        vim.net.request = function(url, _, on_response)
            urls[#urls + 1] = url
            for suffix, body in pairs(remote) do
                if vim.endswith(url, suffix) then
                    if body then return on_response(nil, { body = body }) end
                    break
                end
            end
            on_response("curl: (22) The requested URL returned error: 404\n", nil)
        end
        vim.notify = function(msg) notes[#notes + 1] = msg end
    end)
    after_each(function()
        vim.net.request = net_request
        vim.notify = notify
        helpers.reset_tags_dir()
    end)

    --- @param revision string
    --- @param location string|nil
    local function source(revision, location)
        return function(lang)
            return { url = "https://github.com/owner/tree-sitter-" .. lang, revision = revision, location = location }
        end
    end

    it("is up to date with a matching sidecar", function()
        assert.is_true(helpers.await(function(done) tags.ensure("lua", source("r1"), done) end))
        assert.are.same({}, urls)
    end)

    it("downloads at a changed revision", function()
        remote["/r2/queries/tags.scm"] = "(identifier) @name"
        assert.is_true(helpers.await(function(done) tags.ensure("lua", source("r2"), done) end))
        assert.are.same({
            "https://raw.githubusercontent.com/owner/tree-sitter-lua/r2/tree-sitter.json",
            "https://raw.githubusercontent.com/owner/tree-sitter-lua/r2/queries/tags.scm",
        }, urls)
        assert.are.same({ "(identifier) @name" }, vim.fn.readfile(dir .. "/queries/lua/tags.scm"))
        assert.are.same({ lua = { revision = "r2", tags = true, inherits = {} } },
            vim.json.decode(table.concat(vim.fn.readfile(dir .. "/revisions.json"))))
        assert.are.same({}, notes)
    end)

    it("records a missing tags.scm and does not retry it", function()
        assert.is_false(helpers.await(function(done) tags.ensure("julia", source("r1"), done) end))
        local n_requests = #urls
        assert.is_false(helpers.await(function(done) tags.ensure("julia", source("r1"), done) end))
        assert.are.equal(n_requests, #urls)
        assert.are.same({}, notes)
    end)

    it("gives false with no source", function()
        assert.is_false(helpers.await(function(done) tags.ensure("lua", function() end, done) end))
    end)

    it("follows tree-sitter.json, inheriting node_modules grammars", function()
        remote["/ts/tree-sitter.json"] = vim.json.encode { grammars = {
            { name = "flow", path = "flow", tags = "queries/flow.scm" },
            { name = "ts", path = "ts", tags = { "queries/tags.scm", "node_modules/tree-sitter-base/queries/tags.scm" } },
        } }
        remote["/ts/queries/tags.scm"] = "(type_identifier) @name"
        remote["/b1/queries/tags.scm"] = "(identifier) @name"
        local sources = {
            ts = { url = "https://github.com/owner/tree-sitter-ts", revision = "ts", location = "ts" },
            base = { url = "https://github.com/owner/tree-sitter-base", revision = "b1" },
        }
        assert.is_true(helpers.await(function(done)
            tags.ensure("ts", function(lang) return sources[lang] end, done)
        end))
        assert.are.same({ "; inherits: base", "(type_identifier) @name" }, vim.fn.readfile(dir .. "/queries/ts/tags.scm"))
        assert.are.same({ "(identifier) @name" }, vim.fn.readfile(dir .. "/queries/base/tags.scm"))
    end)

    it("retries an inherited grammar that failed while its heir is up to date", function()
        remote["/ts/tree-sitter.json"] = vim.json.encode { grammars = {
            { name = "ts", tags = { "queries/tags.scm", "node_modules/tree-sitter-base/queries/tags.scm" } },
        } }
        remote["/ts/queries/tags.scm"] = "(type_identifier) @name"
        remote["/b1/tree-sitter.json"] = "<html>not json</html>"
        local sources = {
            ts = { url = "https://github.com/owner/tree-sitter-ts", revision = "ts" },
            base = { url = "https://github.com/owner/tree-sitter-base", revision = "b1" },
        }
        local function source_of(lang) return sources[lang] end
        assert.is_true(helpers.await(function(done) tags.ensure("ts", source_of, done) end))
        assert.are.equal(1, #notes)
        assert.is_nil(vim.uv.fs_stat(dir .. "/queries/base/tags.scm"))

        remote["/b1/tree-sitter.json"] = nil
        remote["/b1/queries/tags.scm"] = "(identifier) @name"
        assert.is_true(helpers.await(function(done) tags.ensure("ts", source_of, done) end))
        assert.are.same({ "(identifier) @name" }, vim.fn.readfile(dir .. "/queries/base/tags.scm"))
    end)

    it("downloads once for calls made while a download runs", function()
        remote["/r2/queries/tags.scm"] = "(identifier) @name"
        local request = vim.net.request
        vim.net.request = function(url, opts, on_response)
            vim.schedule(function() request(url, opts, on_response) end)
        end
        local results = {}
        for _ = 1, 2 do
            tags.ensure("lua", source("r2"), function(ok) results[#results + 1] = ok end)
        end
        assert(vim.wait(3000, function() return #results == 2 end))
        assert.are.same({ true, true }, results)
        assert.are.equal(2, #urls)
    end)

    it("warns on other download errors and records nothing", function()
        vim.net.request = function(_, _, on_response) on_response("curl: (6) Could not resolve host", nil) end
        assert.is_false(helpers.await(function(done) tags.ensure("lua", source("r3"), done) end))
        assert.are.equal(1, #notes)
        assert.are.equal("r1", vim.json.decode(table.concat(vim.fn.readfile(dir .. "/revisions.json"))).lua.revision)
    end)

    it("warns when the request cannot start, and a later call downloads", function()
        local request = vim.net.request
        vim.net.request = function() error("ENOENT: no such file or directory (cmd): 'curl'") end
        assert.is_false(helpers.await(function(done) tags.ensure("lua", source("r3"), done) end))
        assert.are.equal(1, #notes)
        vim.net.request = request
        remote["/r3/queries/tags.scm"] = "(identifier) @name"
        assert.is_true(helpers.await(function(done) tags.ensure("lua", source("r3"), done) end))
    end)
end)
