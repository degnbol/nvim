---@diagnostic disable: undefined-global
local symref = require "symref"
local lsp = require "symref.lsp_util"

--- @param lines string[]
--- @return integer bufnr
local function markdown_buf(lines)
    local bufnr = vim.api.nvim_create_buf(false, true)
    vim.api.nvim_buf_set_lines(bufnr, 0, -1, false, lines)
    vim.treesitter.start(bufnr, "markdown")
    return bufnr
end

describe("symref.code_span_at", function()
    local line = "see `ShellParse.extract_commands` here"
    --         0123456789...: span is bytes 4..32, backticks at 4 and 32
    local bufnr
    before_each(function()
        -- A second paragraph so the parse of row 0 alone is tested.
        bufnr = markdown_buf { line, "", "other `x` row" }
    end)

    it("returns the text inside a span", function()
        local text, range = symref.code_span_at(bufnr, 0, 10)
        assert.are.equal("ShellParse.extract_commands", text)
        assert.are.same({ 0, 4, 0, 33 }, range)
    end)

    it("returns the span on either backtick", function()
        assert.are.equal("ShellParse.extract_commands", (symref.code_span_at(bufnr, 0, 4)))
        assert.are.equal("ShellParse.extract_commands", (symref.code_span_at(bufnr, 0, 32)))
    end)

    it("returns nil outside any span", function()
        assert.is_nil(symref.code_span_at(bufnr, 0, 1))
        assert.is_nil(symref.code_span_at(bufnr, 0, 35))
    end)

    it("finds a span on a later row", function()
        assert.are.equal("x", (symref.code_span_at(bufnr, 2, 7)))
    end)
end)

describe("symref.parse_term", function()
    local cases = {
        ["ShellParse.extract_commands"] = { symbol = "extract_commands", qualifiers = { "ShellParse" } },
        ["Mod:method"] = { symbol = "method", qualifiers = { "Mod" } },
        ["a::b::c"] = { symbol = "c", qualifiers = { "a", "b" } },
        ["Class#method"] = { symbol = "method", qualifiers = { "Class" } },
        ["foo()"] = { symbol = "foo", qualifiers = {} },
        ["foo.bar(x, y)"] = { symbol = "bar", qualifiers = { "foo" } },
        ["foo"] = { symbol = "foo", qualifiers = {} },
    }
    for text, term in pairs(cases) do
        it("parses " .. text, function()
            assert.are.same(term, symref.parse_term(text))
        end)
    end

    for _, text in ipairs { "rg -n", "a b", "" } do
        it("rejects '" .. text .. "'", function()
            assert.is_nil(symref.parse_term(text))
        end)
    end
end)

describe("symref.lsp_util.convert_range", function()
    -- "é" and "ü" are 2 bytes and 1 UTF-16 unit each.
    local fname = vim.fn.tempname() .. ".txt"
    vim.fn.writefile({ "first", "é_ü foo" }, fname)
    local uri = vim.uri_from_fname(fname)
    local utf16 = { start = { line = 1, character = 4 }, ["end"] = { line = 1, character = 7 } }
    local utf8 = { start = { line = 1, character = 6 }, ["end"] = { line = 1, character = 9 } }

    it("converts utf-16 to utf-8 from disk", function()
        assert.are.same(utf8, lsp.convert_range(uri, utf16, "utf-16", "utf-8"))
    end)

    it("converts utf-8 to utf-16 from a loaded buffer", function()
        local bufnr = vim.fn.bufadd(fname)
        vim.fn.bufload(bufnr)
        -- The buffer, not the file, is read.
        vim.api.nvim_buf_set_lines(bufnr, 0, 1, false, {})
        local range = { start = { line = 0, character = 6 }, ["end"] = { line = 0, character = 9 } }
        local converted = lsp.convert_range(uri, range, "utf-8", "utf-16")
        vim.api.nvim_buf_delete(bufnr, { force = true })
        assert.are.same({ start = { line = 0, character = 4 }, ["end"] = { line = 0, character = 7 } }, converted)
    end)

    it("converts an unreadable file's lines as empty", function()
        local missing = vim.uri_from_fname(fname .. ".missing")
        assert.are.same(
            { start = { line = 1, character = 0 }, ["end"] = { line = 1, character = 0 } },
            lsp.convert_range(missing, utf16, "utf-16", "utf-8")
        )
    end)
end)

--- Start an in-process client whose replies `answer` produces. A nil from
--- `answer` means no reply ever.
--- @param name string
--- @param encoding string
--- @param answer fun(method: string, params: table): any
--- @return vim.lsp.Client
local function fake_client(name, encoding, answer)
    local id = vim.lsp.start({
        name = name,
        cmd = function(dispatchers)
            local n = 0
            return {
                request = function(method, params, callback)
                    n = n + 1
                    if method == "initialize" then
                        callback(nil, {
                            capabilities = { positionEncoding = encoding, workspaceSymbolProvider = true },
                        })
                    else
                        local reply = answer(method, params)
                        if reply ~= nil then reply(callback) end
                    end
                    return true, n
                end,
                notify = function() return true end,
                is_closing = function() return false end,
                terminate = function() dispatchers.on_exit(0, 15) end,
            }
        end,
    }, { attach = false })
    return assert(vim.lsp.get_client_by_id(assert(id)))
end

--- @param ... vim.lsp.Client
local function stop(...)
    for _, client in ipairs { ... } do
        client:stop(true)
    end
    vim.wait(1000, function() return #vim.lsp.get_clients() == 0 end)
end

describe("symref.lsp_util.request_all", function()
    local notify
    before_each(function()
        notify = vim.notify
        vim.notify = function() end
    end)
    after_each(function()
        vim.notify = notify
    end)

    it("calls back once, in client order, without excluded, failed or silent clients", function()
        local a = fake_client("a", "utf-16", function()
            return function(cb) vim.schedule(function() cb(nil, "a") end) end
        end)
        local silent = fake_client("silent", "utf-16", function() return nil end)
        local failing = fake_client("failing", "utf-16", function()
            return function(cb) cb({ code = 1, message = "boom" }) end
        end)
        local excluded = fake_client("excluded", "utf-16", function()
            return function(cb) cb(nil, "excluded") end
        end)
        local sync = fake_client("sync", "utf-16", function()
            return function(cb) cb(nil, "sync") end
        end)

        local calls = {}
        lsp.request_all("workspace/symbol", { query = "x" }, { exclude = { "excluded" }, timeout_ms = 50 },
            function(replies) calls[#calls + 1] = replies end)
        vim.wait(200, function() return false end)
        stop(a, silent, failing, excluded, sync)

        assert.are.equal(1, #calls)
        assert.are.same({ "a", "sync" }, vim.tbl_map(function(r) return r.result end, calls[1]))
        assert.are.same({ "a", "sync" }, vim.tbl_map(function(r) return r.client.name end, calls[1]))
    end)

    it("calls back at once with no eligible client", function()
        local replies
        lsp.request_all("workspace/symbol", { query = "x" }, { timeout_ms = 50 }, function(r) replies = r end)
        assert.are.same({}, replies)
    end)
end)

describe("symref.setup", function()
    local bufs = {}
    after_each(function()
        stop(unpack(vim.lsp.get_clients { name = "symref" }))
        vim.api.nvim_del_augroup_by_name("symref")
        for _, bufnr in ipairs(bufs) do
            vim.api.nvim_buf_delete(bufnr, { force = true })
        end
        bufs = {}
    end)

    --- @param filetype string
    --- @param buftype string
    --- @return integer bufnr
    local function buf_of(filetype, buftype)
        local bufnr = vim.api.nvim_create_buf(true, false)
        vim.bo[bufnr].buftype = buftype
        vim.bo[bufnr].filetype = filetype
        bufs[#bufs + 1] = bufnr
        return bufnr
    end

    --- @param bufnr integer
    --- @return boolean
    local function attached(bufnr)
        return #vim.lsp.get_clients { bufnr = bufnr, name = "symref" } > 0
    end

    it("attaches open and later buffers of the filetypes, but not markdown with a buftype", function()
        local open = buf_of("text", "")
        symref.setup { filetypes = { "markdown", "text" } }
        local later = buf_of("markdown", "")
        local float = buf_of("markdown", "nofile")
        local named = buf_of("text", "nofile")
        local other = buf_of("lua", "")
        assert.is_true(attached(open))
        assert.is_true(attached(later))
        assert.is_false(attached(float))
        assert.is_true(attached(named))
        assert.is_false(attached(other))
    end)
end)

describe("symref textDocument/definition", function()
    it("answers with the ranked hits in UTF-8 columns", function()
        local fname = vim.fn.tempname() .. "/shell_parse.lua"
        vim.fn.mkdir(vim.fs.dirname(fname), "p")
        -- "é" is 1 UTF-16 unit and 2 bytes, so `run` is at units 24-27, bytes 25-28.
        vim.fn.writefile({ "", "local é = 1; function M.run()" }, fname)
        local uri = vim.uri_from_fname(fname)
        local code = fake_client("code", "utf-16", function(method, params)
            assert.are.equal("workspace/symbol", method)
            assert.are.equal("run", params.query)
            return function(cb)
                cb(nil, {
                    { name = "run", kind = 12, location = { uri = "file:///elsewhere/other.lua",
                        range = { start = { line = 0, character = 0 }, ["end"] = { line = 0, character = 0 } } } },
                    { name = "run", kind = 12, location = { uri = uri,
                        range = { start = { line = 1, character = 24 }, ["end"] = { line = 1, character = 27 } } } },
                })
            end
        end)

        local bufnr = markdown_buf { "see `ShellParse.run()`" }
        vim.api.nvim_set_current_buf(bufnr)
        symref.attach(bufnr)
        local responses = vim.lsp.buf_request_sync(bufnr, "textDocument/definition", {
            textDocument = { uri = vim.uri_from_bufnr(bufnr) },
            position = { line = 0, character = 8 },
        }, 1000)
        stop(code, unpack(vim.lsp.get_clients { name = "symref" }))

        local _, response = next(assert(responses))
        assert.are.same({
            { uri = uri, range = { start = { line = 1, character = 25 }, ["end"] = { line = 1, character = 28 } } },
        }, assert(response).result)
    end)
end)

describe("symref.rank", function()
    --- @param name string
    --- @param path string
    --- @param container string|nil
    --- @return symref.Hit
    local function hit(name, path, container)
        return {
            info = {
                name = name,
                kind = 8,
                containerName = container,
                location = { uri = vim.uri_from_fname(path), range = {} },
            },
        }
    end
    --- @param hits symref.Hit[]
    --- @return string[]
    local function paths(hits)
        return vim.tbl_map(function(h) return vim.uri_to_fname(h.info.location.uri) end, hits)
    end
    local root = "/repo"

    it("puts the definition before test mocks", function()
        local term = { symbol = "_handle_input_submit_inner", qualifiers = {} }
        local hits = {}
        for _ = 1, 3 do
            table.insert(hits, hit(term.symbol, "/repo/lua/agentic/session_manager.test.lua"))
        end
        table.insert(hits, hit(term.symbol, "/repo/lua/agentic/session_manager.lua"))
        for _ = 1, 3 do
            table.insert(hits, hit(term.symbol, "/repo/lua/agentic/session_manager.test.lua"))
        end
        local ranked = symref.rank(hits, term, root)
        assert.are.equal(7, #ranked)
        assert.are.equal("/repo/lua/agentic/session_manager.lua", paths(ranked)[1])
    end)

    it("chooses the file stem the qualifier names", function()
        local term = { symbol = "extract_commands", qualifiers = { "ShellParse" } }
        local ranked = symref.rank({
            hit("extract_commands", "/repo/lua/other.lua"),
            hit("extract_commands", "/repo/lua/shell_parse.lua"),
        }, term, root)
        assert.are.same({ "/repo/lua/shell_parse.lua" }, paths(ranked))
    end)

    it("keeps every hit when no qualifier matches", function()
        local term = { symbol = "f", qualifiers = { "Nope" } }
        local ranked = symref.rank({ hit("f", "/repo/a.lua"), hit("f", "/repo/b.lua") }, term, root)
        assert.are.same({ "/repo/a.lua", "/repo/b.lua" }, paths(ranked))
    end)

    it("matches containerName", function()
        local term = { symbol = "run", qualifiers = { "my_class" } }
        local ranked = symref.rank({
            hit("run", "/repo/a.py", "Other"),
            hit("run", "/repo/b.py", "MyClass"),
        }, term, root)
        assert.are.same({ "/repo/b.py" }, paths(ranked))
    end)

    it("prefers the last qualifier", function()
        local term = { symbol = "f", qualifiers = { "a", "b" } }
        local ranked = symref.rank({ hit("f", "/repo/a.lua"), hit("f", "/repo/b.lua") }, term, root)
        assert.are.same({ "/repo/b.lua", "/repo/a.lua" }, paths(ranked))
    end)

    it("prefers files under root", function()
        local term = { symbol = "f", qualifiers = {} }
        local ranked = symref.rank({ hit("f", "/elsewhere/a.lua"), hit("f", "/repo/a.lua") }, term, root)
        assert.are.same({ "/repo/a.lua", "/elsewhere/a.lua" }, paths(ranked))
    end)

    it("drops fuzzy non-exact names", function()
        local term = { symbol = "ShellParse", qualifiers = {} }
        local ranked = symref.rank({
            hit("ShellParse", "/repo/a.lua"),
            hit("shell_parse_x", "/repo/b.lua"),
            hit("SHellParser", "/repo/c.lua"),
        }, term, root)
        assert.are.same({ "/repo/a.lua" }, paths(ranked))
    end)
end)
