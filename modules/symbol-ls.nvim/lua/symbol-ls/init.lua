local lsp_util = require "symbol-ls.lsp_util"

--- In-process LSP server that resolves code symbols named in markdown code
--- spans (`` `Foo.bar` ``) through the `workspace/symbol`
--- of the other running servers. Serves definition, hover and workspace symbol.
local M = {}

local NAME = "symbol-ls"

--- @class symbol_ls.Hit
--- @field info lsp.SymbolInformation range in `client`'s position encoding
--- @field client vim.lsp.Client

--- @class symbol_ls.Term
--- @field symbol string
--- @field qualifiers string[]

--- Text of the markdown_inline `code_span` at a position, without its
--- backticks. Finds spans in injected trees too. nil when no span covers the
--- position.
--- @param bufnr integer
--- @param row integer 0-based
--- @param col integer 0-based byte
--- @return string|nil text
--- @return Range4|nil range of the span
function M.code_span_at(bufnr, row, col)
    local parser = vim.treesitter.get_parser(bufnr)
    if not parser then return nil end
    -- The highlighter parses injections only for drawn ranges.
    parser:parse { row, row + 1 }
    local node = parser:named_node_for_range({ row, col, row, col }, { ignore_injections = false })
    while node and node:type() ~= "code_span" do
        node = node:parent()
    end
    if not node then return nil end
    local text = vim.treesitter.get_node_text(node, bufnr):gsub("^`+", ""):gsub("`+$", "")
    return vim.trim(text), { node:range() }
end

--- Split a term into a qualifier path and a symbol. Separators `.`, `:`, `#`,
--- `::`; trailing `()` and call arguments dropped (`foo.bar()` → bar).
--- Returns nil for text that is not an identifier path (`rg -n`, `a b`).
--- @param text string
--- @return symbol_ls.Term|nil
function M.parse_term(text)
    local path = text:gsub("%b()$", "")
    local dotted = path:gsub("::", "."):gsub("[:#]", ".")
    local segments = vim.split(dotted, ".", { plain = true })
    for _, segment in ipairs(segments) do
        if not segment:match("^[%a_][%w_]*$") then return nil end
    end
    local symbol = table.remove(segments)
    return { symbol = symbol, qualifiers = segments }
end

--- Lower-case `name` with `_` and `-` removed, so `ShellParse`, `shell_parse`
--- and `shell-parse` compare equal.
--- @param name string
--- @return string
local function fold_name(name)
    return (name:lower():gsub("[_%-]", ""))
end

local TEST_NAMES = { "%.test%.", "_spec%.", "_test%.", "^test_" }

--- Whether `path` names a test file, by its basename or a `tests/` or `spec/`
--- directory.
--- @param path string
--- @return boolean
local function is_test_file(path)
    local name = vim.fs.basename(path)
    for _, pattern in ipairs(TEST_NAMES) do
        if name:find(pattern) then return true end
    end
    return path:find("/tests/", 1, true) ~= nil or path:find("/spec/", 1, true) ~= nil
end

--- Index of the last of `qualifiers` that matches the symbol's containerName
--- or the stem of its file, case-folded.
--- @param info lsp.SymbolInformation
--- @param path string
--- @param qualifiers string[]
--- @return integer 0 when none matches
local function last_matching_qualifier(info, path, qualifiers)
    local stem = fold_name(vim.fn.fnamemodify(path, ":t:r"))
    local container = info.containerName and fold_name(info.containerName)
    for i = #qualifiers, 1, -1 do
        local qualifier = fold_name(qualifiers[i])
        if qualifier == stem or qualifier == container then return i end
    end
    return 0
end

--- Rank exact-name hits. Tiers, best first:
---   1. a qualifier matches containerName or the file stem (case-folded,
---      `_`/`-` removed; the last qualifier is tried first),
---   2. file under `root`,
---   3. not a test file (`*.test.*`, `*_spec.*`, `*_test.*`, `test_*`,
---      `tests/`, `spec/`).
--- Ties keep input order. Hits whose `name` differs from `symbol` are dropped,
--- and with qualifiers present, hits that match none are dropped only when
--- some hit matches one.
--- @param hits symbol_ls.Hit[]
--- @param term symbol_ls.Term
--- @param root string|nil
--- @return symbol_ls.Hit[]
function M.rank(hits, term, root)
    local entries = {}
    for i, hit in ipairs(hits) do
        if hit.info.name == term.symbol then
            local path = vim.uri_to_fname(hit.info.location.uri)
            entries[#entries + 1] = {
                hit = hit,
                i = i,
                qualifier = last_matching_qualifier(hit.info, path, term.qualifiers),
                under_root = root ~= nil and vim.fs.relpath(root, path) ~= nil,
                test = is_test_file(path),
            }
        end
    end
    if vim.iter(entries):any(function(e) return e.qualifier > 0 end) then
        entries = vim.tbl_filter(function(e) return e.qualifier > 0 end, entries)
    end
    table.sort(entries, function(a, b)
        if a.qualifier ~= b.qualifier then return a.qualifier > b.qualifier end
        if a.under_root ~= b.under_root then return a.under_root end
        if a.test ~= b.test then return b.test end
        return a.i < b.i
    end)
    return vim.tbl_map(function(e) return e.hit end, entries)
end

--- workspace/symbol for `symbol` from every running server except symbol-ls,
--- each hit tagged with its client. Waits at most 800 ms. Drops hits that
--- have no range.
--- @param symbol string
--- @param on_done fun(hits: symbol_ls.Hit[])
function M.query(symbol, on_done)
    -- Below the 1000 ms that 'tagfunc' waits for this server.
    local opts = { exclude = { NAME }, timeout_ms = 800 }
    lsp_util.request_all("workspace/symbol", { query = symbol }, opts, function(replies)
        --- @type symbol_ls.Hit[]
        local hits = {}
        for _, reply in ipairs(replies) do
            for _, info in ipairs(reply.result or {}) do
                -- A WorkspaceSymbol may hold a location without a range,
                -- which only workspaceSymbol/resolve fills in.
                if info.location.range then
                    hits[#hits + 1] = { info = info, client = reply.client }
                end
            end
        end
        on_done(hits)
    end)
end

--- Project root of `bufnr`: the `.git` root for a file buffer, else the cwd
--- of a window showing the buffer, or of the current window.
--- @param bufnr integer
--- @return string|nil nil for a file buffer outside a git repository
function M.root(bufnr)
    if vim.bo[bufnr].buftype == "" and vim.api.nvim_buf_get_name(bufnr) ~= "" then
        return vim.fs.root(bufnr, ".git")
    end
    return vim.fn.getcwd(vim.fn.win_findbuf(bufnr)[1] or 0)
end

--- Forward textDocument/hover to the client that reported `hit`, at the
--- range start in that client's encoding. No buffer is loaded or attached
--- (lua_ls answers for unopened files). nil when the client has no hover or
--- answers nothing.
--- @param hit symbol_ls.Hit
--- @param on_done fun(result: lsp.Hover|nil)
function M.hover_at(hit, on_done)
    local method = "textDocument/hover"
    local client, location = hit.client, hit.info.location
    if not client:supports_method(method) then return on_done(nil) end
    local params = { textDocument = { uri = location.uri }, position = location.range.start }
    local sent = client:request(method, params, function(err, result)
        if err then
            vim.notify(("%s (%s): %s"):format(method, client.name, err.message), vim.log.levels.WARN)
        end
        on_done(result)
    end)
    if not sent then
        vim.notify(("%s (%s): request not sent"):format(method, client.name), vim.log.levels.WARN)
        on_done(nil)
    end
end

--- Location of a hit in UTF-8 columns.
--- @param hit symbol_ls.Hit
--- @return lsp.Location
local function utf8_location(hit)
    local location = hit.info.location
    local range = lsp_util.convert_range(location.uri, location.range, hit.client.offset_encoding, "utf-8")
    return { uri = location.uri, range = range }
end

--- Ranked hits for the code span at a text document position.
--- @param params lsp.TextDocumentPositionParams
--- @param on_done fun(hits: symbol_ls.Hit[], span: Range4|nil)
local function span_hits(params, on_done)
    local uri = params.textDocument.uri
    local current = vim.api.nvim_get_current_buf()
    -- Unnamed buffers share the URI "file://", which uri_to_bufnr cannot map back.
    local bufnr = vim.uri_from_bufnr(current) == uri and current or vim.uri_to_bufnr(uri)
    local text, span = M.code_span_at(bufnr, params.position.line, params.position.character)
    local term = text and M.parse_term(text)
    if not term then return on_done({}, nil) end
    local root = M.root(bufnr)
    M.query(term.symbol, function(hits)
        on_done(M.rank(hits, term, root), span)
    end)
end

--- Request handlers, each answering through `reply`.
--- @type table<string, fun(params: table, reply: fun(result: any))>
local HANDLERS = {
    initialize = function(_, reply)
        reply {
            capabilities = {
                -- The request position is a byte column, as code_span_at takes.
                -- No textDocumentSync: the handlers read the buffer directly.
                positionEncoding = "utf-8",
                definitionProvider = true,
                hoverProvider = true,
                workspaceSymbolProvider = true,
            },
        }
    end,
    -- All ranked hits: `<C-]>` jumps to the first, `g]` lists them. nil lets 'tagfunc' fall back to the tags file.
    ["textDocument/definition"] = function(params, reply)
        span_hits(params, function(hits)
            reply(#hits > 0 and vim.tbl_map(utf8_location, hits) or nil)
        end)
    end,
    -- The reply's range is the code span, which vim.lsp.buf.hover highlights
    -- in the current buffer.
    ["textDocument/hover"] = function(params, reply)
        span_hits(params, function(hits, span)
            if not hits[1] or not span then return reply(nil) end
            M.hover_at(hits[1], function(result)
                if result then
                    result.range = {
                        start = { line = span[1], character = span[2] },
                        ["end"] = { line = span[3], character = span[4] },
                    }
                end
                reply(result)
            end)
        end)
    end,
    -- Live-search pickers send empty and partial queries; those get `[]`.
    ["workspace/symbol"] = function(params, reply)
        local term = M.parse_term(params.query)
        if not term then return reply({}) end
        local root = M.root(vim.api.nvim_get_current_buf())
        M.query(term.symbol, function(hits)
            local ranked = M.rank(hits, term, root)
            reply(vim.tbl_map(function(hit)
                return vim.tbl_extend("force", hit.info, { location = utf8_location(hit) })
            end, ranked))
        end)
    end,
}

--- In-process server (`:h vim.lsp.rpc.PublicClient`) for definition, hover
--- and workspace symbol.
--- @param dispatchers vim.lsp.rpc.Dispatchers
--- @return vim.lsp.rpc.PublicClient
function M._rpc(dispatchers)
    local request_id = 0
    local closing = false
    return {
        request = function(method, params, callback, on_reply)
            request_id = request_id + 1
            local id = request_id
            local handler = HANDLERS[method] or function(_, reply) reply(nil) end
            handler(params, function(result)
                callback(nil, result, id)
                if on_reply then on_reply(id) end
            end)
            return true, id
        end,
        notify = function(method)
            if method == "exit" then
                closing = true
                dispatchers.on_exit(0, 0)
            end
            return true
        end,
        is_closing = function() return closing end,
        terminate = function()
            closing = true
            dispatchers.on_exit(0, 15)
        end,
    }
end

--- Attach `bufnr` to the symbol-ls client, starting it on first use. One client
--- serves every buffer.
--- @param bufnr integer
function M.attach(bufnr)
    -- No root_dir, so reuse_client_default reuses the client for any buffer.
    vim.lsp.start({ name = NAME, cmd = M._rpc }, { bufnr = bufnr })
end

--- Whether `attach` should run on a buffer of a configured filetype.
--- @param bufnr integer
--- @return boolean
local function is_attachable(bufnr)
    -- Hover and preview floats, and many plugins' scratch buffers, are
    -- markdown with a 'buftype'. Other filetypes are named on purpose.
    return vim.bo[bufnr].filetype ~= "markdown" or vim.bo[bufnr].buftype == ""
end

--- Attach every buffer of `opts.filetypes`, open now or later, except
--- markdown buffers with a 'buftype'.
--- @param opts { filetypes: string[]|nil }|nil filetypes default `{ "markdown" }`
function M.setup(opts)
    local filetypes = opts and opts.filetypes or { "markdown" }
    vim.api.nvim_create_autocmd("FileType", {
        group = vim.api.nvim_create_augroup(NAME, { clear = true }),
        pattern = filetypes,
        callback = function(ev)
            if is_attachable(ev.buf) then M.attach(ev.buf) end
        end,
    })
    for _, bufnr in ipairs(vim.api.nvim_list_bufs()) do
        if vim.api.nvim_buf_is_loaded(bufnr)
            and vim.list_contains(filetypes, vim.bo[bufnr].filetype)
            and is_attachable(bufnr) then
            M.attach(bufnr)
        end
    end
end

return M
