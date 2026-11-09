local lsp_util = require "symbol-ls.lsp_util"
local tags = require "symbol-ls.tags"

--- In-process LSP server that resolves code symbols named in markdown code
--- spans (`` `Foo.bar` ``) through the `workspace/symbol` of the other running
--- servers, else through tree-sitter `tags` queries on the git-tracked files
--- that contain the symbol. Serves definition, hover and workspace symbol.
local M = {}

local NAME = "symbol-ls"
-- Byte caps on files read per lookup, since parsing runs on the main thread.
local READ_FILE_MAX = 1e6
local READ_TOTAL_MAX = 5e6

--- @class symbol_ls.Hit
--- @field info lsp.SymbolInformation range in `client`'s position encoding, else UTF-8
--- @field client? vim.lsp.Client nil for a hit from a tags query
--- @field definition? symbol_ls.Definition set for a hit from a tags query

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

--- Files tracked by git under `root`, submodules included, that contain
--- `word` as a whole word. Binary files are skipped.
--- @param root string
--- @param word string
--- @param on_done fun(paths: string[]) absolute. `{}` outside a git repository or on failure (notified at WARN)
function M.files_with_word(root, word, on_done)
    if not vim.fs.root(root, ".git") then return on_done({}) end
    local cmd = { "git", "grep", "-z", "-l", "-I", "-w", "-F", "--recurse-submodules", "-e", word }
    vim.system(cmd, { cwd = root, text = true }, vim.schedule_wrap(function(result)
        -- Exit code 1: no match.
        if result.code == 1 and result.stderr == "" then return on_done({}) end
        if result.code ~= 0 then
            vim.notify(("git grep (%s): %s"):format(root, result.stderr), vim.log.levels.WARN)
            return on_done({})
        end
        local paths = vim.split(result.stdout, "\0", { trimempty = true })
        on_done(vim.tbl_map(function(path) return vim.fs.joinpath(root, path) end, paths))
    end))
end

--- @class symbol_ls.TaggedFile
--- @field path string absolute
--- @field lang string tree-sitter language

--- Languages whose tags query raised an error, skipped until restart.
--- @type table<string, true>
local failed_langs = {}

--- `files` with those whose stem matches one of `qualifiers` under
--- `fold_name` first. Order is kept otherwise.
--- @param files symbol_ls.TaggedFile[]
--- @param qualifiers string[]
--- @return symbol_ls.TaggedFile[]
local function qualified_first(files, qualifiers)
    local folded = vim.tbl_map(fold_name, qualifiers)
    local matching, rest = {}, {}
    for _, file in ipairs(files) do
        local stem = fold_name(vim.fn.fnamemodify(file.path, ":t:r"))
        table.insert(vim.list_contains(folded, stem) and matching or rest, file)
    end
    return vim.list_extend(matching, rest)
end

--- The git-tracked files under `root` that contain `word` as a word and whose
--- language has a parser and a tree-sitter `tags` query. Missing tags queries
--- are downloaded from the grammar repositories in nvim-treesitter's registry.
--- @param root string
--- @param word string
--- @param on_done fun(files: symbol_ls.TaggedFile[]) in `files_with_word` order
function M.tagged_files(root, word, on_done)
    M.files_with_word(root, word, function(paths)
        local files = {} --- @type symbol_ls.TaggedFile[]
        local langs = {} --- @type table<string, true>
        for _, path in ipairs(paths) do
            local filetype = vim.filetype.match { filename = path }
            local lang = filetype and vim.treesitter.language.get_lang(filetype)
            if lang and not failed_langs[lang] and vim.treesitter.language.add(lang) then
                files[#files + 1] = { path = path, lang = lang }
                langs[lang] = true
            end
        end
        local with_query = {} --- @type table<string, boolean>
        local function finish()
            on_done(vim.tbl_filter(function(file) return with_query[file.lang] end, files))
        end
        local n_pending = vim.tbl_count(langs)
        if n_pending == 0 then return finish() end
        for lang in pairs(langs) do
            tags.ensure(lang, tags.registry_source, function(ok)
                with_query[lang] = ok
                n_pending = n_pending - 1
                if n_pending == 0 then finish() end
            end)
        end
    end)
end

--- Definitions of `term.symbol` in `files` by their tags queries. Files whose
--- stem matches a qualifier of `term` are read first. Files over
--- `READ_FILE_MAX` bytes are skipped, and reading stops at `READ_TOTAL_MAX`
--- bytes in total. A language whose query raises an error is reported at WARN
--- and skipped until restart.
--- @param files symbol_ls.TaggedFile[]
--- @param term symbol_ls.Term
--- @return symbol_ls.Hit[] unranked
function M.tag_hits(files, term)
    local hits, n_bytes = {}, 0
    for _, file in ipairs(qualified_first(files, term.qualifiers)) do
        local stat = vim.uv.fs_stat(file.path)
        if stat and stat.size <= READ_FILE_MAX and not failed_langs[file.lang] then
            n_bytes = n_bytes + stat.size
            if n_bytes > READ_TOTAL_MAX then break end
            local ok, result = pcall(tags.definitions, vim.fn.readblob(file.path), file.lang, term.symbol, file.path)
            if ok then
                vim.list_extend(hits, result)
            else
                failed_langs[file.lang] = true
                vim.notify(("tags query for %s: %s"):format(file.lang, result), vim.log.levels.WARN)
            end
        end
    end
    return hits
end

--- Markdown hover for a hit from a tags query: its path relative to `root`,
--- its `@doc` text, and in a code block the comments directly above the
--- definition and the definition's first line.
--- @param hit symbol_ls.Hit
--- @param root string|nil nil, or a root not above the file, shows the absolute path
--- @return lsp.Hover
function M.snippet(hit, root)
    local path = vim.uri_to_fname(hit.info.location.uri)
    local definition = assert(hit.definition)
    local rows = vim.list_slice(vim.fn.readfile(path, "", definition.range[1] + 1), definition.first_row + 1)
    local lines = { ("`%s`"):format(root and vim.fs.relpath(root, path) or path) }
    if definition.doc then vim.list_extend(lines, { "", definition.doc, "" }) end
    table.insert(lines, "```" .. (vim.filetype.match { filename = path } or ""))
    vim.list_extend(lines, rows)
    table.insert(lines, "```")
    return { contents = { kind = "markdown", value = table.concat(lines, "\n") } }
end

--- Hover for `hit`. For a hit from a tags query, the definition's source and
--- doc text. Else textDocument/hover forwarded to the client that reported
--- it, at the range start in that client's encoding. No buffer is loaded or
--- attached (lua_ls answers for unopened files). nil when the client has no
--- hover or answers nothing.
--- @param hit symbol_ls.Hit
--- @param root string|nil paths of tags query hits are shown relative to it
--- @param on_done fun(result: lsp.Hover|nil)
function M.hover_at(hit, root, on_done)
    local method = "textDocument/hover"
    local client, location = hit.client, hit.info.location
    if not client then return on_done(M.snippet(hit, root)) end
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
    if not hit.client then return location end
    local range = lsp_util.convert_range(location.uri, location.range, hit.client.offset_encoding, "utf-8")
    return { uri = location.uri, range = range }
end

--- Ranked hits for the code span at a text document position: those of the
--- running servers, else those of tags queries. The tagged files are found
--- while the servers answer, but parsed only when no server has a hit.
--- @param params lsp.TextDocumentPositionParams
--- @param on_done fun(hits: symbol_ls.Hit[], span: Range4|nil, root: string|nil)
local function span_hits(params, on_done)
    local uri = params.textDocument.uri
    local current = vim.api.nvim_get_current_buf()
    -- Unnamed buffers share the URI "file://", which uri_to_bufnr cannot map back.
    local bufnr = vim.uri_from_bufnr(current) == uri and current or vim.uri_to_bufnr(uri)
    local text, span = M.code_span_at(bufnr, params.position.line, params.position.character)
    local term = text and M.parse_term(text)
    if not term then return on_done({}, nil, nil) end
    local root = M.root(bufnr)
    local server_hits, files
    local function finish()
        if not server_hits or not files then return end
        local ranked = M.rank(server_hits, term, root)
        if #ranked == 0 then ranked = M.rank(M.tag_hits(files, term), term, root) end
        on_done(ranked, span, root)
    end
    if root then
        M.tagged_files(root, term.symbol, function(tagged)
            files = tagged
            finish()
        end)
    else
        files = {}
    end
    M.query(term.symbol, function(hits)
        server_hits = hits
        finish()
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
        span_hits(params, function(hits, span, root)
            if not hits[1] or not span then return reply(nil) end
            M.hover_at(hits[1], root, function(result)
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
--- markdown buffers with a 'buftype'. Add `tags.dir()` to 'runtimepath'.
--- @param opts { filetypes: string[]|nil }|nil filetypes default `{ "markdown" }`
function M.setup(opts)
    local filetypes = opts and opts.filetypes or { "markdown" }
    vim.opt.runtimepath:append(tags.dir())
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
