--- Tree-sitter `tags` queries (https://tree-sitter.github.io/tree-sitter/4-code-navigation.html):
--- downloading them from grammar repositories and finding definitions with them.
local M = {}

local SymbolKind = vim.lsp.protocol.SymbolKind

--- Symbol kind of each `@definition.<kind>` suffix. Others are Variable.
--- @type table<string, lsp.SymbolKind>
local KINDS = {
    class = SymbolKind.Class,
    constant = SymbolKind.Constant,
    ["function"] = SymbolKind.Function,
    interface = SymbolKind.Interface,
    macro = SymbolKind.Function,
    method = SymbolKind.Method,
    module = SymbolKind.Module,
    type = SymbolKind.Class,
}

--- @class symbol_ls.Source
--- @field url string repository
--- @field revision string commit or tag
--- @field location? string directory of the grammar in the repository

--- @class symbol_ls.Definition
--- @field range Range4 the `@definition.<kind>` node
--- @field first_row integer 0-based first row of the comments directly above the node, else its start row
--- @field doc? string `@doc` text

--- Directory holding downloaded queries as `queries/<lang>/tags.scm` and the
--- revision each was downloaded at in `revisions.json`.
--- @return string
function M.dir()
    return vim.fn.stdpath("data") .. "/symbol-ls"
end

--- Repository of `lang`'s parser in nvim-treesitter's registry.
--- @param lang string
--- @return symbol_ls.Source|nil nil without nvim-treesitter or a registry entry
function M.registry_source(lang)
    if not vim.api.nvim_get_runtime_file("lua/nvim-treesitter/parsers.lua", false)[1] then return nil end
    local entry = require("nvim-treesitter.parsers")[lang]
    return entry and entry.install_info
end

--- @param lang string
--- @return string
local function query_path(lang)
    return M.dir() .. "/queries/" .. lang .. "/tags.scm"
end

--- @class symbol_ls.Revision
--- @field revision string
--- @field tags boolean whether the repository has a tags query
--- @field inherits? string[] languages the query inherits from

--- `revisions.json`. One that does not decode is notified at WARN and read as
--- empty.
--- @return table<string, symbol_ls.Revision> by language
local function read_revisions()
    local path = M.dir() .. "/revisions.json"
    if not vim.uv.fs_stat(path) then return {} end
    local ok, revisions = pcall(vim.json.decode, vim.fn.readblob(path))
    if ok then return revisions end
    vim.notify(("%s: %s"):format(path, revisions), vim.log.levels.WARN)
    return {}
end

--- @param lang string
--- @param revision symbol_ls.Revision
local function write_revision(lang, revision)
    local revisions = read_revisions()
    revisions[lang] = revision
    vim.fn.mkdir(M.dir(), "p")
    vim.fn.writefile({ vim.json.encode(revisions) }, M.dir() .. "/revisions.json")
end

--- URL of raw files in a GitHub repository at a revision.
--- @param source symbol_ls.Source
--- @return string|nil nil for a repository not on GitHub
local function raw_base(source)
    local repo = source.url:gsub("%.git$", ""):match("^https://github%.com/([^/]+/[^/]+)/?$")
    return repo and ("https://raw.githubusercontent.com/%s/%s"):format(repo, source.revision)
end

--- GET `url`, calling back on the main loop.
--- @param url string
--- @param on_done fun(body: string|nil, err: string|nil) both nil for a 404
local function fetch(url, on_done)
    -- vim.system raises when curl is missing.
    local sent, send_err = pcall(vim.net.request, url, nil, function(err, response)
        vim.schedule(function()
            if err and err:find("returned error: 404", 1, true) then return on_done(nil, nil) end
            on_done(response and response.body, err)
        end)
    end)
    if not sent then vim.schedule(function() on_done(nil, send_err) end) end
end

--- Tags query files of the grammar at `location` per a `tree-sitter.json`.
--- @param json string|nil content of `tree-sitter.json`, nil for a repository without one
--- @param location string|nil
--- @return string[]|nil paths relative to the repository root, nil when `json` does not decode
--- @return string|nil err
local function tags_paths(json, location)
    local default = { "queries/tags.scm" }
    if not json then return default end
    local ok, config = pcall(vim.json.decode, json)
    if not ok then return nil, "tree-sitter.json: " .. config end
    for _, grammar in ipairs(config.grammars or {}) do
        if (grammar.path or ".") == (location or ".") then
            local paths = grammar.tags or default
            return type(paths) == "string" and { paths } or paths
        end
    end
    return default
end

--- Tags query of the grammar at `location` in the repository at `base`: its
--- own query files concatenated, after an `; inherits:` line naming the
--- grammars whose query files it takes from `node_modules/tree-sitter-<lang>/`.
--- @param base string URL of the repository's raw files at a revision
--- @param location string|nil
--- @param on_done fun(text: string|nil, inherited: string[], err: string|nil) text nil when there is no query
local function download(base, location, on_done)
    fetch(base .. "/tree-sitter.json", function(json, json_err)
        if json_err then return on_done(nil, {}, json_err) end
        local paths, paths_err = tags_paths(json, location)
        if not paths then return on_done(nil, {}, paths_err) end
        local inherited, own = {}, {}
        for _, path in ipairs(paths) do
            local lang = path:match("^node_modules/tree%-sitter%-([%w_%-]+)/")
            if lang then
                inherited[#inherited + 1] = lang:gsub("%-", "_")
            else
                own[#own + 1] = path
            end
        end
        local lines = #inherited > 0 and { "; inherits: " .. table.concat(inherited, ",") } or {}
        local function fetch_own(i)
            if i > #own then
                return on_done(#lines > 0 and table.concat(lines, "\n") or nil, inherited, nil)
            end
            fetch(base .. "/" .. own[i], function(body, err)
                if err then return on_done(nil, {}, err) end
                lines[#lines + 1] = body
                fetch_own(i + 1)
            end)
        end
        fetch_own(1)
    end)
end

--- Store `lang`'s tags query and its revision.
--- @param lang string
--- @param text string|nil the query, nil for none
--- @param revision symbol_ls.Revision
local function store(lang, text, revision)
    if text then
        vim.fn.mkdir(vim.fs.dirname(query_path(lang)), "p")
        vim.fn.writefile(vim.split(text, "\n"), query_path(lang))
        -- query.get memoizes, nil included, until its OptionSet autocmd clears it.
        vim.o.runtimepath = vim.o.runtimepath
    end
    write_revision(lang, revision)
end

--- Callbacks of the `ensure` calls waiting on a download, by language.
--- @type table<string, fun(ok: boolean)[]>
local waiting = {}

--- `ensure` each of `langs`, then call back.
--- @param langs string[]
--- @param source_of fun(lang: string): symbol_ls.Source|nil
--- @param on_done fun()
local function ensure_all(langs, source_of, on_done)
    local n_pending = #langs
    if n_pending == 0 then return on_done() end
    for _, lang in ipairs(langs) do
        M.ensure(lang, source_of, function()
            n_pending = n_pending - 1
            if n_pending == 0 then on_done() end
        end)
    end
end

--- Make `lang`'s tags query available to `vim.treesitter.query.get` from
--- `dir()`, downloading it at the revision `source_of(lang)` names unless
--- `revisions.json` records that revision. A grammar repository without a tags
--- query is recorded too. Grammars the query inherits from are ensured with
--- the same `source_of`, also when `lang` is up to date. Download and write
--- errors are notified at WARN and recorded nowhere.
--- @param lang string
--- @param source_of fun(lang: string): symbol_ls.Source|nil
--- @param on_done fun(ok: boolean) whether `lang` has a tags query; false without a source on GitHub
function M.ensure(lang, source_of, on_done)
    local source = source_of(lang)
    local base = source and raw_base(source)
    if not source or not base then return on_done(false) end
    local known = read_revisions()[lang]
    if known and known.revision == source.revision and (not known.tags or vim.uv.fs_stat(query_path(lang))) then
        return ensure_all(known.inherits or {}, source_of, function() on_done(known.tags) end)
    end
    if waiting[lang] then
        table.insert(waiting[lang], on_done)
        return
    end
    waiting[lang] = { on_done }
    local function finish(ok)
        local callbacks = waiting[lang]
        waiting[lang] = nil
        for _, callback in ipairs(callbacks) do callback(ok) end
    end

    download(base, source.location, function(text, inherited, err)
        local ok = err == nil
        if ok then
            ok, err = pcall(store, lang, text, { revision = source.revision, tags = text ~= nil, inherits = inherited })
        end
        if not ok then
            vim.notify(("tags query for %s: %s"):format(lang, err), vim.log.levels.WARN)
            return finish(false)
        end
        ensure_all(inherited, source_of, function() finish(text ~= nil) end)
    end)
end

--- Vim regex matching what the Rust regex `pattern` matches, for the
--- constructs tags queries use: alternation, `+`, anchors, and `\s` in and
--- outside `[]`.
--- @param pattern string
--- @return string
local function vim_regex(pattern)
    local sets_expanded = pattern:gsub("%b[]", function(set) return (set:gsub("\\s", " \\t\\n\\r")) end)
    return "\\v" .. sets_expanded:gsub("\\s", "\\_s")
end

--- Store the text of the nodes `match[id]`, joined by newlines and each with
--- `metadata[id].strip` removed, as `metadata[id].text`.
--- @param match table<integer, TSNode[]>
--- @param id integer
--- @param source integer|string
--- @param metadata vim.treesitter.query.TSMetadata
local function set_doc_text(match, id, source, metadata)
    metadata[id] = metadata[id] or {}
    local strip = metadata[id].strip --[[@as string|nil]]
    metadata[id].text = table.concat(vim.tbl_map(function(node)
        local text = vim.treesitter.get_node_text(node, source)
        return strip and vim.fn.substitute(text, strip, "", "g") or text
    end, match[id] or {}), "\n")
end

-- `(#strip! @doc "regex")`: `@doc`'s text without the matches of a Rust regex.
vim.treesitter.query.add_directive("strip!", function(match, _, source, predicate, metadata)
    local id = predicate[2]
    metadata[id] = metadata[id] or {}
    metadata[id].strip = vim_regex(predicate[3])
    set_doc_text(match, id, source, metadata)
end, {})

--- `(#select-adjacent! @doc @definition)`: keep only the `@doc` nodes in an
--- unbroken run of rows ending on the row before `@definition`.
--- @param match table<integer, TSNode[]>
--- @param source integer|string
--- @param predicate any[]
--- @param metadata vim.treesitter.query.TSMetadata
local function select_adjacent(match, _, source, predicate, metadata)
    local id, adjacent = predicate[2], match[predicate[3]]
    local nodes = match[id] or {}
    if adjacent then
        local first = #nodes + 1
        local row = adjacent[1]:start()
        while first > 1 and nodes[first - 1]:end_() + 1 >= row do
            first = first - 1
            row = nodes[first]:start()
        end
        match[id] = vim.list_slice(nodes, first)
    end
    set_doc_text(match, id, source, metadata)
end
vim.treesitter.query.add_directive("select-adjacent!", select_adjacent, {})
-- tree-sitter-go spells it this way.
vim.treesitter.query.add_directive("set-adjacent!", select_adjacent, {})

--- First row of the comment nodes directly above `node`, each ending on or
--- after the row before the next; `node`'s start row when there are none.
--- @param node TSNode
--- @return integer 0-based
local function comments_start(node)
    local row = node:start()
    local sibling = node:prev_sibling()
    while sibling and sibling:type():find("comment") and sibling:end_() + 1 >= row do
        row = sibling:start()
        sibling = sibling:prev_sibling()
    end
    return row
end

--- Definitions named `symbol` in `src` per the tags query of `lang`: matches
--- with a `@definition.<kind>` capture and a `@name` capture whose text is
--- `symbol`.
--- @param src string
--- @param lang string
--- @param symbol string
--- @param path string file `src` is the content of
--- @return symbol_ls.Hit[] with UTF-8 ranges, `definition` and no `client`; `{}` without a tags query
function M.definitions(src, lang, symbol, path)
    local query = vim.treesitter.query.get(lang, "tags")
    if not query then return {} end
    local root = vim.treesitter.get_string_parser(src, lang):parse()[1]:root()
    local hits = {} --- @type symbol_ls.Hit[]
    for _, match, metadata in query:iter_matches(root, src) do
        local name, definition, kind, doc
        for id, nodes in pairs(match) do
            local capture = query.captures[id]
            local suffix = capture:match("^definition%.(.+)$")
            if suffix then
                kind, definition = suffix, nodes[1]
            elseif capture == "name" then
                name = nodes[1]
            elseif capture == "doc" and metadata[id] and metadata[id].text ~= "" then
                doc = metadata[id].text
            end
        end
        if definition and name and vim.treesitter.get_node_text(name, src) == symbol then
            local row, col, end_row, end_col = name:range()
            hits[#hits + 1] = {
                info = {
                    name = symbol,
                    kind = KINDS[kind] or SymbolKind.Variable,
                    location = {
                        uri = vim.uri_from_fname(path),
                        range = { start = { line = row, character = col }, ["end"] = { line = end_row, character = end_col } },
                    },
                },
                definition = {
                    range = { definition:range() },
                    first_row = doc and definition:start() or comments_start(definition),
                    doc = doc,
                },
            }
        end
    end
    return hits
end

return M
