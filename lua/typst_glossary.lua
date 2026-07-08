---Resolve typst glossary `@term` refs to their definition in the source dict.
---
---A glossy `@ts` ref points at a dict pair `ts: ( short: ..., ... )` living in a
---`.typ` file imported (transitively) by the document. tinymist resolves such
---refs to the glossy package's runtime `label()` placeholder, not the entry, so
---`grd` needs this custom path for glossary keys and defers to the LSP otherwise.
---All parsing is treesitter over the typst grammar — no line patterns, no `rg`.
local M = {}

---Strip a leading and trailing `"` from a treesitter string-node text.
---@param s string
---@return string
local function dequote(s)
    return (s:gsub('^"(.*)"$', "%1"))
end

---Glossary key referenced by a ref node's text.
---`@ts`/`@ts:pl`/`@ts:short` all key `ts`; `@fig:f` keys `fig`. Robust across
---hyphenated/unicode keys (`@wagner-meerwein`, `@π-bond`).
---@param text string ref node text (leading `@`)
---@return string key
function M.ref_key(text)
    return text:gsub("^@", ""):match("^[^:]+")
end

---Relative-import file paths declared in typst source, resolved against `dir`.
---Package imports (`@preview/…`, `@local/…`) are dropped; file imports (bare or
---`./`-prefixed) are kept and normalised to absolute paths. Pure path math — no
---filesystem access, so callers get deterministic results for symlink handling.
---@param content string typst source
---@param dir string directory the source lives in (absolute)
---@return string[] paths
function M.imports(content, dir)
    local root = vim.treesitter.get_string_parser(content, "typst"):parse()[1]:root()
    local query = vim.treesitter.query.parse("typst", "(import (string) @path)")
    local paths = {}
    for _, node in query:iter_captures(root, content) do
        local path = dequote(vim.treesitter.get_node_text(node, content))
        if not path:match("^@") then
            table.insert(paths, vim.fs.normalize(vim.fs.joinpath(dir, path)))
        end
    end
    return paths
end

---Glossary dict entries in typst source, mapping key → its 0-indexed `{row, col}`.
---Matches entry-level pairs `key: (…)` (bare or quoted key) via the `(group)`
---value constraint, which excludes string/content fields (`short: "TS"`).
---@param content string typst source
---@return table<string, [integer, integer]> positions
function M.entries(content)
    local root = vim.treesitter.get_string_parser(content, "typst"):parse()[1]:root()
    local query = vim.treesitter.query.parse("typst", "(tagged [(ident) (string)] @key (group))")
    local positions = {}
    for _, node in query:iter_captures(root, content) do
        local row, col = node:range()
        positions[dequote(vim.treesitter.get_node_text(node, content))] = { row, col }
    end
    return positions
end

---Glossary field values of the entry keyed `key`, or nil if no such entry.
---Scopes to the matched entry's `(group)` and pulls its `field: value` pairs
---whose value is a string or content block (the fields a hover shows), keyed by
---field name. Content keeps its raw typst source (brackets stripped) so callers
---can render it as typst; strings are dequoted.
---@param content string typst source
---@param key string entry key
---@return table<string, { text: string, kind: "string"|"content" }>|nil fields
function M.fields(content, key)
    local root = vim.treesitter.get_string_parser(content, "typst"):parse()[1]:root()
    local entry_q = vim.treesitter.query.parse("typst", "(tagged [(ident) (string)] @key (group) @group)")
    local group
    for _, match in entry_q:iter_matches(root, content) do
        local key_node = match[1][1] -- @key is the first capture
        if dequote(vim.treesitter.get_node_text(key_node, content)) == key then
            group = match[2][1] -- @group
            break
        end
    end
    if not group then return nil end

    local field_q = vim.treesitter.query.parse("typst", "(tagged [(ident) (string)] @field [(content) (string)] @value)")
    local fields = {}
    for _, match in field_q:iter_matches(group, content) do
        local name = dequote(vim.treesitter.get_node_text(match[1][1], content))
        local value_node = match[2][1]
        local text = vim.treesitter.get_node_text(value_node, content)
        local kind = value_node:type()
        text = kind == "content" and text:gsub("^%[(.*)%]$", "%1") or dequote(text)
        fields[name] = { text = vim.trim(text), kind = kind }
    end
    return fields
end

---Hover markdown for glossary entry `fields`, or nil if there's nothing to show.
---Shows the short form (`short-fmt` over `short`, which can differ from the ref
---key) — em-dash — long form (`long-fmt` over `long`, omitted when neither
---exists), then the description. Content-typed values render in a ```typst```
---fence for typst highlighting via the standard markdown-hover treesitter path;
---plain-string descriptions render as a markdown paragraph.
---@param fields table<string, { text: string, kind: "string"|"content" }>
---@return string|nil markdown
function M.render(fields)
    local short = fields["short-fmt"] or fields["short"]
    local long = fields["long-fmt"] or fields["long"]
    local desc = fields["description"]
    if not (short or long) then return nil end

    local signature = short and short.text or ""
    if long then signature = signature .. " — " .. long.text end
    local lines = { "```typst", signature, "```" }
    if desc then
        table.insert(lines, "")
        if desc.kind == "content" then
            vim.list_extend(lines, { "```typst", desc.text, "```" })
        else
            table.insert(lines, desc.text)
        end
    end
    return table.concat(lines, "\n")
end

---Walk the relative-import graph from `bufnr`, collecting `extract`'s non-nil
---results. The start buffer is read from memory (unsaved edits), imports from
---disk; the buffer is included so refs between entries inside a glossary file
---resolve (descriptions cite sibling terms). BFS follows only relative file
---imports, so sibling docs and package deps stay unreachable.
---@param bufnr integer
---@param extract fun(path: string, content: string, lines: string[]): any|nil
---@return any[] results in BFS order
function M.walk(bufnr, extract)
    local buf_path = vim.fn.resolve(vim.api.nvim_buf_get_name(bufnr))
    local start_lines = vim.api.nvim_buf_get_lines(bufnr, 0, -1, false)
    local results = {}
    local visited = {}
    local queue = { buf_path }
    while #queue > 0 do
        local path = vim.fn.resolve(table.remove(queue, 1))
        if not visited[path] then
            visited[path] = true
            local lines = path == buf_path and start_lines
                or (vim.uv.fs_stat(path) or {}).type == "file" and vim.fn.readfile(path)
            if lines then
                local content = table.concat(lines, "\n")
                local result = extract(path, content, lines)
                if result ~= nil then table.insert(results, result) end
                vim.list_extend(queue, M.imports(content, vim.fs.dirname(path)))
            end
        end
    end
    return results
end

---Key of the glossary `ref` node covering position `(row, col)`, or nil.
---@param bufnr integer
---@param row integer 0-indexed
---@param col integer 0-indexed
---@return string|nil key
local function ref_key_at(bufnr, row, col)
    local ok, parser = pcall(vim.treesitter.get_parser, bufnr, "typst")
    if not ok or not parser then return nil end
    parser:parse({ row, row }) -- get_node doesn't parse on its own
    local node = vim.treesitter.get_node({ bufnr = bufnr, pos = { row, col } })
    while node and node:type() ~= "ref" do node = node:parent() end
    if not node then return nil end
    return M.ref_key(vim.treesitter.get_node_text(node, bufnr))
end

---Definition qf items for the glossary ref under the cursor, or nil to defer.
---Returns nil when the cursor isn't on a ref (LSP owns `@fig:`/`@sec:`/…) or the
---key isn't a glossary entry in the buffer or its import graph.
---@param bufnr integer
---@return table[]|nil items `:h setqflist-what` items
function M.resolve(bufnr)
    local row, col = unpack(vim.api.nvim_win_get_cursor(0))
    local key = ref_key_at(bufnr, row - 1, col)
    if not key then return nil end
    local items = M.walk(bufnr, function(path, content, lines)
        local pos = M.entries(content)[key]
        if pos then
            return {
                filename = path,
                lnum = pos[1] + 1,
                col = pos[2] + 1,
                text = vim.trim(lines[pos[1] + 1] or ""),
            }
        end
    end)
    if #items > 0 then return items end
    return nil
end

---Hover markdown for the glossary ref at position `(row, col)`, or nil to defer.
---nil when the position isn't on a ref, or its key isn't a glossary entry in the
---buffer or its import graph — leaving the hover to other LSP clients.
---@param bufnr integer
---@param row integer 0-indexed
---@param col integer 0-indexed
---@return string|nil markdown
function M.hover(bufnr, row, col)
    local key = ref_key_at(bufnr, row, col)
    if not key then return nil end
    local fields = M.walk(bufnr, function(_, content) return M.fields(content, key) end)[1]
    return fields and M.render(fields)
end

return M
