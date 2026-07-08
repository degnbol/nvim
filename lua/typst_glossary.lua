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

---Definition qf items for the glossary ref under the cursor, or nil to defer.
---Returns nil when the cursor isn't on a ref (LSP owns `@fig:`/`@sec:`/…) or the
---key isn't a glossary entry in the buffer or its import graph. The buffer is
---searched too so refs between entries inside a glossary file itself resolve
---(descriptions reference sibling terms). BFS walks only the relative imports
---the document declares, so sibling docs and package deps stay unreachable.
---@param bufnr integer
---@return table[]|nil items `:h setqflist-what` items
function M.resolve(bufnr)
    local node = vim.treesitter.get_node()
    while node and node:type() ~= "ref" do node = node:parent() end
    if not node then return nil end
    local key = M.ref_key(vim.treesitter.get_node_text(node, bufnr))

    -- Unsaved buffer edits: read the start buffer from memory, imports from disk.
    local buf_path = vim.fn.resolve(vim.api.nvim_buf_get_name(bufnr))
    local start_lines = vim.api.nvim_buf_get_lines(bufnr, 0, -1, false)

    local items = {}
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
                local pos = M.entries(content)[key]
                if pos then
                    table.insert(items, {
                        filename = path,
                        lnum = pos[1] + 1,
                        col = pos[2] + 1,
                        text = vim.trim(lines[pos[1] + 1] or ""),
                    })
                end
                vim.list_extend(queue, M.imports(content, vim.fs.dirname(path)))
            end
        end
    end
    if #items > 0 then return items end
    return nil
end

return M
