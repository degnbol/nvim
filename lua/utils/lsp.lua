local qf = require "utils/qf"

local M = {}

--- Normalize a `WorkspaceEdit` to a shape-agnostic list of its per-file edit
--- arrays. A `WorkspaceEdit` carries edits either as `changes` (a `uri → edits`
--- map) or as `documentChanges` (a list of `{ textDocument, edits }`, possibly
--- interleaved with resource operations that have no `edits`). This returns the
--- same `{ uri, edits }` view for both shapes, where each `edits` is the *live*
--- array from the `WorkspaceEdit` — appending to it mutates the original.
--- @param workspace_edit lsp.WorkspaceEdit
--- @return { uri: string, edits: lsp.TextEdit[] }[] files
function M.workspace_edit_files(workspace_edit)
    --- @type { uri: string, edits: lsp.TextEdit[] }[]
    local files = {}
    if workspace_edit.changes then
        for uri, edits in pairs(workspace_edit.changes) do
            files[#files + 1] = { uri = uri, edits = edits }
        end
    elseif workspace_edit.documentChanges then
        for _, change in ipairs(workspace_edit.documentChanges) do
            if change.edits and change.textDocument then
                files[#files + 1] = { uri = change.textDocument.uri, edits = change.edits }
            end
        end
    end
    return files
end

--- Buffer text spanned by a single-line LSP `range`, resolving its
--- offset-encoding character columns to byte columns against the (loaded)
--- buffer. `range` must not span lines (rename targets never do).
--- @param bufnr integer
--- @param range lsp.Range
--- @param encoding string  client offset_encoding
--- @return string
function M.range_text(bufnr, range, encoding)
    local line = vim.api.nvim_buf_get_lines(bufnr, range.start.line, range.start.line + 1, false)[1]
    local scol = vim.str_byteindex(line, encoding, range.start.character)
    local ecol = vim.str_byteindex(line, encoding, range["end"].character)
    return line:sub(scol + 1, ecol)
end

--- Whether LSP position `a` comes before `b`, or is equal to it.
--- @param a lsp.Position
--- @param b lsp.Position
--- @return boolean
local function position_le(a, b)
    return a.line < b.line or (a.line == b.line and a.character <= b.character)
end

--- Target of the LSP document link covering a buffer position, from any
--- attached client that serves `textDocument/documentLink`. Links without a
--- resolved `target` are skipped. Blocks for up to 1 s while clients respond;
--- a timeout or a server error notifies at WARN.
--- @param buf integer
--- @param row integer 0-indexed
--- @param col integer 0-indexed byte column
--- @return string|nil target URI of the first covering link, nil if none
function M.document_link_at(buf, row, col)
    local method = "textDocument/documentLink"
    -- buf_request_sync without a client waits out its whole timeout.
    if #vim.lsp.get_clients { bufnr = buf, method = method } == 0 then return nil end
    local line = vim.api.nvim_buf_get_lines(buf, row, row + 1, true)[1]
    local params = { textDocument = vim.lsp.util.make_text_document_params(buf) }
    local responses, err = vim.lsp.buf_request_sync(buf, method, params)
    if not responses then
        vim.notify(method .. ": " .. tostring(err), vim.log.levels.WARN)
        return nil
    end
    for client_id, response in pairs(responses) do
        local client = assert(vim.lsp.get_client_by_id(client_id))
        if response.err then
            vim.notify(("%s (%s): %s"):format(method, client.name, response.err.message), vim.log.levels.WARN)
        end
        --- @type lsp.Position
        local pos = { line = row, character = vim.str_utfindex(line, client.offset_encoding, col, false) }
        for _, link in ipairs(response.result or {}) do
            if link.target and position_le(link.range.start, pos) and not position_le(link.range["end"], pos) then
                return link.target
            end
        end
    end
    return nil
end

--- LSP list options whose `on_list` drops the items that match any predicate
--- and passes the rest to `qf.jump_or_load`.
--- @param ... fun(item: table): boolean Predicates on a `:h setqflist-what` item.
--- True drops the item.
--- @return vim.lsp.LocationOpts opts
local function filtered_list_opts(...)
    local rejects = { ... }
    return {
        on_list = function(options)
            options.items = vim.tbl_filter(function(item)
                return not vim.iter(rejects):any(function(reject) return reject(item) end)
            end, options.items)
            qf.jump_or_load(options)
        end
    }
end

--- Go to the LSP definition of the symbol under the cursor, with the results
--- handled by `qf.jump_or_load`.
function M.definition()
    vim.lsp.buf.definition { on_list = qf.jump_or_load }
end

--- Go to the LSP references of the symbol under the cursor, without the items
--- that match any predicate, with the rest handled by `qf.jump_or_load`.
--- @param ... fun(item: table): boolean Predicates on a `:h setqflist-what` item.
--- True drops the item.
function M.references(...)
    vim.lsp.buf.references(nil, filtered_list_opts(...))
end

---Create a root_dir function for vim.lsp.config that resolves symlinks before searching.
---Needed for symlinked dotfiles where .git may not be in the apparent ancestor chain.
---@param markers string|string[] Root markers to search for (default: { '.git' })
---@return function root_dir_fn Function compatible with vim.lsp.config root_dir
function M.symlink_root_dir(markers)
    markers = markers or { '.git' }
    if type(markers) == 'string' then markers = { markers } end
    return function(bufnr, on_dir)
        local fname = vim.api.nvim_buf_get_name(bufnr)
        local resolved = vim.fn.resolve(fname)
        local root = vim.fs.root(resolved, markers)
        if root then
            on_dir(root)
        else
            -- Fallback to file's directory for single-file support
            on_dir(vim.fs.dirname(resolved))
        end
    end
end

return M
