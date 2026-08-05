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

return M
