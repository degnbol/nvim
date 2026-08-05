-- Augment LSP rename with structural string references the server treats as
-- opaque text (e.g. argparse/Tap `add_argument` option strings). A per-client
-- wrapper around the `textDocument/rename` handler adds `TextEdit`s to the
-- server's `WorkspaceEdit` *before* it is applied, so string and identifier
-- edits land in one undo block. Only files the server already touches are
-- augmented — see notes/PLAN-rename-string-refs.md.

local lsp = require "utils.lsp"

local M = {}

--- Filetype → matchers. Each matcher exposes `edits(src, old, new) ->
--- lsp_rename.Edit[]`, pure over the source string.
local matchers = {
    python = { require "lsp_rename.argparse" },
}

--- Convert a matcher's byte-column edit to an `lsp.TextEdit` in `encoding`.
--- @param bufnr integer
--- @param edit lsp_rename.Edit
--- @param encoding string
--- @return lsp.TextEdit
local function to_text_edit(bufnr, edit, encoding)
    local util = vim.lsp.util
    --- @type lsp.TextEdit
    local text_edit = {
        range = {
            start = { line = edit.srow, character = util.character_offset(bufnr, edit.srow, edit.scol, encoding) },
            ["end"] = { line = edit.erow, character = util.character_offset(bufnr, edit.erow, edit.ecol, encoding) },
        },
        newText = edit.newText,
    }
    return text_edit
end

--- Add matcher string edits to a rename `WorkspaceEdit`, in place. The old name
--- is the buffer text under a server edit; the new name is that edit's newText.
--- @param workspace_edit lsp.WorkspaceEdit
--- @param ctx lsp.HandlerContext
local function augment(workspace_edit, ctx)
    local client = assert(vim.lsp.get_client_by_id(ctx.client_id))
    local encoding = client.offset_encoding
    local files = lsp.workspace_edit_files(workspace_edit)

    local old, new
    for _, file in ipairs(files) do
        local edit = file.edits[1]
        if edit then
            local bufnr = vim.uri_to_bufnr(file.uri)
            vim.fn.bufload(bufnr)
            old = lsp.range_text(bufnr, edit.range, encoding)
            new = edit.newText
            break
        end
    end
    if not old or old == new then return end

    for _, file in ipairs(files) do
        local bufnr = vim.uri_to_bufnr(file.uri)
        vim.fn.bufload(bufnr)
        local fts = matchers[vim.bo[bufnr].filetype]
        if fts then
            local src = table.concat(vim.api.nvim_buf_get_lines(bufnr, 0, -1, false), "\n")
            for _, matcher in ipairs(fts) do
                for _, edit in ipairs(matcher.edits(src, old, new)) do
                    file.edits[#file.edits + 1] = to_text_edit(bufnr, edit, encoding)
                end
            end
        end
    end
end

--- Clients whose rename handler is already wrapped. Weak keys so entries drop
--- when a stopped client is collected.
--- @type table<vim.lsp.Client, true>
local installed = setmetatable({}, { __mode = "k" })

--- Whether any matcher handles `filetype` (gate for `install`).
--- @param filetype string
--- @return boolean
function M.supports(filetype)
    return matchers[filetype] ~= nil
end

--- Wrap a client's `textDocument/rename` handler to augment its `WorkspaceEdit`
--- with structural string edits. Idempotent per client: `client.handlers` is a
--- single shared object, so re-wrapping would run the augment twice.
--- @param client vim.lsp.Client
function M.install(client)
    if installed[client] then return end
    installed[client] = true

    local method = "textDocument/rename"
    local default = client.handlers[method] or vim.lsp.handlers[method]
    client.handlers[method] = function(err, result, ctx, config)
        if not err and result then
            local ok, msg = pcall(augment, result, ctx)
            if not ok then
                vim.notify("lsp_rename: " .. tostring(msg), vim.log.levels.ERROR)
            end
        end
        return default(err, result, ctx, config)
    end
end

return M
