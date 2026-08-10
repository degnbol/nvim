-- The `:Docstring` consumer: reconcile the docstring of the function under the
-- cursor against its signature, in the buffer. Ties the pieces together —
-- treesitter locators (lang), parse/render (style), and the reconcile algorithm
-- — and writes the result back, re-indented to the function body.

local python = require "docstring.lang.python"
local google = require "docstring.style.google"
local reconcile = require "docstring.reconcile"
local model_mod = require "docstring.model"

local M = {}

--- The leading whitespace of a buffer line.
--- @param bufnr integer
--- @param row integer  0-indexed
--- @return string
local function line_indent(bufnr, row)
    local line = vim.api.nvim_buf_get_lines(bufnr, row, row + 1, true)[1]
    return line:match("^%s*")
end

--- The buffer lines of a Python `"""…"""` literal for dedented docstring `text`.
--- The first line carries no indent (it lands after the existing indent, whether
--- placed by `set_text` at the opening-quote column or by the caller); each
--- continuation and the closing line carry `indent`; blank lines stay empty to
--- avoid trailing whitespace.
--- @param text string  dedented Google-style docstring text
--- @param indent string
--- @return string[]
local function literal_lines(text, indent)
    local lines = vim.split(text, "\n", { plain = true })
    if #lines == 1 then
        return { '"""' .. lines[1] .. '"""' }
    end
    local out = { '"""' .. lines[1] }
    for i = 2, #lines do
        out[#out + 1] = lines[i] == "" and "" or indent .. lines[i]
    end
    out[#out + 1] = indent .. '"""'
    return out
end

--- Reconcile the docstring of the function under the cursor with its signature,
--- editing the current buffer in place. Replaces an existing docstring or, if
--- there is none, inserts a fresh one before the first body statement.
function M.reconcile_at_cursor()
    local bufnr = vim.api.nvim_get_current_buf()
    local target = python.locate_target()
    if not target then
        vim.notify("Docstring: cursor is not in a function", vim.log.levels.WARN)
        return
    end

    local params = python.signature_params(target, bufnr)
    local node, text = python.locate_docstring(target, bufnr)
    local model = text and google.parse(text) or model_mod.empty()
    model = reconcile.reconcile(params, model)
    local rendered = google.render(model)

    if node then
        local srow, scol, erow, ecol = node:range()
        vim.api.nvim_buf_set_text(bufnr, srow, scol, erow, ecol, literal_lines(rendered, line_indent(bufnr, srow)))
    else
        local row = python.body_row(target, bufnr)
        if not row then
            vim.notify("Docstring: cannot insert into a one-line function body", vim.log.levels.WARN)
            return
        end
        local indent = line_indent(bufnr, row)
        local lines = literal_lines(rendered, indent)
        lines[1] = indent .. lines[1]
        vim.api.nvim_buf_set_lines(bufnr, row, row, false, lines)
    end
end

return M
