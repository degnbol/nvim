--- Surface conform formatter errors as notifications and inline diagnostics.
--- For execution errors conform mutes the formatter's stderr, notifying only
--- "Formatter failed. See :ConformInfo" (conform.nvim init.lua:471). This
--- restores the real message and, when the stderr carries a compiler-style
--- line:col, pins it as a diagnostic on that line — the only structured
--- feedback available for formatters without a matching linter/LSP.
local M = {}

local ns = vim.api.nvim_create_namespace("conform_error")

--- Extract a compiler-style position from an error message.
--- Matches the first `:line:col:` (or `:line:`) run, so it is agnostic to the
--- leading filename (absolute path, `<standard input>`, etc.).
--- @param message string
--- @return number|nil line 1-based; nil when the message carries no position
--- @return number|nil col 1-based; 1 when the message carries no column
--- @return string|nil message Text following the position, trailing space trimmed
function M.parse(message)
    local line, col, rest = message:match(":(%d+):(%d+): ?(.*)")
    if not line then
        line, rest = message:match(":(%d+): ?(.*)")
        col = "1"
    end
    if not line then return nil end
    local text = (rest or ""):gsub("%s+$", "")
    return tonumber(line), tonumber(col), text
end

--- Report a conform.format callback error, if any.
--- Notifies the real stderr and pins an in-range position as a diagnostic;
--- resets prior diagnostics so a subsequent clean format clears them. Only acts
--- on execution errors — conform already notifies timeouts and the like, and
--- both execution-error messages start with "Formatter '<name>' error"
--- (conform.nvim runner.lua:462,471).
--- @param message string|nil The err passed to conform.format's callback: a
---   message string on failure, nil on success
--- @param bufnr number
function M.report(message, bufnr)
    if not vim.api.nvim_buf_is_valid(bufnr) then return end
    vim.diagnostic.reset(ns, bufnr)
    if not message or not message:match("^Formatter '.-' error") then return end
    -- Notify "<formatter>: <text>" — the position lives on the line as the
    -- diagnostic below, so it stays out of the notification. Fall back to the
    -- trimmed raw message when there is no position to parse.
    local line, col, text = M.parse(message)
    local name = message:match("^Formatter '(.-)'")
    if line and text ~= "" then
        vim.notify(("%s: %s"):format(name, text), vim.log.levels.ERROR)
    else
        vim.notify((message:gsub("%s+$", "")), vim.log.levels.ERROR)
    end
    if not line or line > vim.api.nvim_buf_line_count(bufnr) then return end
    vim.diagnostic.set(ns, bufnr, {
        {
            lnum = line - 1,
            col = col - 1,
            message = text ~= "" and text or message,
            severity = vim.diagnostic.severity.ERROR,
            source = name,
        },
    })
end

return M
