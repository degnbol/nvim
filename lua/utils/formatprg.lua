local M = {}

---'formatexpr' that filters the lines through 'formatprg', as the built-in `gq`
---does, but notifies the command's stderr instead of inserting it into the buffer
---(a built-in filter merges stderr into the output). The command runs with
---`pipefail`, and on a non-zero exit the lines are left unchanged — otherwise a
---failing first stage would replace them with the empty output of the next.
---Falls back to internal formatting when 'formatprg' is empty, and in Insert and
---Replace mode, where 'formatexpr' is called for auto-wrap at 'textwidth'.
---Use as `formatexpr=v:lua.require'utils.formatprg'.formatexpr()` (`v:lua` does
---not parse a `/` in the module name).
---@return integer 0 when handled, 1 to fall back to internal formatting
function M.formatexpr()
    local prg = vim.o.formatprg
    if prg == "" or vim.fn.mode():find("^[iR]") then return 1 end
    local first = vim.v.lnum
    local last = first + vim.v.count - 1
    local lines = vim.api.nvim_buf_get_lines(0, first - 1, last, false)
    local cmd = { vim.o.shell, vim.o.shellcmdflag, "set -o pipefail; " .. prg }
    local obj = vim.system(cmd, { stdin = table.concat(lines, "\n") .. "\n", text = true }):wait()
    if obj.code ~= 0 then
        require("utils/init").notify_failure({ prg }, obj)
        return 0
    end
    local stderr = vim.trim(obj.stderr)
    if stderr ~= "" then vim.notify(stderr, vim.log.levels.WARN) end
    local formatted = vim.split(obj.stdout, "\n")
    if formatted[#formatted] == "" then table.remove(formatted) end
    if not vim.deep_equal(formatted, lines) then
        vim.api.nvim_buf_set_lines(0, first - 1, last, false, formatted)
    end
    return 0
end

return M
