#!/usr/bin/env lua
local M = {}

---Is this terminal kitty?
function M.term()
    return os.execute("kitty @ ls 2> /dev/null > /dev/null")
end

---Enable/diable ligatures for the current kitty window.
---@param enable boolean
---@return boolean success
function M.ligatures(enable)
    if enable then
        return os.execute("kitty @ disable-ligatures never") or false
    else
        return os.execute("kitty @ disable-ligatures always") or false
    end
end

---Enable/disable ligatures for the current kitty window when a buffer with given filename pattern is in focus.
---E.g. set it for a filetype.
---@param enable boolean
---@param pattern string
---@return boolean success
function M.ligatures_pattern(enable, pattern)
    if not M.term() then return false end
    local grp = vim.api.nvim_create_augroup("kitty_ligatures", {clear=true})
    vim.api.nvim_create_autocmd("BufEnter", {
        pattern = pattern,
        group = grp,
        callback = function () M.ligatures(enable) end
    })
    vim.api.nvim_create_autocmd({"BufLeave", "BufWinLeave", "BufDelete"}, {
        pattern = pattern,
        group = grp,
        callback = function () M.ligatures(not enable) end
    })
    return true
end

return M
