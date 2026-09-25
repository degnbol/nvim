
local M = {}

---Is this terminal kitty?
function M.term()
    return vim.env.KITTY_PID ~= nil
end

---Cell layout. All numbers except `aspect` are fractions of the cell height.
---@class kitty.CellFractions
---@field postscript_name string medium face
---@field aspect number cell height / cell width
---@field below number baseline to cell bottom
---@field x number ink height of "x"
---@field cap number ink height of "H"
---@field desc number ink height of "p" minus that of "x"

local fractions ---@type kitty.CellFractions|nil

---Kitty's cell layout for the configured font. Independent of zoom.
---Runs `kitty +launch` synchronously, and errors if `kitty` is not on PATH.
---A success is cached for the session, a failure is not.
---@return kitty.CellFractions|nil fractions
---@return string|nil err kitty's stderr on failure
function M.cell_fractions()
    if fractions then return fractions end
    local script = vim.fn.stdpath("config") .. "/scripts/kitty_cell_fractions.py"
    local res = vim.system({ "kitty", "+launch", script }, { text = true }):wait()
    if res.code ~= 0 then return nil, res.stderr end
    fractions = vim.json.decode(res.stdout) -- Cache
    return fractions
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
