local util = require "utils/init"

local M = {}

-- This allows for importing M as `map` and directly calling `map("n", ...)` for convenience.
setmetatable(M, {
    __call = function(_, ...)
        vim.keymap.set(...)
    end
})

---Make a map function for a fixed mode.
---@param mode string|string[]
---@return fun(lhs: string, rhs: string|function, desc: string|nil, opts: table|nil)
local function mode_map(mode)
    return function(lhs, rhs, desc, opts)
        vim.keymap.set(mode, lhs, rhs, vim.tbl_extend("force", opts or {}, { desc = desc }))
    end
end

---Normal map
M.n = mode_map 'n'
---Insert map
M.i = mode_map 'i'
---Commandline map
M.c = mode_map 'c'
---Visual map including select mode.
M.v = mode_map 'v'
---Visual map excluding select mode.
M.x = mode_map 'x'
---Operator pending map.
M.o = mode_map 'o'
---Normal and visual map excluding select mode.
M.nx = mode_map { 'n', 'x' }
---Operator pending mode and for visual.
M.ox = mode_map { 'o', 'x' }
---Normal, operator pending mode and for visual.
M.nox = mode_map { 'n', 'o', 'x' }

---Blockwise visual map, i.e. <C-v> mappings. The key does nothing in other
---visual modes and in select mode.
---@param lhs string
---@param rhs string|function
---@param desc string|nil
---@param opts table|nil
function M.cv(lhs, rhs, desc, opts)
    if type(rhs) == "string" then
        M.v(lhs, function()
            if not util.is_visual_blockwise() then return "" end
            return rhs
        end, desc, vim.tbl_extend("force", opts or {}, { expr = true }))
    else
        M.v(lhs, function()
            if not util.is_visual_blockwise() then return end
            return rhs()
        end, desc, opts)
    end
end

---Buffer-local map.
---@param mode string|string[]
---@param lhs string
---@param rhs string|function
---@param desc string|nil
---@param opts table|nil
function M.buf(mode, lhs, rhs, desc, opts)
    mode_map(mode)(lhs, rhs, desc, vim.tbl_extend("force", opts or {}, { buf = 0 }))
end

---Add desc(ription) to an already-defined keymap, so a plugin's `after`
---callback can annotate its own keymaps without having to set them up
---inside mini.clue's main config. Returns false when the mapping doesn't
---exist for at least one of the given modes.
---@param mode string|table
---@param lhs string
---@param desc string
---@return boolean ok true iff `desc` was applied for every mode
function M.desc(mode, lhs, desc)
    --- @type string[]
    local modes = type(mode) == "table" and mode or { mode }
    local all_ok = true
    for _, m in ipairs(modes) do
        local map_data = vim.fn.maparg(lhs, m, false, true)
        if vim.tbl_count(map_data) > 0 then
            map_data.desc = desc
            vim.fn.mapset(m, false, map_data)
        else
            all_ok = false
        end
    end
    return all_ok
end

return M
