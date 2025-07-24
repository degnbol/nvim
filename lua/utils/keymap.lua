local util = require "utils/init"

local M = {}

M.map = vim.keymap.set
setmetatable(M, {
    __call = function(self, ...)
        vim.keymap.set(...)
    end
})

---Normal map
---@param lhs string
---@param rhs string|function
---@param desc? string
---@param opts? table
function M.n(lhs, rhs, desc, opts)
    opts = opts or {}
    opts.desc = desc
    vim.keymap.set('n', lhs, rhs, opts)
end

---Insert map
---@param lhs string
---@param rhs string|function
---@param desc? string
---@param opts? table
function M.i(lhs, rhs, desc, opts)
    opts = opts or {}
    opts.desc = desc
    vim.keymap.set('i', lhs, rhs, opts)
end

---Commandline map
---@param lhs string
---@param rhs string|function
---@param desc? string
---@param opts? table
function M.c(lhs, rhs, desc, opts)
    opts = opts or {}
    opts.desc = desc
    vim.keymap.set('c', lhs, rhs, opts)
end

---Visual map excluding select mode.
---@param lhs string
---@param rhs string|function
---@param desc? string
---@param opts? table
function M.v(lhs, rhs, desc, opts)
    opts = opts or {}
    opts.desc = desc
    vim.keymap.set('v', lhs, rhs, opts)
end

---Visual map including select mode.
---@param lhs string
---@param rhs string|function
---@param desc? string
---@param opts? table
function M.x(lhs, rhs, desc, opts)
    opts = opts or {}
    opts.desc = desc
    vim.keymap.set('x', lhs, rhs, opts)
end

---Blockwise visual map, i.e. <C-v> mappings.
---@param lhs string
---@param rhs string|function
---@param desc? string
---@param opts? table
function M.cv(lhs, rhs, desc, opts)
    opts = opts or {}
    opts.desc = desc
    if type(rhs) == "string" then
        opts.expr = true
        vim.keymap.set('v', lhs, function()
            if not util.is_visual_blockwise() then return "" end
            return rhs
        end, opts)
    else
        vim.keymap.set('v', lhs, function()
            if not util.is_visual_blockwise() then return end
            return rhs()
        end, opts)
    end
end

---Add desc(ription) to an already defined keymap.
---@param mode string
---@param lhs string
---@param desc string
M.desc = function(mode, lhs, desc)
    pcall(require "mini.clue".set_keymap_desc, mode, lhs, desc)
end

return M
