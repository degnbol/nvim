local util = require "utils/init"

local M = {}

---Add global user command.
---Convenience link to nvim_create_user_command and doesn't require opts to be at least {}.
---@param name string
---@param command string|fun(args: vim.api.keyset.create_user_command.command_args)
---@param opts? vim.api.keyset.user_command
M.cmd = function (name, command, opts)
    opts = opts or {}
    vim.api.nvim_create_user_command(name, command, opts)
end

M.map = vim.keymap.set
-- This allows for importing M as `map` and directly calling `map("n", ...)` for convenience.
setmetatable(M, {
    __call = function(_, ...)
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

---Visual map including select mode.
---@param lhs string
---@param rhs string|function
---@param desc? string
---@param opts? table
function M.v(lhs, rhs, desc, opts)
    opts = opts or {}
    opts.desc = desc
    vim.keymap.set('v', lhs, rhs, opts)
end

---Visual map excluding select mode.
---@param lhs string
---@param rhs string|function
---@param desc? string
---@param opts? table
function M.x(lhs, rhs, desc, opts)
    opts = opts or {}
    opts.desc = desc
    vim.keymap.set('x', lhs, rhs, opts)
end

---Operator pending map.
---@param lhs string
---@param rhs string|function
---@param desc? string
---@param opts? table
function M.o(lhs, rhs, desc, opts)
    opts = opts or {}
    opts.desc = desc
    vim.keymap.set('o', lhs, rhs, opts)
end

---Normal and visual map excluding select mode.
---@param lhs string
---@param rhs string|function
---@param desc? string
---@param opts? table
function M.nx(lhs, rhs, desc, opts)
    opts = opts or {}
    opts.desc = desc
    vim.keymap.set({ 'n', 'x' }, lhs, rhs, opts)
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

---Operator pending mode and for visual.
---@param lhs string
---@param rhs string|function
---@param desc? string
---@param opts? table
function M.ox(lhs, rhs, desc, opts)
    opts = opts or {}
    opts.desc = desc
    vim.keymap.set({ 'o', 'x' }, lhs, rhs, opts)
end

---Normal, operator pending mode and for visual.
---@param lhs string
---@param rhs string|function
---@param desc? string
---@param opts? table
function M.nox(lhs, rhs, desc, opts)
    opts = opts or {}
    opts.desc = desc
    vim.keymap.set({ 'n', 'o', 'x' }, lhs, rhs, opts)
end

---Buffer-local map. Mode can be string or table.
---@param mode string|table
---@param lhs string
---@param rhs string|function
---@param desc? string
---@param opts? table
function M.buf(mode, lhs, rhs, desc, opts)
    opts = opts or {}
    opts.desc = desc
    opts.buffer = true
    vim.keymap.set(mode, lhs, rhs, opts)
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
