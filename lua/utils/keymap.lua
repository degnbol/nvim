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

---Add desc(ription) to an already defined keymap.
---@param mode string
---@param lhs string
---@param desc string
function M.desc(mode, lhs, desc)
    pcall(require "mini.clue".set_keymap_desc, mode, lhs, desc)
end

---Return whether an item from with the dict format described by
---:h setqflist-what
---is referring to the line we are currently on with the cursor.
---Useful for filtering lsp results, e.g. goto references and goto definition.
---@param item table
---@return boolean
function M.qf_item_is_self(item)
    return item.filename == vim.api.nvim_buf_get_name(0) and item.lnum == vim.api.nvim_win_get_cursor(0)[1]
end

---Open quickfix with height set to number of entries if less than default 
---height (10) and if only one entry, jump to it without opening qf.
---@param options vim.fn.setqflist.what `:h setqflist-what`
function M.qf_mini(options)
    if #options.items == 1 then
        local item = options.items[1]
        vim.api.nvim_win_set_cursor(0, {item.lnum, item.col-1})
    else
        vim.fn.setqflist({}, ' ', options)
        local default_qf_height = 10
        local height = math.min(default_qf_height, #options.items)
        vim.cmd('botright copen ' .. height)
    end
end

---Get the ListOpts which can be given to e.g. vim.lsp.buf.references or other
---lsp function.
---Filters the results placed in qf using the given `fun`.
---`fun` gets one argument `item`, see `:h setqflist-what`.
---It returns a boolean, indicating if the given `item` should be kept.
---@param fun function
---@return vim.lsp.ListOpts
function M.filter_lsp_items(fun)
    return {
        on_list = function(options)
            local items = {}
            for _, item in ipairs(options.items) do
                if fun(item) then
                    table.insert(items, item)
                end
            end
            options.items = items
            M.qf_mini(options)
        end
    }
end

return M
