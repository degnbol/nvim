local M = {}

---Iterate through parent treesitter nodes until reaching one of a given type.
---Useful to get e.g. the surrounding calling function, with type "call_expression".
---@param type string
---@return TSNode|nil
function M.get_parent(type)
    local node = vim.treesitter.get_node()
    while node ~= nil do
        if node:type() == type then
            return node
        end
        node = node:parent()
    end
end

return M
