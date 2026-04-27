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

---Query directive `(#trim! @cap PREFIX_BYTES SUFFIX_BYTES)`.
---Skips PREFIX_BYTES from the start (newline-aware: increments row, resets col)
---and SUFFIX_BYTES from the end (assumed not to cross a newline) of @cap, then
---sets metadata.range to a (row, col, byte) triple consistent across all three
---coordinates. `#offset!` does naive (row+drow, col+dcol) arithmetic and breaks
---when col_delta would push col past the end of a short line — e.g. injecting
---lua into `nvim -c 'lua\nCODE\n'` where the first line of the raw_string ends
---right after `'lua`.
---@param match table<integer, TSNode[]>
---@param _ integer pattern index (unused)
---@param source integer|string buffer or string
---@param pred any[]
---@param metadata vim.treesitter.query.TSMetadata
function M.trim_directive(match, _, source, pred, metadata)
    local capture_id = pred[2]
    local prefix_bytes = tonumber(pred[3]) or 0
    local suffix_bytes = tonumber(pred[4]) or 0
    local nodes = match[capture_id]
    if not nodes or #nodes == 0 then return end
    local node = nodes[1]
    local sr, sc, sb, er, ec, eb = node:range(true)
    local text = vim.treesitter.get_node_text(node, source)
    local new_sr, new_sc, new_sb = sr, sc, sb
    for i = 1, prefix_bytes do
        if text:byte(i) == 0x0a then
            new_sr = new_sr + 1
            new_sc = 0
        else
            new_sc = new_sc + 1
        end
        new_sb = new_sb + 1
    end
    if not metadata[capture_id] then metadata[capture_id] = {} end
    metadata[capture_id].range = {
        new_sr, new_sc, new_sb,
        er, ec - suffix_bytes, eb - suffix_bytes,
    }
end

return M
