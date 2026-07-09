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

---Advance `n` bytes forward from a start position through `text`, tracking
---(row, col, byte). Newline-aware: `\n` increments the row and resets the
---byte-column to 0. All coordinates are 0-indexed, matching `TSNode:range(true)`.
---@param text string
---@param n integer bytes to advance
---@param row integer
---@param col integer byte-column
---@param byte integer byte offset
---@return integer row
---@return integer col
---@return integer byte
local function advance_bytes(text, n, row, col, byte)
    for i = 1, n do
        if text:byte(i) == 0x0a then
            row = row + 1
            col = 0
        else
            col = col + 1
        end
        byte = byte + 1
    end
    return row, col, byte
end

---Query directive `(#trim! @cap PREFIX_BYTES SUFFIX_BYTES)`.
---Trims PREFIX_BYTES from the start and SUFFIX_BYTES from the end of @cap, then
---sets metadata.range to a (row, col, byte) triple consistent across all three
---coordinates. Both ends are computed by a newline-aware forward walk, so the
---range stays valid when the node spans lines — or when error recovery of an
---unterminated string ends the node at column 0, where `#offset!`'s naive
---(row+drow, col+dcol) arithmetic would emit a negative column and crash
---`set_included_ranges` ("Range value out of bounds"). A degenerate result
---(end before start) drops the injection.
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
    local sr, sc, sb = node:range(true)
    local text = vim.treesitter.get_node_text(node, source)
    local keep = #text - suffix_bytes
    if keep < prefix_bytes then return end -- nothing left after trimming
    local new_sr, new_sc, new_sb = advance_bytes(text, prefix_bytes, sr, sc, sb)
    local new_er, new_ec, new_eb = advance_bytes(text, keep, sr, sc, sb)
    if not metadata[capture_id] then metadata[capture_id] = {} end
    metadata[capture_id].range = { new_sr, new_sc, new_sb, new_er, new_ec, new_eb }
end

---Query directive `(#inject-by-ext! @dest)` — set `injection.language` from the
---file extension of @dest (a heredoc redirect destination). Resolves extension
---→ filetype (`vim.filetype.match`) → parser language
---(`vim.treesitter.language.get_lang`, which honours the sh/bash → zsh
---registrations). Surrounding quotes on the destination are stripped first.
---Unknown extensions leave `injection.language` unset, so the capture is
---ignored and the base heredoc-tag injection (`<<LUA ... LUA`) still applies.
---@param match table<integer, TSNode[]>
---@param _ integer pattern index (unused)
---@param source integer|string buffer or string
---@param pred any[]
---@param metadata vim.treesitter.query.TSMetadata
function M.inject_by_ext_directive(match, _, source, pred, metadata)
    local capture_id = pred[2]
    local nodes = match[capture_id]
    if not nodes or #nodes == 0 then return end
    local dest = vim.treesitter.get_node_text(nodes[1], source)
    dest = dest:gsub("^['\"]", ""):gsub("['\"]$", "")
    local ft = vim.filetype.match({ filename = dest })
    if not ft then return end
    metadata["injection.language"] = vim.treesitter.language.get_lang(ft) or ft
end

---Query predicate `(#any-basename-of? @cap "name" ...)` — like `#any-of?` but
---compares the *basename* of @cap, so a path-prefixed interpreter such as
---`.venv/bin/python` or `/usr/bin/python3` matches the bare name (`python`,
---`python3`). Any leading `dir/` components are stripped before the test.
---@param match table<integer, TSNode[]>
---@param _ integer pattern index (unused)
---@param source integer|string buffer or string
---@param pred any[]
---@return boolean
function M.any_basename_of(match, _, source, pred)
    local nodes = match[pred[2]]
    if not nodes or #nodes == 0 then return true end
    for _, node in ipairs(nodes) do
        local basename = vim.treesitter.get_node_text(node, source):gsub(".*/", "")
        for i = 3, #pred do
            if basename == pred[i] then return true end
        end
    end
    return false
end

---Query directive `(#head! @cap N)` — narrow @cap to its first N bytes.
---Assumes the first N bytes do not span a newline.
---@param match table<integer, TSNode[]>
---@param _ integer pattern index (unused)
---@param _source integer|string buffer or string (unused)
---@param pred any[]
---@param metadata vim.treesitter.query.TSMetadata
function M.head_directive(match, _, _source, pred, metadata)
    local capture_id = pred[2]
    local n = tonumber(pred[3]) or 1
    local nodes = match[capture_id]
    if not nodes or #nodes == 0 then return end
    local sr, sc, sb = nodes[1]:range(true)
    if not metadata[capture_id] then metadata[capture_id] = {} end
    metadata[capture_id].range = { sr, sc, sb, sr, sc + n, sb + n }
end

---Query directive `(#tail! @cap N)` — narrow @cap to its last N bytes.
---Assumes the last N bytes do not span a newline.
---@param match table<integer, TSNode[]>
---@param _ integer pattern index (unused)
---@param _source integer|string buffer or string (unused)
---@param pred any[]
---@param metadata vim.treesitter.query.TSMetadata
function M.tail_directive(match, _, _source, pred, metadata)
    local capture_id = pred[2]
    local n = tonumber(pred[3]) or 1
    local nodes = match[capture_id]
    if not nodes or #nodes == 0 then return end
    local _, _, _, er, ec, eb = nodes[1]:range(true)
    if not metadata[capture_id] then metadata[capture_id] = {} end
    metadata[capture_id].range = { er, ec - n, eb - n, er, ec, eb }
end

return M
