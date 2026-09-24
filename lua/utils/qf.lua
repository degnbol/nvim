local M = {}

--- Whether position `a` comes strictly before position `b`.
--- @param a integer[] `{ lnum, col }`
--- @param b integer[] `{ lnum, col }`, indexed like `a`
--- @return boolean
local function pos_lt(a, b)
    return a[1] < b[1] or (a[1] == b[1] and a[2] < b[2])
end

--- Whether quickfix item `a` sorts before item `b` by filename, then lnum, then col.
--- @param a table `:h setqflist-what` item with `filename`
--- @param b table `:h setqflist-what` item with `filename`
--- @return boolean
local function item_lt(a, b)
    return a.filename < b.filename or (a.filename == b.filename and pos_lt({ a.lnum, a.col }, { b.lnum, b.col }))
end

--- Quickfix item at the cursor of the current window.
--- @return table item `:h setqflist-what` item with the buffer name as `filename`,
--- the cursor as 1-based `lnum` and `col`, and the cursor line as `text`
local function cursor_item()
    local lnum, col = unpack(vim.api.nvim_win_get_cursor(0))
    return {
        filename = vim.api.nvim_buf_get_name(0),
        lnum = lnum,
        col = col + 1,
        text = vim.api.nvim_get_current_line(),
    }
end

--- Whether the cursor of the current window is in the range of a quickfix item.
--- The range is from (`lnum`, `col`) inclusive to (`end_lnum`, `end_col`)
--- exclusive, and always holds the start position, so zero-width ranges and
--- items without an end hold their start.
--- @param item table `:h setqflist-what` item. Its `filename` must equal the
--- buffer name exactly, so a symlinked path does not match.
--- @return boolean
function M.contains_cursor(item)
    if item.filename ~= vim.api.nvim_buf_get_name(0) then return false end
    local lnum, col = unpack(vim.api.nvim_win_get_cursor(0))
    local cursor = { lnum, col + 1 }
    local start = { item.lnum, item.col }
    local stop = { item.end_lnum or item.lnum, item.end_col or item.col }
    return not pos_lt(cursor, start) and (pos_lt(cursor, stop) or not pos_lt(start, cursor))
end

--- Load `what` as a new quickfix list whose current entry is the first item that
--- contains the cursor, then jump to the other item if there is exactly one.
--- With more, print their count and keep the cursor. The quickfix window is not
--- opened. With no other item, print a message and keep the previous list.
---
--- If no item contains the cursor, an item at the cursor is inserted before the
--- first item that sorts after it by filename, lnum, col. This keeps the order
--- only if `items` is already in that order.
--- @param what vim.fn.setqflist.what `:h setqflist-what` whose `items` all have
--- `filename`. Not mutated.
function M.jump_or_load(what)
    local items = assert(what.items)
    if not vim.iter(items):any(M.contains_cursor) then
        local here = cursor_item()
        local i_after = vim.iter(ipairs(items)):find(function(_, item) return item_lt(here, item) end)
        items = vim.list_extend({}, items)
        table.insert(items, i_after or #items + 1, here)
    end
    local i_self, i_other, n_other = nil, nil, 0
    for i, item in ipairs(items) do
        if M.contains_cursor(item) then
            i_self = i_self or i
        else
            i_other, n_other = i, n_other + 1
        end
    end
    if n_other == 0 then
        print("No other items.")
        return
    end
    vim.fn.setqflist({}, ' ', vim.tbl_extend('force', what, { items = items, idx = i_self }))
    if n_other == 1 then
        vim.cmd.cc(i_other)
    else
        print(("%s: %d items"):format(what.title or "Quickfix", n_other))
    end
end

return M
