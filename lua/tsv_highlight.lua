-- Highlights of TSV values and comment lines, as extmarks set once per line
-- rather than by a treesitter query, which runs again on every redraw.
local M = {}

local util = require "utils/init"

local ns = vim.api.nvim_create_namespace("tsv_highlight")

---Rows whose marks are up to date, per attached buffer.
---@type table<integer, table<integer, true>>
local done_rows = {}

---True for unsigned digits with a decimal point and digits on at least one side.
---@param text string
---@return boolean
local function is_decimal(text)
    return text:match("^%d+%.%d*$") ~= nil or text:match("^%.%d+$") ~= nil
end

---The highlight group of a TSV cell's value.
---@param text string the cell, without its tabs
---@return string|nil group `@number`, `@number.float` or `@boolean`; nil for text
function M.cell_group(text)
    if text == "true" or text == "false" then return "@boolean" end
    if text:match("^0[xX]%x+$") then return "@number" end
    local unsigned = text:gsub("^[-+]", "")
    local mantissa = unsigned:match("^(.-)[eE][-+]?%d+$")
    if mantissa then
        if mantissa:match("^%d+$") or is_decimal(mantissa) then return "@number.float" end
    elseif unsigned:match("^%d+$") then
        return "@number"
    elseif is_decimal(unsigned) then
        return "@number.float"
    end
end

---@alias tsv_highlight.Mark { [1]: integer, [2]: integer, [3]: string } col, end_col, group

---The highlights of a TSV line.
---@param line string
---@param leader string what a comment opens with, e.g. "#"
---@return tsv_highlight.Mark[] marks 0-based, end-exclusive byte ranges
function M.line_marks(line, leader)
    if util.is_comment_line(line, leader) then return { { 0, #line, "@comment" } } end
    local marks = {}
    local col = 0
    for cell in vim.gsplit(line, "\t", { plain = true }) do
        local group = M.cell_group(cell)
        if group then marks[#marks + 1] = { col, col + #cell, group } end
        col = col + #cell + 1
    end
    return marks
end

vim.api.nvim_set_decoration_provider(ns, {
    on_win = function(_, _, buf, toprow, botrow)
        local done = done_rows[buf]
        if not done then return false end
        local leader = util.comment_leader(buf)
        for row = toprow, math.min(botrow, vim.api.nvim_buf_line_count(buf) - 1) do
            if not done[row] then
                vim.api.nvim_buf_clear_namespace(buf, ns, row, row + 1)
                local line = vim.api.nvim_buf_get_lines(buf, row, row + 1, true)[1]
                for _, mark in ipairs(M.line_marks(line, leader)) do
                    vim.api.nvim_buf_set_extmark(buf, ns, row, mark[1], {
                        end_col = mark[2],
                        hl_group = mark[3],
                        priority = vim.hl.priorities.treesitter,
                    })
                end
                done[row] = true
            end
        end
        return false
    end,
})

vim.api.nvim_create_autocmd("OptionSet", {
    pattern = "commentstring",
    group = vim.api.nvim_create_augroup("tsv_highlight", {}),
    callback = function()
        -- OptionSet gives no <abuf>; the option is set with its buffer current.
        local buf = vim.api.nvim_get_current_buf()
        if not done_rows[buf] then return end
        done_rows[buf] = {}
        vim.api.nvim__redraw { buf = buf, valid = false }
    end,
})

---Buffers with a live `nvim_buf_attach` callback. It outlives `detach` until
---the next change, so `attach` before then must not add a second one.
---@type table<integer, true>
local listening = {}

---Highlight the values and comment lines of a TSV buffer. A line is classified
---when a window first draws it, and again after it changes.
---@param buf integer 0 for the current buffer
function M.attach(buf)
    if buf == 0 then buf = vim.api.nvim_get_current_buf() end
    done_rows[buf] = {}
    if listening[buf] then return end
    listening[buf] = true
    vim.api.nvim_buf_attach(buf, false, {
        on_lines = function(_, _, _, first, last, new_last)
            local done = done_rows[buf]
            if not done then
                listening[buf] = nil
                return true
            end
            if new_last ~= last then
                -- The marks moved with their lines, the rows did not.
                done_rows[buf] = {}
                return
            end
            -- Through the row after the change: the start of a replaced
            -- line's mark moves to the next row.
            for row = first, new_last do done[row] = nil end
        end,
        on_reload = function() done_rows[buf] = {} end,
        on_detach = function()
            done_rows[buf] = nil
            listening[buf] = nil
        end,
    })
end

---Remove the highlights `attach` gave a buffer.
---@param buf integer 0 for the current buffer
function M.detach(buf)
    if buf == 0 then buf = vim.api.nvim_get_current_buf() end
    done_rows[buf] = nil
    vim.api.nvim_buf_clear_namespace(buf, ns, 0, -1)
end

return M
