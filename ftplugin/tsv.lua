#!/usr/bin/env lua
local util = require "utils/init"

local defaults = {
    maxwidth = 10, -- excl gap
    checklines = 100, -- max lines to look at to detect column widths
    checkevents = {"BufEnter", "BufWritePost"},
}

local widths -- set in BufEnter. Record of unhidden column widths.
-- abililty to have different column widths and used for indicating if column is currently hidden
local maxwidths = {}

local function getCommentChar()
    return vim.opt_local.commentstring:get():sub(1,1)
end
local function isComment(line, commentchar)
    return line:match("^[\t%s]*" .. commentchar)
end

---@row number 1-indexed
---@col number 1-indexed
---@return (c1, c2) 0-indexed, inclusive
local function getCellRange(row, col)
    local line = util.get_line(row-1)
    local fields = vim.split(line, '\t')
    local c = 0
    for i = 1, col-1 do c=c+#fields[i]+1 end
    return c, c + #fields[col]
end

-- a namespace per column so we can have id=row
local ns = {}

---@row 1-index
---@col 1-index
---@c 0-index
local function setExtmark(row, col, c)
    return vim.api.nvim_buf_set_extmark(0, ns[col], row-1, c, {
        id=row,
        virt_text={{"…", "NonText"}},
        virt_text_pos="overlay",
        -- virt_text_hide=true,
        -- hide if area is rm and "invalid" key in "details" of nvim_buf_get_extmarks
        invalidate=true,
    })
end

---@row 1-index
---@col 1-index
---@return bool found Whether the extmark was found
local function delExtmark(row, col)
    return vim.api.nvim_buf_del_extmark(0, ns[col], row)
end

local function getExtmark(row, col)
    local r, c, details = unpack(vim.api.nvim_buf_get_extmark_by_id(0, ns[col], row, {details=true}))
    return r, c, details.invalid
end

local hidden = {}

-- use with large maxwidth to unhide
local function hideCell(row, col)
    -- get range and text of full cell
    local c1, c2 = getCellRange(row, col)
    local text = util.get_text(row-1, c1, c2)
    -- if already hidden
    local hiddenText = hidden[col][row]
    if hiddenText ~= nil then
        -- reset buffer text
        text=text..hiddenText
        vim.api.nvim_buf_set_text(0, row-1, c1, row-1, c2, {text})
        c2=c2+#hiddenText
    end
    
    local cHidden = c1 + maxwidths[col]
    -- only hide if there is something to hide
    if cHidden >= c2 then
        hidden[col][row] = nil
        return delExtmark(row, col)
    end

    -- will make new or update the current one if already hidden
    setExtmark(row, col, cHidden+1) -- +1 to place in gap

    -- store hidden text
    hidden[col][row] = text:sub(cHidden - c1+1)
    -- remove text from buffer
    vim.api.nvim_buf_set_text(0, row-1, cHidden, row-1, c2, {})
end

-- update vartabstop
local function updateVartabstop()
    local vartabstop = {}
    for col, width in ipairs(widths) do
        -- use "hidden" to indicate if a column is currently hiding content 
        -- since we want to have maxwidths remember last used maxwidth for the 
        -- column
        if hidden[col] ~= nil and maxwidths[col] < width then
            -- 2 spaces gap
            table.insert(vartabstop, maxwidths[col]+2)
        else
            table.insert(vartabstop, width)
        end
    end
    vim.opt_local.vartabstop = vartabstop
end

local function hide(cols, maxwidth)
    -- 0 or less not valid maxwidth
    if maxwidth and maxwidth <= 0 then maxwidth = nil end

    local commentchar = getCommentChar()
    local commentlines = {}
    local lines = vim.api.nvim_buf_get_lines(0, 0, -1, true)
    for i, line in ipairs(lines) do
        if isComment(line, commentchar) then
            commentlines[i] = true
        end
    end

    for _, col in ipairs(cols) do
        -- mark column for hiding
        if ns[col] == nil then
            ns[col] = vim.api.nvim_create_namespace("hide-" .. col)
        end
        -- if no explicit maxwidth is provided then use one remembered for the 
        -- column or fallback to the default.
        maxwidths[col] = maxwidth or maxwidths[col] or defaults.maxwidth
        if hidden[col] == nil then hidden[col] = {} end
        -- for now hide all lines (except comments)
        for row = 1, #lines do
            if not commentlines[row] then
                hideCell(row, col)
            end
        end
    end
    updateVartabstop()
end

local function unhide(cols)
    for _, col in ipairs(cols) do
        if hidden[col] ~= nil then
            for row, hiddenText in pairs(hidden[col]) do
                local r, c, invalid = getExtmark(row, col)
                if not invalid then
                    vim.api.nvim_buf_set_text(0, r, c, r, c, {hiddenText})
                end
                delExtmark(row, col)
                hidden[col][row] = nil
            end
        end
        hidden[col] = nil
    end
    updateVartabstop()
end

---@r 0-indexed win column
---@c 0-indexed win column
---@return int 1-indexed for table column of the cursor.
local function getCol(r, c)
    local line = util.get_line(r)
    local _, nsub = line:sub(1,c):gsub('\t', '')
    return nsub+1
end

---@return int 1-indexed for current table column of the cursor.
local function getCurCol()
    local r, c = util.get_cursor()
    return getCol(r, c)
end

local function allCols()
    local cols = {}
    -- #widths indicates number of columns
    for i = 1, #widths do
        table.insert(cols, i)
    end
    return cols
end

---@return {int} 1-indexed for current table columns of visual selection range or cursor if normal mode.
local function getCurCols()
    local mode = util.get_mode()
    if mode == 'n' then
        return {getCurCol()}
    elseif mode == 'V' then
        return allCols()
    else
        local cols = {}
        local r1, c1, r2, c2 = util.get_visual_range()
        for i = getCol(r1-1, c1), getCol(r2-1, c2) do
            table.insert(cols, i)
        end
        util.gv()
        return cols
    end
end

vim.keymap.set( {'n', 'v'}, 'zc',
    function () hide(getCurCols(), vim.v.count) end,
    { desc="Hide current column(s)" }
)

vim.keymap.set( {'n', 'v'}, 'zo',
    function () unhide(getCurCols()) end,
    { desc="Unhide current column(s)" }
)

local grp = vim.api.nvim_create_augroup("hide", {clear=true})
-- this autocmd sets vartabstop on file save based on longest cell in each column.
vim.api.nvim_create_autocmd(defaults.checkevents, {
    buffer=0,
    group=grp,
    callback = function()
        local commentchar = getCommentChar()
        widths = {}

        local lines = vim.api.nvim_buf_get_lines(0, 0, defaults.checklines, false)
        for _, line in ipairs(lines) do
            -- ignore comment lines
            if not isComment(line, commentchar) then
                local fields = vim.split(line, '\t', true)
                -- when the line is shorter than or equal in length to a line seen so far
                for i = 1, math.min(#widths, #fields) do
                    widths[i] = math.max(widths[i], fields[i]:len())
                end
                -- when the line is longer than any line seen so far
                for i = #widths+1, #fields do
                    table.insert(widths, fields[i]:len())
                end
            end
        end

        -- min 2 visual spaces gap between columns
        for i, v in ipairs(widths) do
            widths[i] = v+2
        end

        -- apply
        vim.opt_local.vartabstop = widths
    end
})
-- unhide and rehide when saving as to always save the full text to file
vim.api.nvim_create_autocmd("BufWritePre", {
    buffer = 0,
    group = grp,
    callback = function () unhide(allCols()) end
})
vim.api.nvim_create_autocmd("BufWritePost", {
    buffer = 0,
    group = grp,
    callback = function ()
        -- rehide all
        for col, maxwidth in pairs(maxwidths) do
            hide({col}, maxwidth)
        end
    end
})


vim.keymap.set('n', ']]', function ()
    local r, c = unpack(vim.api.nvim_win_get_cursor(0))
    local line = vim.api.nvim_get_current_line()
    local nextTab = line:sub(c+1):match('()\t')
    if nextTab ~= nil then
        vim.api.nvim_win_set_cursor(0, {r, c+nextTab})
    end
end, { desc="Goto next start of cell" })
vim.keymap.set('n', '[[', function ()
    local r, c = unpack(vim.api.nvim_win_get_cursor(0))
    if c == 0 then return end
    local line = vim.api.nvim_get_current_line()
    local prevTab = line:sub(1, c-1):match('\t?()[^\t]*$')
    vim.api.nvim_win_set_cursor(0, {r, prevTab-1})
end, { desc="Goto previous start of cell" })
vim.keymap.set('n', '][', function ()
    local r, c = unpack(vim.api.nvim_win_get_cursor(0))
    local line = vim.api.nvim_get_current_line()
    -- TODO
end, { desc="Goto next end of cell" })
vim.keymap.set('n', '[]', function ()
    local r, c = unpack(vim.api.nvim_win_get_cursor(0))
    if c == 0 then return end
    local line = vim.api.nvim_get_current_line()
    -- TODO
end, { desc="Goto previous end of cell" })



