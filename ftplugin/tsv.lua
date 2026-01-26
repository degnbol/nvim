local util = require "utils/init"

-- NOTE Limitations:
-- - Assumes tab is field sep and newline is record sep, no multiline cells.
-- - Hiding cuts through multi-byte chars like €, might be a TODO.
-- - Undo/redo: u and <C-r> are mapped to smart versions that unhide before
--   the operation and re-hide after, making hide transparent. Use zo/za to
--   actually unhide. Users can remap u/<C-r> or use <Plug>(TsvUndo/Redo).

local defaults = {
    maxwidth = 10,     -- excl gap
    checklines = 1000, -- max lines to look at to detect column widths
    checkevents = { "BufEnter", "FileType", "BufWritePost" },
}

-- Buffer-local state accessors (stored in vim.b.tsv_*)
-- vim.b converts sparse tables by padding with vim.NIL, which is truthy but not indexable.
-- We must convert vim.NIL back to nil for safe usage.
local function denilify(t)
    if t == nil then return {} end
    local result = {}
    for k, v in pairs(t) do
        if v ~= vim.NIL then
            if type(v) == "table" then
                result[k] = denilify(v)
            else
                result[k] = v
            end
        end
    end
    return result
end
local function get_widths() return denilify(vim.b.tsv_widths) end
local function set_widths(w) vim.b.tsv_widths = w end
local function get_maxwidths() return denilify(vim.b.tsv_maxwidths) end
local function set_maxwidths(m) vim.b.tsv_maxwidths = m end
local function get_hidden() return denilify(vim.b.tsv_hidden) end
local function set_hidden(h) vim.b.tsv_hidden = h end

local function getCommentChar()
    return vim.opt_local.commentstring:get():sub(1, 1)
end
local function isComment(line, commentchar)
    return commentchar and commentchar ~= "" and line:match("^[\t%s]*" .. commentchar)
end

--- Takes a count, to override default number of lines to check for detecting column max lengths.
local function updateWidths()
    local commentchar = getCommentChar()
    local widths = {}

    local checklines = vim.v.count
    if checklines == 0 then checklines = defaults.checklines end
    local lines = vim.api.nvim_buf_get_lines(0, 0, checklines, false)

    for _, line in ipairs(lines) do
        -- ignore comment lines
        if not isComment(line, commentchar) then
            local fields = vim.split(line, '\t', true)
            -- when the line is shorter than or equal in length to a line seen so far
            for i = 1, math.min(#widths, #fields) do
                widths[i] = math.max(widths[i], fields[i]:len())
            end
            -- when the line is longer than any line seen so far
            for i = #widths + 1, #fields do
                table.insert(widths, fields[i]:len())
            end
        end
    end

    set_widths(widths)

    -- min 2 visual spaces gap between columns
    local vartabstop = {}
    for i, width in ipairs(widths) do
        vartabstop[i] = width + 2
    end
    vim.opt_local.vartabstop = vartabstop
end

--- Not adding hidden text lengths.
---@param row integer 1-indexed table row
---@param col integer 1-indexed table column
---@return integer? c1 0-indexed, inclusive or nil if the cell doesn't exist.
---@return integer? c2 0-indexed, inclusive, i.e. location of following tab/newline
local function getCellRange(row, col)
    local line = util.get_line(row - 1)
    local fields = vim.split(line, '\t')
    if #fields < col then return end
    local c = 0
    for i = 1, col - 1 do c = c + #fields[i] + 1 end
    return c, c + #fields[col]
end

-- a namespace per column so we can have id=row
local ns = {}

---@param row integer 1-indexed
---@param col integer 1-indexed
---@param c integer 0-indexed
local function setExtmark(row, col, c)
    return vim.api.nvim_buf_set_extmark(0, ns[col], row - 1, c, {
        id = row,
        virt_text = { { "…", "NonText" } },
        virt_text_pos = "overlay",
        -- virt_text_hide=true,
        -- hide if area is rm and "invalid" key in "details" of nvim_buf_get_extmarks
        invalidate = true,
    })
end

---@param row integer 1-indexed
---@param col integer 1-indexed
---@return boolean found Whether the extmark was found
local function delExtmark(row, col)
    return vim.api.nvim_buf_del_extmark(0, ns[col], row)
end

---@param row integer 1-indexed
---@param col integer 1-indexed
---@return integer r 0-indexed
---@return integer c 0-indexed
---@return boolean? invalid Whether the extmark is invalid
local function getExtmark(row, col)
    local r, c, details = unpack(vim.api.nvim_buf_get_extmark_by_id(0, ns[col], row, { details = true }))
    return r, c, details.invalid
end

local function updateVartabstop()
    local widths = get_widths()
    local maxwidths = get_maxwidths()
    local hidden = get_hidden()
    local vartabstop = {}
    for col, width in ipairs(widths) do
        -- use "hidden" to indicate if a column is currently hiding content
        -- since we want to have maxwidths remember last used maxwidth for the
        -- column
        if hidden[col] ~= nil and maxwidths[col] < width then
            -- 2 spaces gap
            table.insert(vartabstop, maxwidths[col] + 2)
        else
            -- 2 spaces gap
            table.insert(vartabstop, width + 2)
        end
    end
    vim.opt_local.vartabstop = vartabstop
end


local function unhideCell(row, col, hiddenText)
    local r, c, invalid = getExtmark(row, col)
    if not invalid then
        vim.api.nvim_buf_set_text(0, r, c, r, c, { hiddenText })
    end
    delExtmark(row, col)
    local hidden = get_hidden()
    if hidden[col] then hidden[col][row] = nil end
    set_hidden(hidden)
    return r, c
end

-- use with large maxwidth to unhide
local function hideCell(row, col)
    local hidden = get_hidden()
    local maxwidths = get_maxwidths()
    -- get range and text of full cell
    local c1, c2 = getCellRange(row, col)
    -- in case the cell doesn't exist
    if c1 == nil then return end
    local text = util.get_text(row - 1, c1, c2)
    -- if already hidden
    local hiddenText = hidden[col] and hidden[col][row]
    if hiddenText ~= nil then
        -- reset
        local rUnhidden, _ = unhideCell(row, col, hiddenText)
        -- if lines have been added or removed the rows may be shifted
        if rUnhidden ~= row - 1 then
            -- rehide the unhidden cell with corrected id
            hideCell(rUnhidden + 1, col)
        else
            -- recalculate
            c1, c2 = getCellRange(row, col)
            text = util.get_text(row - 1, c1, c2)
        end
        -- re-fetch hidden after unhideCell modified it
        hidden = get_hidden()
    end

    local cHidden = c1 + maxwidths[col]
    -- only hide if there is something to hide
    if cHidden >= c2 then
        if hidden[col] then hidden[col][row] = nil end
        set_hidden(hidden)
        return delExtmark(row, col)
    end

    -- will make new or update the current one if already hidden
    setExtmark(row, col, cHidden + 1) -- +1 to place in gap

    -- store hidden text
    if not hidden[col] then hidden[col] = {} end
    hidden[col][row] = text:sub(cHidden - c1 + 1)
    set_hidden(hidden)
    -- remove text from buffer
    vim.api.nvim_buf_set_text(0, row - 1, cHidden, row - 1, c2, {})
end

local function unhide(cols)
    local was_modified = vim.bo.modified
    local hidden = get_hidden()
    for _, col in ipairs(cols) do
        if hidden[col] ~= nil then
            for row, hiddenText in pairs(hidden[col]) do
                unhideCell(row, col, hiddenText)
            end
        end
        hidden[col] = nil
    end
    set_hidden(hidden)
    updateVartabstop()
    vim.bo.modified = was_modified
end


local function hide(cols, maxwidth)
    local was_modified = vim.bo.modified
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

    local maxwidths = get_maxwidths()
    local hidden = get_hidden()
    for _, col in ipairs(cols) do
        -- mark column for hiding
        if ns[col] == nil then
            ns[col] = vim.api.nvim_create_namespace("hide-" .. col)
        end
        -- if no explicit maxwidth is provided then use one remembered for the
        -- column or fallback to the default.
        maxwidths[col] = maxwidth or maxwidths[col] or defaults.maxwidth
        if hidden[col] == nil then hidden[col] = {} end
    end
    set_maxwidths(maxwidths)
    set_hidden(hidden)
    -- for now hide all lines (except comments)
    for _, col in ipairs(cols) do
        for row = 1, #lines do
            if not commentlines[row] then
                hideCell(row, col)
            end
        end
    end
    updateVartabstop()
    vim.bo.modified = was_modified
end


---@param r integer 0-indexed line number
---@param c integer? 0-indexed character column, or nil to get number of columns (max)
---@return integer col 1-indexed for table column of the cursor.
local function getCol(r, c)
    local line = util.get_line(r)
    local _, nsub = line:sub(1, c):gsub('\t', '')
    return nsub + 1
end

---Get current cell row and column where the cursor is.
---@return integer 1-indexed for current table column.
local function getCurCol()
    local r, c = util.get_cursor()
    return getCol(r, c)
end
---@return integer r  0-indexed line number
---@return integer c1 0-indexed char column
---@return integer c2 0-indexed char column
local function getCurCellRange()
    local r, c = util.get_cursor()
    local col = getCol(r, c)
    local c1, c2 = getCellRange(r + 1, col)
    return r, c1, c2
end

local function allCols()
    local cols = {}
    local widths = get_widths()
    -- #widths indicates number of columns
    for i = 1, #widths do
        table.insert(cols, i)
    end
    return cols
end

---@return integer[] cols 1-indexed columns of visual selection range or cursor if normal mode
local function getCurCols()
    local mode = util.get_mode()
    if mode == 'n' then
        return { getCurCol() }
    elseif mode == 'V' then
        return allCols()
    else
        local cols = {}
        local r1, c1, r2, c2 = util.get_visual_range()
        for i = getCol(r1 - 1, c1), getCol(r2 - 1, c2) do
            table.insert(cols, i)
        end
        util.gv()
        return cols
    end
end

vim.keymap.set({ 'n', 'v' }, 'zc',
    function() hide(getCurCols(), vim.v.count) end,
    { desc = "Hide current column(s)" }
)

vim.keymap.set({ 'n', 'v' }, 'zo',
    function() unhide(getCurCols()) end,
    { desc = "Unhide current column(s)" }
)

vim.keymap.set({ 'n', 'v' }, 'za',
    function()
        local cols = getCurCols()
        local hidden = get_hidden()
        local to_hide = {}
        local to_unhide = {}
        for _, col in ipairs(cols) do
            if hidden[col] ~= nil then
                table.insert(to_unhide, col)
            else
                table.insert(to_hide, col)
            end
        end
        if #to_unhide > 0 then unhide(to_unhide) end
        if #to_hide > 0 then hide(to_hide, vim.v.count) end
    end,
    { desc = "Toggle hiding current column(s)" }
)

-- vim.keymap.set( {'n', 'v'}, '<Plug>TsvHideMore',
vim.keymap.set({ 'n', 'v' }, '{', function()
        local cols = getCurCols()
        local hidden = get_hidden()
        local maxwidths = get_maxwidths()
        local widths = get_widths()
        for _, col in ipairs(cols) do
            -- use maxwidth if currently hiding
            local width = hidden[col] and maxwidths[col] or widths[col]
            hide({ col }, math.max(1, width - vim.v.count1))
        end
    end,
    { desc = "Hide more of current column(s)", buffer = true }
)

-- vim.keymap.set( {'n', 'v'}, '<Plug>TsvHideMore',
vim.keymap.set({ 'n', 'v' }, '}', function()
        local cols = getCurCols()
        local hidden = get_hidden()
        local maxwidths = get_maxwidths()
        local widths = get_widths()
        for _, col in ipairs(cols) do
            -- hiding less is only relevant if currently hiding
            if hidden[col] ~= nil then
                local maxwidth = maxwidths[col] + vim.v.count1
                if maxwidth < widths[col] then
                    hide({ col }, maxwidth)
                else
                    unhide({ col })
                end
            end
        end
    end,
    { desc = "Hide less of current column(s)", buffer = true }
)

local grp = vim.api.nvim_create_augroup("hide", { clear = true })
-- this autocmd sets vartabstop on file save based on longest cell in each column.
vim.api.nvim_create_autocmd(defaults.checkevents, {
    buffer = 0, -- since the extension is not just .tsv but can also be .tab or .bed as defined in ftdetect/tsv.vim
    group = grp,
    callback = updateWidths,
})

-- also set a manual call in case we don't want to save
vim.keymap.set('n', '<localleader>a', updateWidths, { desc = "Align columns" })

-- unhide and rehide when saving as to always save the full text to file
local cols_hidden_before_write = {}
vim.api.nvim_create_autocmd("BufWritePre", {
    buffer = 0,
    group = grp,
    callback = function()
        -- Remember which columns are hidden before we unhide them
        cols_hidden_before_write = {}
        local hidden = get_hidden()
        for col, _ in pairs(hidden) do
            table.insert(cols_hidden_before_write, col)
        end
        unhide(allCols())
    end
})
vim.api.nvim_create_autocmd("BufWritePost", {
    buffer = 0,
    group = grp,
    callback = function()
        -- Re-hide only columns that were hidden before save
        local maxwidths = get_maxwidths()
        for _, col in ipairs(cols_hidden_before_write) do
            if maxwidths[col] then
                hide({ col }, maxwidths[col])
            end
        end
        cols_hidden_before_write = {}
    end
})

-- FIXME: this causes StackOverflow in large file, e.g. yank and paste a line with hidden cells in
-- /Users/cdmadsen/Documents/Topology/Chromatin/Pub/Hinch_2019/GSE124991_DMC1_and_H3K4me3_B6CASTF1.PRDM9hc.txt
-- yank hidden text as well
vim.api.nvim_create_autocmd("TextYankPost", {
    buffer = 0,
    group = grp,
    callback = function()
        local hidden = get_hidden()
        if vim.tbl_isempty(hidden) then return end
        local lines = vim.v.event.regcontents
        local linewise = vim.v.event.regtype == 'V'
        -- get range of yank 0-ind
        local r1, c1, r2, _ = util.last_changeyank_range()
        local row1 = r1 + 1
        local row2 = r2 + 1
        -- shift for each line. zero for other lines until we start adding in hidden text
        local cShift = { -c1 }
        for row = row1, row2 do
            for col, h in pairs(hidden) do
                local text = h[row]
                if text ~= nil then
                    local r, c, invalid = getExtmark(row, col)
                    if not invalid then
                        local i = r - r1 + 1 -- 1-ind
                        local line = lines[i]
                        if line ~= nil then
                            local shift = cShift[i] or 0
                            c = c + shift
                            if linewise or (c >= 0 and (c < #line or row < row2)) then
                                lines[i] = line:sub(1, c) .. text .. line:sub(c + 1)
                                cShift[i] = shift + #text
                            end
                        end
                    end
                end
            end
        end
        -- update register
        local regname = vim.v.event.regname
        if regname == "" then regname = "*" end
        vim.fn.setreg(regname, table.concat(lines, '\n'), vim.v.event.regtype)
    end
})

--- Check if any hidden column has cells longer than maxwidth (hidden text restored).
local function hiddenTextRestored()
    local hidden = get_hidden()
    local maxwidths = get_maxwidths()
    for col, rows in pairs(hidden) do
        if type(rows) == "table" then
            local maxw = maxwidths[col] or defaults.maxwidth
            for row, _ in pairs(rows) do
                if type(row) == "number" then
                    local c1, c2 = getCellRange(row, col)
                    if c1 and (c2 - c1) > maxw then
                        return true
                    end
                end
            end
        end
    end
    return false
end

--- Smart undo: if undo restores hidden text (undoes a hide), skip it.
local function smartUndo()
    local hidden = get_hidden()
    local maxwidths = get_maxwidths()

    -- Collect which columns are hidden
    local hidden_cols = {}
    for col, rows in pairs(hidden) do
        if type(rows) == "table" and not vim.tbl_isempty(rows) then
            hidden_cols[col] = maxwidths[col] or defaults.maxwidth
        end
    end

    if vim.tbl_isempty(hidden_cols) then
        vim.cmd("silent! undo")
        return
    end

    -- Do undo
    local seq_before = vim.fn.undotree().seq_cur
    vim.cmd("silent! undo")
    local seq_after = vim.fn.undotree().seq_cur

    -- If nothing changed, we're at oldest state
    if seq_before == seq_after then
        return
    end

    -- Check if undo restored hidden text (undid a hide operation)
    if hiddenTextRestored() then
        -- Try another undo to skip the hide
        local seq_before2 = vim.fn.undotree().seq_cur
        vim.cmd("silent! undo")
        local seq_after2 = vim.fn.undotree().seq_cur

        if seq_before2 ~= seq_after2 then
            -- Second undo worked, clear old hidden state and re-hide fresh
            for col, _ in pairs(hidden_cols) do
                hidden[col] = nil
                if ns[col] then
                    vim.api.nvim_buf_clear_namespace(0, ns[col], 0, -1)
                end
            end
            set_hidden(hidden)
            for col, maxw in pairs(hidden_cols) do
                hide({ col }, maxw)
            end
        else
            -- At oldest edit, redo to restore hidden state
            vim.cmd("silent! redo")
        end
    end
end

--- Smart redo: if redo redoes a hide, skip it.
local function smartRedo()
    local hidden = get_hidden()
    local maxwidths = get_maxwidths()

    -- Collect which columns are hidden
    local hidden_cols = {}
    for col, rows in pairs(hidden) do
        if type(rows) == "table" and not vim.tbl_isempty(rows) then
            hidden_cols[col] = maxwidths[col] or defaults.maxwidth
        end
    end

    if vim.tbl_isempty(hidden_cols) then
        vim.cmd("silent! redo")
        return
    end

    -- Do redo
    local seq_before = vim.fn.undotree().seq_cur
    vim.cmd("silent! redo")
    local seq_after = vim.fn.undotree().seq_cur

    -- If nothing changed, we're at newest state
    if seq_before == seq_after then
        return
    end

    -- Check if redo made cells shorter (redid a hide operation)
    -- If cells are now short but we have hidden state, the redo was a hide
    -- In this case, do another redo to get the next edit
    local all_short = true
    for col, maxw in pairs(hidden_cols) do
        local rows = hidden[col]
        if type(rows) == "table" then
            for row, _ in pairs(rows) do
                if type(row) == "number" then
                    local c1, c2 = getCellRange(row, col)
                    if c1 and (c2 - c1) > maxw then
                        all_short = false
                        break
                    end
                end
            end
        end
        if not all_short then break end
    end

    -- If all cells are short after redo, we redid a hide - try another redo
    if all_short then
        vim.cmd("silent! redo")
    end
end

-- Plugin mappings for smart undo/redo
vim.keymap.set('n', '<Plug>(TsvUndo)', smartUndo, { buffer = true })
vim.keymap.set('n', '<Plug>(TsvRedo)', smartRedo, { buffer = true })

-- Default mappings (users can override in their config)
vim.keymap.set('n', 'u', '<Plug>(TsvUndo)', { buffer = true, remap = true })
vim.keymap.set('n', '<C-r>', '<Plug>(TsvRedo)', { buffer = true, remap = true })


vim.keymap.set('n', ']]', function()
    local r, c = util.get_cursor()
    local col = getCol(r, c)
    local nCol = getCol(r, nil)
    local c1, _ = getCellRange(r + 1, math.min(nCol, col + vim.v.count1))
    util.set_cursor(r, c1)
end, { desc = "Goto next start of cell", buffer = true })
vim.keymap.set('n', '[[', function()
    local r, c = util.get_cursor()
    local col = getCol(r, c)
    local c1, _ = getCellRange(r + 1, math.max(1, col - vim.v.count1))
    util.set_cursor(r, c1)
end, { desc = "Goto previous start of cell", buffer = true })
vim.keymap.set('n', '][', function()
    local r, c = util.get_cursor()
    local col = getCol(r, c)
    local nCol = getCol(r, nil)
    local c1, c2 = getCellRange(r + 1, math.min(nCol, col + vim.v.count1))
    util.set_cursor(r, math.max(c2 - 1, c1))
end, { desc = "Goto next end of cell" })
vim.keymap.set('n', '[]', function()
    local r, c = util.get_cursor()
    local col = getCol(r, c)
    local c1, c2 = getCellRange(r + 1, math.max(1, col - vim.v.count1))
    util.set_cursor(r, math.max(c2 - 1, c1))
end, { desc = "Goto previous end of cell", buffer = true })

vim.keymap.set('x', 'ax', function()
    local mode = util.get_mode()
    local r, c1, c2 = getCurCellRange()
    if c1 == nil then return end
    util.end_visual()
    util.set_cursor(r, c1)
    util.set_mode(mode)
    util.set_cursor(r, c2)
end, { desc = "In cell" })
vim.keymap.set('x', 'ix', function()
    local mode = util.get_mode()
    local r, c1, c2 = getCurCellRange()
    if c1 == nil or c1 >= c2 then return end
    util.end_visual()
    util.set_cursor(r, c1)
    util.set_mode(mode)
    util.set_cursor(r, c2 - 1)
end, { desc = "In cell" })
vim.keymap.set('o', 'ax', function()
    local r, c1, c2 = getCurCellRange()
    if c1 == nil then return end
    util.end_visual()
    util.set_cursor(r, c1)
    util.set_mode('v')
    util.set_cursor(r, c2)
end, { desc = "In cell" })
vim.keymap.set('o', 'ix', function()
    local r, c1, c2 = getCurCellRange()
    if c1 == nil or c1 >= c2 then return end
    util.end_visual()
    util.set_cursor(r, c1)
    util.set_mode('v')
    util.set_cursor(r, c2 - 1)
end, { desc = "In cell" })




-- PERSISTENT HEADER

local getopt = vim.api.nvim_get_option_value
local setopt = vim.api.nvim_set_option_value

local grp = vim.api.nvim_create_augroup("TSV-header", { clear = true })
local grp_float = vim.api.nvim_create_augroup("TSV-header-float", { clear = true })
local winid, winid_float

local function close_header()
    vim.api.nvim_win_close(winid_float, false)
    vim.api.nvim_clear_autocmds({ group = grp_float })
    winid_float = nil
end

---Open a float showing table header when scrolled away from top of file.
---@param row integer? 1-indexed row of header, default=1
---@param height integer? height of header, default=1 + number of comment lines.
local function open_header(row, height)
    row = row or 1
    if height == nil then
        -- height detected as number of comment lines plus 1
        local commentchar = getCommentChar()
        -- consider max 20 lines
        local lines = vim.api.nvim_buf_get_lines(0, 0, 20, false)
        for i, line in ipairs(lines) do
            if not isComment(line, commentchar) then
                height = i
                break
            end
        end
    end

    -- open minimal float
    winid_float = vim.api.nvim_open_win(0, false, {
        relative = 'win',
        width = vim.api.nvim_win_get_width(0),
        height = height,
        row = 0,
        col = 0,
        bufpos = { row - 1, 0 },
        focusable = false,
        noautocmd = true,
    })

    winid = vim.api.nvim_get_current_win()
    local bufid = vim.api.nvim_win_get_buf(winid)

    -- scrollbind horizontally
    setopt("scrollbind", true, { win = winid })
    setopt("scrollbind", true, { win = winid_float })
    setopt("scrollopt", "hor", {})
    -- if we have multiple header rows the wrong lines will be shown if scrolloff > 0
    setopt("scrolloff", 0, { win = winid_float })
    -- WARNING: hacky solution here
    -- for multiline header when some lines aren't visible the float will show the wrong lines if we simply move its cursor to the final header line.
    -- Therefore we have to scroll the view, but there is no lua function for this currently.
    -- We scroll by moving cursor twice, first to the first row then to the last row.
    -- There are other ways to scroll, e.g. set win with api then call zb or <C-y> but they end up moving the main buffer window.
    vim.api.nvim_win_set_cursor(winid_float, { row, vim.api.nvim_win_get_cursor(winid_float)[2] })
    vim.api.nvim_win_set_cursor(winid_float, { row + height - 1, vim.api.nvim_win_get_cursor(winid_float)[2] })
    -- Due to moving the cursor to the first row for the sake of scrolling we now have to potentially fix wrong horizontal scroll, if the first row was shorter than the last.
    -- We do this by triggering a response to update the horizontal scroll but the trick here is we don't actually want to change the view.
    -- I tried with `vim.cmd.doautocmd "WinScrolled"` and tried setting vim.v.event but didn't trigger any update to header scroll.
    -- I trigger scroll event by moving down and up. Was the best option among left/right or up then down for no movement of view in edge cases.
    -- The cursor may still move so mark h, then jump to mark h makes sure the cursor doesn't move.
    util.press("mh<C-e><C-y>`h", { remap = false })

    -- map CursorLineNr to regular LineNr since there is no cursor in float
    setopt("winhighlight", "CursorLineNr:LineNr", { win = winid_float })
    -- hide LineNrs since they are inconsistent with relativenumber=on
    vim.api.nvim_win_set_hl_ns(winid_float, 1)
    vim.defer_fn(function()
        local bg = vim.api.nvim_get_hl(0, { name = "NormalFloat", link = false })["bg"]
        vim.api.nvim_set_hl(1, "CursorLineNr", { fg = bg })
        vim.api.nvim_set_hl(1, "LineNr", { fg = bg })
    end, 500)

    vim.api.nvim_create_autocmd("OptionSet", {
        pattern = { "*number", "signcolumn" },
        group = grp_float,
        callback = function(ev)
            local setto
            if ev.match:match(".*number$") then
                -- relativenumber makes little sense for the header since it would
                -- show relative to non-existing cursor in the float, plus it
                -- shifts the text incorrectly since there might only be 1 column
                -- for numbering in header while there is usually two or more in a
                -- normally sized buffer.
                ev.match = "number"
                setto = getopt("number", { win = winid }) or getopt("relativenumber", { win = winid })
            else
                setto = getopt(ev.match, { win = winid })
            end
            setopt(ev.match, setto, { win = winid_float })
        end
    })
    vim.api.nvim_create_autocmd("BufLeave", {
        buffer = bufid,
        group = grp_float,
        callback = function()
            close_header()
            -- open again if we come back
            vim.api.nvim_create_autocmd("BufEnter", {
                buffer = bufid,
                group = grp,
                callback = function()
                    open_header(row, height)
                    -- remove this aucmd after running it, since we make a new in the open_header call
                    return true
                end
            })
        end
    })
end

-- TODO: make config and have an if statement whether we want header by default
if false then
    vim.api.nvim_create_autocmd("BufEnter", {
        pattern = "*.tsv",
        group = grp,
        callback = open_header,
    })
end

vim.keymap.set('n', '<LocalLeader>h', function()
    if vim.v.count == 0 then
        -- toggle
        if winid_float == nil then
            open_header()
        else
            close_header()
        end
    else
        -- set size using count
        if winid_float == nil then
            open_header(1, vim.v.count)
        else
            vim.api.nvim_win_set_config(winid_float, { height = vim.v.count })
        end
    end
end, { buffer = true, desc = "Toggle header or set its size with count" })
vim.keymap.set('v', '<LocalLeader>h', function()
    local r1, _, r2, _ = util.get_visual_range()
    -- close then reopen if already open
    if winid_float ~= nil then close_header() end
    open_header(r1, r2 - r1 + 1)
end, { buffer = true, desc = "Set header to the selected range of lines" })


-- FIXME: bug with header row where comment rows before it are excluded.
-- when scrolling to side it jumps up one line.
