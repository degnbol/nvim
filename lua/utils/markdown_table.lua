--- Markdown (GFM) pipe table alignment, parsed with treesitter.
---
--- A cell's width is the `nvim_strwidth` of its raw text, so the aligned file
--- lines up wherever markup is shown verbatim (conceallevel 0, `cat`, raw view).
---
--- @class utils.markdown_table
local M = {}

--- @alias utils.markdown_table.Align "left"|"centre"|"right"

--- @class utils.markdown_table.Cell
--- @field text string raw cell text, trimmed
--- @field col integer 0-indexed start byte column of the trimmed text
--- @field end_col integer 0-indexed end byte column, exclusive
--- @field lpipe? integer 0-indexed byte column of the `|` before the cell; nil when the row has no leading pipe
--- @field rpipe? integer 0-indexed byte column of the `|` after the cell; nil when the row has no trailing pipe

--- @class utils.markdown_table.Row
--- @field row integer 0-indexed buffer row
--- @field prefix string text before the table on this line (indent, `> ` of a blockquote, list continuation)
--- @field cells utils.markdown_table.Cell[]
--- @field delimiter boolean

--- @class utils.markdown_table.Table
--- @field node TSNode the `pipe_table` node
--- @field rows utils.markdown_table.Row[] header first, then delimiter, then body
--- @field align utils.markdown_table.Align[] per column, from the delimiter row
--- @field ncols integer max cell count over rows

--- @class utils.markdown_table.Opts
--- @field max_width? integer cap on a cell's contribution to its column width; nil = no cap
--- @field align_cells? boolean pad per alignment colons (right/centre); default false = all left

local ROW_TYPES = {
    pipe_table_header = true,
    pipe_table_delimiter_row = true,
    pipe_table_row = true,
}

--- `pipe_table` nodes that intersect rows [start_row, end_row], including
--- tables in injected markdown.
--- @param buf integer
--- @param start_row integer 0-indexed, inclusive
--- @param end_row integer 0-indexed, inclusive
--- @return TSNode[] nodes sorted by start row; empty when the buffer has no parser
function M.tables(buf, start_row, end_row)
    local parser = vim.treesitter.get_parser(buf)
    if not parser then
        return {}
    end
    parser:parse({ start_row, 0, end_row + 1, 0 })
    local query = vim.treesitter.query.parse("markdown", "(pipe_table) @table")
    local nodes = {}
    parser:for_each_tree(function(tree, ltree)
        if ltree:lang() ~= "markdown" then
            return
        end
        for _, node in query:iter_captures(tree:root(), buf, start_row, end_row + 1) do
            local sr, _, er, ec = node:range()
            -- The range ends at column 0 of the row below the table.
            local last_row = ec == 0 and er - 1 or er
            if sr <= end_row and last_row >= start_row then
                nodes[#nodes + 1] = node
            end
        end
    end)
    table.sort(nodes, function(a, b) return a:start() < b:start() end)
    return nodes
end

--- Byte columns of the `|` tokens of a table row, in order.
--- @param row_node TSNode
--- @return integer[] pipes 0-indexed byte columns
local function row_pipes(row_node)
    local pipes = {}
    local function add(node)
        if node:type() == "|" then
            pipes[#pipes + 1] = select(2, node:start())
        end
    end
    for child in row_node:iter_children() do
        -- The parser emits `ERROR "|"` in place of a cell for an empty cell `||`.
        if child:type() == "ERROR" then
            for token in child:iter_children() do
                add(token)
            end
        else
            add(child)
        end
    end
    return pipes
end

--- Trimmed cell between byte columns [from, to) of `line`.
--- @param line string
--- @param from integer 0-indexed
--- @param to integer 0-indexed, exclusive
--- @return utils.markdown_table.Cell
local function trimmed_cell(line, from, to)
    local raw = line:sub(from + 1, to)
    local lead = #raw:match("^%s*")
    local text = vim.trim(raw)
    --- @type utils.markdown_table.Cell
    local cell = { text = text, col = from + lead, end_col = from + lead + #text }
    return cell
end

--- Cells of a table row, split at its `|` tokens.
--- @param line string the row's buffer line
--- @param row_node TSNode
--- @return utils.markdown_table.Cell[]
local function row_cells(line, row_node)
    local _, row_start, _, row_end = row_node:range()
    local bounds = row_pipes(row_node)
    local leading = bounds[1] == row_start
    local trailing = bounds[#bounds] == row_end - 1
    if not leading then
        table.insert(bounds, 1, row_start - 1)
    end
    if not trailing then
        bounds[#bounds + 1] = row_end
    end
    local cells = {}
    for i = 1, #bounds - 1 do
        local cell = trimmed_cell(line, bounds[i] + 1, bounds[i + 1])
        cell.lpipe = (i > 1 or leading) and bounds[i] or nil
        cell.rpipe = (i < #bounds - 1 or trailing) and bounds[i + 1] or nil
        cells[#cells + 1] = cell
    end
    return cells
end

--- Alignment of a delimiter cell from its colons.
--- @param text string trimmed delimiter cell, e.g. `:--:`
--- @return utils.markdown_table.Align
local function delimiter_align(text)
    local left, right = text:sub(1, 1) == ":", text:sub(-1) == ":"
    return left and right and "centre" or right and "right" or "left"
end

--- Parse a `pipe_table` node.
--- @param buf integer
--- @param node TSNode the `pipe_table` node
--- @return utils.markdown_table.Table
function M.parse(buf, node)
    local sr, _, er = node:range()
    local lines = vim.api.nvim_buf_get_lines(buf, sr, er + 1, false)
    local rows = {}
    local align = {}
    local ncols = 0
    for child in node:iter_children() do
        if ROW_TYPES[child:type()] then
            local row, col = child:start()
            local line = lines[row - sr + 1]
            local delimiter = child:type() == "pipe_table_delimiter_row"
            local cells = row_cells(line, child)
            if delimiter then
                for c, cell in ipairs(cells) do
                    align[c] = delimiter_align(cell.text)
                end
            end
            ncols = math.max(ncols, #cells)
            rows[#rows + 1] = {
                row = row,
                prefix = line:sub(1, col),
                cells = cells,
                delimiter = delimiter,
            }
        end
    end
    --- @type utils.markdown_table.Table
    local tbl = { node = node, rows = rows, align = align, ncols = ncols }
    return tbl
end

--- Delimiter cell of `width` keeping the colons of `text`.
--- @param text string original delimiter cell, empty for a missing cell
--- @param width integer
--- @return string
local function delimiter_cell(text, width)
    local left = text:sub(1, 1) == ":" and ":" or ""
    local right = text:sub(-1) == ":" and ":" or ""
    return left .. string.rep("-", width - #left - #right) .. right
end

--- Pad `text` to display `width` per `align`. Text wider than `width` is not padded.
--- @param text string
--- @param width integer display cells
--- @param align utils.markdown_table.Align
--- @return string
local function pad(text, width, align)
    local n = math.max(0, width - vim.api.nvim_strwidth(text))
    local left = align == "right" and n or align == "centre" and math.floor(n / 2) or 0
    return string.rep(" ", left) .. text .. string.rep(" ", n - left)
end

--- Lines of an aligned table, one per row. Column width is the widest
--- non-delimiter cell, each capped at `max_width`, and at least 3.
--- @param tbl utils.markdown_table.Table
--- @param opts utils.markdown_table.Opts|nil
--- @return string[]
function M.format(tbl, opts)
    opts = opts or {}
    local widths = {}
    for c = 1, tbl.ncols do
        widths[c] = 3
        for _, row in ipairs(tbl.rows) do
            local cell = row.cells[c]
            if cell and not row.delimiter then
                local w = vim.api.nvim_strwidth(cell.text)
                widths[c] = math.max(widths[c], math.min(w, opts.max_width or w))
            end
        end
    end
    local lines = {}
    for _, row in ipairs(tbl.rows) do
        local parts = {}
        for c = 1, tbl.ncols do
            local text = row.cells[c] and row.cells[c].text or ""
            parts[c] = row.delimiter and delimiter_cell(text, widths[c])
                or pad(text, widths[c], opts.align_cells and tbl.align[c] or "left")
        end
        lines[#lines + 1] = row.prefix .. "| " .. table.concat(parts, " | ") .. " |"
    end
    return lines
end

--- True when a table has a parse error between its rows, e.g. the stray `|`
--- of a row starting with `||`. Errors inside a row are empty cells.
--- @param node TSNode the `pipe_table` node
--- @return boolean
local function has_row_level_error(node)
    for child in node:iter_children() do
        if child:type() == "ERROR" or child:missing() then
            return true
        end
    end
    return false
end

--- Replace `old` with `new` on buffer row `row` by setting only the span
--- between their common prefix and suffix, so extmarks outside it stay put.
--- @param buf integer
--- @param row integer 0-indexed
--- @param old string
--- @param new string
local function set_line_diff(buf, row, old, new)
    local p = 0
    while p < #old and p < #new and old:byte(p + 1) == new:byte(p + 1) do
        p = p + 1
    end
    local s = 0
    while s < #old - p and s < #new - p and old:byte(#old - s) == new:byte(#new - s) do
        s = s + 1
    end
    vim.api.nvim_buf_set_text(buf, row, p, row, #old - s, { new:sub(p + 1, #new - s) })
end

--- Align every table that intersects rows [start_row, end_row], in place.
--- Tables with parse errors other than empty cells are skipped with a warning.
--- Only changed rows are written.
--- @param buf integer
--- @param start_row integer 0-indexed, inclusive
--- @param end_row integer 0-indexed, inclusive
--- @param opts utils.markdown_table.Opts|nil
--- @return integer n_tables tables found, including skipped ones
function M.align(buf, start_row, end_row, opts)
    local nodes = M.tables(buf, start_row, end_row)
    local tables = {}
    for _, node in ipairs(nodes) do
        if not has_row_level_error(node) then
            tables[#tables + 1] = M.parse(buf, node)
        else
            vim.notify(("markdown table at line %d has parse errors, not aligned"):format(node:start() + 1),
                vim.log.levels.WARN)
        end
    end
    -- Formatting keeps the row count, so rows of later tables stay valid.
    for _, tbl in ipairs(tables) do
        for i, new in ipairs(M.format(tbl, opts)) do
            local row = tbl.rows[i].row
            local old = vim.api.nvim_buf_get_lines(buf, row, row + 1, false)[1]
            if old ~= new then
                set_line_diff(buf, row, old, new)
            end
        end
    end
    return #nodes
end

return M
