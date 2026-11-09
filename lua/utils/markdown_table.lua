--- Markdown (GFM) pipe table alignment, parsed with treesitter.
---
--- Widths are of the raw text, so the aligned file lines up wherever markup
--- is shown verbatim (conceallevel 0, `cat`, raw view). Markup that conceal
--- hides at a cell edge takes the place of the space next to the pipe.
--- `layout` and `pads` give the virtual spaces that line a table up on screen
--- under conceal, without changing the file.
---
--- @class utils.markdown_table
local M = {}

local Conceal = require "agentic.utils.conceal"
local Ts = require "agentic.utils.treesitter"

--- Filetypes whose buffers hold markdown tables.
M.FILETYPES = { "markdown", "pandoc", "AgenticInput" }
--- Default cap on a cell's contribution to its column width.
M.MAX_WIDTH = 40

--- @alias utils.markdown_table.Align "left"|"centre"|"right"

--- @class utils.markdown_table.Cell
--- @field text string raw cell text, trimmed
--- @field col integer 0-indexed start byte column of the trimmed text
--- @field end_col integer 0-indexed end byte column, exclusive
--- @field from integer 0-indexed byte column after the `|` before the cell, or of the row start
--- @field to integer 0-indexed byte column of the `|` after the cell, or of the row end
--- @field piped boolean whether a `|` is on both sides of the cell

--- @class utils.markdown_table.Row
--- @field row integer 0-indexed buffer row
--- @field line string the row's buffer line
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
--- @param piped boolean whether a `|` is on both sides
--- @return utils.markdown_table.Cell
local function trimmed_cell(line, from, to, piped)
    local raw = line:sub(from + 1, to)
    local lead = #raw:match("^%s*")
    local text = vim.trim(raw)
    --- @type utils.markdown_table.Cell
    local cell = { text = text, col = from + lead, end_col = from + lead + #text, from = from, to = to, piped = piped }
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
        local piped = (i > 1 or leading) and (i < #bounds - 1 or trailing)
        cells[#cells + 1] = trimmed_cell(line, bounds[i] + 1, bounds[i + 1], piped)
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
                line = line,
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

--- Tab stop for measuring rows. Tabs inside a table are not supported: pads
--- move the tab stops after them.
local TABSTOP = 8

--- Whether conceal hides `byte`: of the ranges covering it, the one with the
--- highest priority (the last on a tie) has an empty replacement.
--- @param ranges agentic.utils.Conceal.Range[] in insertion order
--- @param byte integer 0-indexed
--- @return boolean
local function hidden(ranges, byte)
    local decider --- @type agentic.utils.Conceal.Range|nil
    for _, range in ipairs(ranges) do
        if range.col <= byte and byte < range.end_col and (not decider or range.priority >= decider.priority) then
            decider = range
        end
    end
    return decider ~= nil and decider.text == ""
end

--- Bytes from the first to after the last character of a cell's text that
--- conceal does not hide.
--- @param cell utils.markdown_table.Cell
--- @param ranges agentic.utils.Conceal.Range[] conceal ranges on the cell's row
--- @return integer vis_col 0-indexed; `cell.end_col` when every character is hidden
--- @return integer vis_end 0-indexed, exclusive; `vis_col` when every character is hidden
local function visible_span(cell, ranges)
    local vis_col, vis_end = cell.col, cell.end_col
    while vis_col < cell.end_col and hidden(ranges, vis_col) do
        vis_col = vis_col + 1
    end
    while vis_end > vis_col and hidden(ranges, vis_end - 1) do
        vis_end = vis_end - 1
    end
    return vis_col, vis_end
end

--- Display widths of the parts of a cell, from its `from` to its `to`.
--- @class utils.markdown_table.CellParts
--- @field ws_l integer whitespace before the text
--- @field mark_l integer hidden markup from the text start to the first visible character
--- @field body integer first to last visible character; 0 without one
--- @field mark_r integer hidden markup from the last visible character to the text end
--- @field ws_r integer whitespace after the text

--- Parts of a cell, measured in `cols`.
--- @param cols table<integer, integer> vcol at each character start byte of the row
--- @param cell utils.markdown_table.Cell
--- @param vis_col integer 0-indexed byte of the first visible character
--- @param vis_end integer 0-indexed byte after the last visible character
--- @return utils.markdown_table.CellParts
local function cell_parts(cols, cell, vis_col, vis_end)
    --- @type utils.markdown_table.CellParts
    local parts = {
        ws_l = cols[cell.col] - cols[cell.from],
        mark_l = cols[vis_col] - cols[cell.col],
        body = cols[vis_end] - cols[vis_col],
        mark_r = cols[cell.end_col] - cols[vis_end],
        ws_r = cols[cell.to] - cols[cell.end_col],
    }
    return parts
end

--- Display widths of a cell around its body.
--- @class utils.markdown_table.Sides
--- @field pre integer before the body
--- @field body integer 0 without a body
--- @field post integer after the body

--- Sides of a cell as it is: whitespace and markup count.
--- @param parts utils.markdown_table.CellParts
--- @return utils.markdown_table.Sides
local function layout_sides(parts)
    --- @type utils.markdown_table.Sides
    local sides = { pre = parts.ws_l + parts.mark_l, body = parts.body, post = parts.mark_r + parts.ws_r }
    return sides
end

--- Sides of a cell as formatted: markup at an edge is the gap next to the
--- pipe, else the gap is 1 where the cell has whitespace there, else 0.
--- @param parts utils.markdown_table.CellParts
--- @return utils.markdown_table.Sides
local function format_sides(parts)
    --- @type utils.markdown_table.Sides
    local sides = {
        pre = parts.mark_l > 0 and parts.mark_l or math.min(parts.ws_l, 1),
        body = parts.body,
        post = parts.mark_r > 0 and parts.mark_r or math.min(parts.ws_r, 1),
    }
    return sides
end

--- Where a column places its bodies.
--- @class utils.markdown_table.Column
--- @field align utils.markdown_table.Align
--- @field width integer display cells between the column's pipes
--- @field lead integer widest `pre` of a cell with a body
--- @field trail integer widest `post` of a cell with a body

--- Narrowest column that fits every cell of `sides` with its body placed per
--- `align`: left bodies start at `lead`, right bodies end at `width - trail`,
--- centred bodies are centred. A cell without a body needs `pre + post`.
--- @param sides utils.markdown_table.Sides[]
--- @param align utils.markdown_table.Align
--- @param min_width integer
--- @return utils.markdown_table.Column
local function column_geometry(sides, align, min_width)
    --- @type utils.markdown_table.Column
    local column = { align = align, width = min_width, lead = 0, trail = 0 }
    for _, s in ipairs(sides) do
        if s.body > 0 then
            column.lead = math.max(column.lead, s.pre)
            column.trail = math.max(column.trail, s.post)
        end
    end
    for _, s in ipairs(sides) do
        local width = s.body == 0 and s.pre + s.post
            or align == "left" and column.lead + s.body + s.post
            or align == "right" and s.pre + s.body + column.trail
            or math.max(s.body + 2 * s.pre, s.body + 2 * s.post - 1)
        column.width = math.max(column.width, width)
    end
    return column
end

--- Display cells from the start of `column` to a cell's body: the column's
--- place for it, or after its `pre` when that is further. For a cell
--- without a body, its `pre`.
--- @param column utils.markdown_table.Column
--- @param sides utils.markdown_table.Sides
--- @return integer
local function body_start(column, sides)
    if sides.body == 0 then
        return sides.pre
    end
    local anchor = column.align == "left" and column.lead
        or column.align == "right" and column.width - column.trail - sides.body
        or math.floor((column.width - sides.body) / 2)
    return math.max(sides.pre, anchor)
end

--- Spaces to add before and after a cell in its slot of `column`, for a
--- cell that starts `excess` cells past the column start. The spaces put
--- its body at `body_start` and end the slot at the column end. Where the
--- slot is too narrow, the spaces after the cell go first, then those before.
--- @param column utils.markdown_table.Column
--- @param sides utils.markdown_table.Sides that place the body
--- @param before integer display cells the cell has before its body
--- @param width integer display cells of the cell
--- @param excess integer display cells
--- @return integer left
--- @return integer right
local function slot_spaces(column, sides, before, width, excess)
    local room = column.width - excess
    local left = math.max(0, math.min(body_start(column, sides) - excess - before, room - width))
    return left, math.max(0, room - left - width)
end

--- Conceal ranges of the cells of a row, each cell's text parsed alone as
--- markdown_inline.
--- @param row utils.markdown_table.Row
--- @return agentic.utils.Conceal.Range[] byte columns on `row.line`
local function cell_conceals(row)
    local ranges = {}
    for _, cell in ipairs(row.cells) do
        local trees = Ts.string_highlight_trees(cell.text, "markdown_inline")
        for _, range in ipairs(Ts.conceal_ranges(cell.text, trees, 0, 1)) do
            ranges[#ranges + 1] = vim.tbl_extend("force", range, {
                col = cell.col + range.col,
                end_col = cell.col + math.min(range.end_col, #cell.text),
            })
        end
    end
    return ranges
end

--- Lines of an aligned table, one per row. Each column is just wide enough
--- for every body to keep its gap to each pipe: markup at an edge, else 1
--- space where the input has whitespace there, else none. The delimiter row
--- keeps its gaps and gets at least 3 dashes. It has no other effect on
--- column widths. A cell wider than `max_width` makes its column fit
--- `max_width` within the widest gaps of the column's other cells. Where it
--- does not fit its column, it loses its whitespace, then overruns.
--- The rest of its row loses whitespace until it is back on the column grid.
--- @param tbl utils.markdown_table.Table
--- @param opts utils.markdown_table.Opts|nil
--- @return string[]
function M.format(tbl, opts)
    opts = opts or {}
    local parts = {} -- per row, per column
    local sides, long_parts, min_widths = {}, {}, {} -- per column
    for c = 1, tbl.ncols do
        sides[c], long_parts[c], min_widths[c] = {}, {}, 3
    end
    for i, row in ipairs(tbl.rows) do
        local conceals = row.delimiter and {} or cell_conceals(row)
        local cols = Conceal.columns(row.line, {}, {}, TABSTOP, 0).raw
        parts[i] = {}
        for c, cell in ipairs(row.cells) do
            local p = cell_parts(cols, cell, visible_span(cell, conceals))
            local s = format_sides(p)
            local text_width = p.mark_l + p.body + p.mark_r
            parts[i][c] = p
            if row.delimiter then
                min_widths[c] = math.max(min_widths[c], 3 + s.pre + s.post)
            elseif text_width <= (opts.max_width or text_width) then
                table.insert(sides[c], s)
            else
                table.insert(long_parts[c], p)
            end
        end
    end
    local columns = {}
    for c = 1, tbl.ncols do
        local column = column_geometry(sides[c], opts.align_cells and tbl.align[c] or "left", min_widths[c])
        -- The gaps of the other cells, not the long cell's own whitespace,
        -- which formatting removes.
        for _, p in ipairs(long_parts[c]) do
            local width = math.max(column.lead, p.mark_l) + opts.max_width + math.max(column.trail, p.mark_r)
            column.width = math.max(column.width, width)
        end
        columns[c] = column
    end
    local lines = {}
    for i, row in ipairs(tbl.rows) do
        local slots = {}
        local excess = 0 -- how far the row so far overruns the column grid
        for c, column in ipairs(columns) do
            local cell, p = row.cells[c], parts[i][c]
            if row.delimiter then
                local s = p and format_sides(p) or { pre = 0, post = 0 }
                slots[c] = string.rep(" ", s.pre) .. delimiter_cell(cell and cell.text or "", column.width - s.pre - s.post)
                    .. string.rep(" ", s.post)
            elseif cell then
                local text_width = p.mark_l + p.body + p.mark_r
                local left, right = slot_spaces(column, format_sides(p), p.mark_l, text_width, excess)
                slots[c] = string.rep(" ", left) .. cell.text .. string.rep(" ", right)
            else
                slots[c] = string.rep(" ", math.max(0, column.width - excess))
            end
            excess = math.max(0, excess + vim.api.nvim_strwidth(slots[c]) - column.width)
        end
        lines[#lines + 1] = row.prefix .. "|" .. table.concat(slots, "|") .. "|"
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

--- Conceal ranges on each row of a table, from the active treesitter
--- highlighter and from extmarks, in drawing order. Without a highlighter
--- on `buf`, only the extmark ranges, as on screen.
--- @param buf integer
--- @param tbl utils.markdown_table.Table
--- @return table<integer, agentic.utils.Conceal.Range[]> conceals per index into `tbl.rows`
function M.conceals(buf, tbl)
    local srow, erow = tbl.rows[1].row, tbl.rows[#tbl.rows].row + 1
    -- The highlighter parses injections of drawn rows only.
    vim.treesitter.get_parser(buf):parse({ srow, 0, erow, 0 })
    local ranges = Ts.conceal_ranges(buf, Ts.highlight_trees(buf, srow, erow), srow, erow)
    vim.list_extend(ranges, Conceal.extmark_ranges(buf, srow, erow))
    local index, conceals = {}, {}
    for i, row in ipairs(tbl.rows) do
        index[row.row] = i
        conceals[i] = {}
    end
    for _, range in ipairs(ranges) do
        local i = index[range.row]
        if i then
            table.insert(conceals[i], range)
        end
    end
    return conceals
end

--- @class utils.markdown_table.CellLayout
--- @field vis_col integer 0-indexed byte of the first character conceal does not hide; `end_col` when there is none
--- @field vis_end integer 0-indexed byte after the last character conceal does not hide; `vis_col` when there is none
--- @field slotted boolean false when the cell's raw width is over the cap

--- On-screen layout of a table, from its raw text. Markup is text conceal
--- hides.
--- @class utils.markdown_table.Layout
--- @field cells table<integer, utils.markdown_table.CellLayout>[] per row, per column; no entry for a cell without both pipes
--- @field columns utils.markdown_table.Column[] at conceallevel 0

--- Display layout of a table. Each column is the narrowest that fits its
--- cells as they are at conceallevel 0, whitespace and markup included,
--- with every body placed per the column's alignment. The delimiter row
--- counts only by its width. Cells over the cap do not count.
--- @param tbl utils.markdown_table.Table
--- @param conceals table<integer, agentic.utils.Conceal.Range[]> per index into `tbl.rows`
--- @param opts utils.markdown_table.Opts|nil
--- @return utils.markdown_table.Layout
function M.layout(tbl, conceals, opts)
    opts = opts or {}
    --- @type utils.markdown_table.Layout
    local layout = { cells = {}, columns = {} }
    local sides, min_widths = {}, {} -- per column
    for c = 1, tbl.ncols do
        sides[c], min_widths[c] = {}, 0
    end
    for i, row in ipairs(tbl.rows) do
        local raw = Conceal.columns(row.line, conceals[i], {}, TABSTOP, 0).raw
        layout.cells[i] = {}
        for c, cell in ipairs(row.cells) do
            if cell.piped then
                local vis_col, vis_end = visible_span(cell, conceals[i])
                local s = layout_sides(cell_parts(raw, cell, vis_col, vis_end))
                local text_width = raw[cell.end_col] - raw[cell.col]
                --- @type utils.markdown_table.CellLayout
                local cl = {
                    vis_col = vis_col,
                    vis_end = vis_end,
                    slotted = text_width <= (opts.max_width or text_width),
                }
                layout.cells[i][c] = cl
                if row.delimiter then
                    min_widths[c] = math.max(min_widths[c], s.pre + s.body + s.post)
                elseif cl.slotted then
                    table.insert(sides[c], s)
                end
            end
        end
    end
    for c = 1, tbl.ncols do
        layout.columns[c] = column_geometry(sides[c], opts.align_cells and tbl.align[c] or "left", min_widths[c])
    end
    return layout
end

--- Virtual spaces that put every cell of a table in its column on screen,
--- for a window at `level`: each body at its place in the column, measured
--- at `level`, and each cell to the column end. A cell wider than its column
--- shifts the rest of its row right. Later cells in that row get less space
--- until the row is back on the column grid.
--- @param tbl utils.markdown_table.Table
--- @param layout utils.markdown_table.Layout
--- @param conceals table<integer, agentic.utils.Conceal.Range[]> per index into `tbl.rows`
--- @param level integer 'conceallevel'
--- @param unconcealed_row integer|nil 0-indexed buffer row drawn without conceal
--- @return { row: integer, col: integer, width: integer }[] pads `width` spaces drawn before 0-indexed `row`, `col`
function M.pads(tbl, layout, conceals, level, unconcealed_row)
    local pads = {}
    local function add(row, col, width)
        if width > 0 then
            pads[#pads + 1] = { row = row, col = col, width = width }
        end
    end
    for i, row in ipairs(tbl.rows) do
        local row_level = row.row == unconcealed_row and 0 or level
        local shown = Conceal.columns(row.line, conceals[i], {}, TABSTOP, row_level).shown
        local excess = 0 -- how far the row so far overruns the column grid
        for c = 1, tbl.ncols do
            local cl = layout.cells[i][c]
            if cl then
                local cell, column = row.cells[c], layout.columns[c]
                local s = layout_sides(cell_parts(shown, cell, cl.vis_col, cl.vis_end))
                local total = s.pre + s.body + s.post
                local left, right = 0, column.width - excess - total
                if not row.delimiter then
                    left, right = slot_spaces(column, s, s.pre, total, excess)
                end
                add(row.row, cell.col, left)
                add(row.row, cell.end_col, right)
                excess = math.max(0, excess + left + total - column.width)
            end
        end
    end
    return pads
end

return M
