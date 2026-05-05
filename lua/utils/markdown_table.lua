--- Markdown table alignment.
---
--- Pads cells with spaces so that columns line up when read at conceallevel 1+
--- (which conceals emphasis/code/link syntax). Width is measured with
--- vim.api.nvim_strwidth after stripping common concealed delimiters, so the
--- aligned text in the buffer is deterministic — re-aligning is a no-op and
--- diffs stay clean. With conceallevel=0 the visual width of cells containing
--- emphasis or links will exceed the column; that's a deliberate trade-off
--- (one stable on-disk format vs. matching every possible runtime conceal).
---
--- @class utils.markdown_table
local M = {}

--- @class utils.markdown_table.Opts
--- @field max_width? integer Cap on per-cell measured width. Cells longer
---   than this don't expand the column; their content overflows past the
---   trailing space and the row's later cells shift right. nil = no cap.

--- Check if a line is a markdown table row (starts with optional whitespace then `|`).
--- @param line string
--- @return boolean
function M.is_table_line(line)
    return line:match("^%s*|") ~= nil
end

--- Split a string on unescaped `|` delimiters.
--- `\|` is a literal pipe, `\\` is a literal backslash.
--- @param s string
--- @return string[]
local function split_on_pipes(s)
    local parts = {}
    local cur = ""
    local i = 1
    local len = #s
    while i <= len do
        local ch = s:sub(i, i)
        if ch == "\\" and i < len then
            local next_ch = s:sub(i + 1, i + 1)
            if next_ch == "|" or next_ch == "\\" then
                cur = cur .. ch .. next_ch
                i = i + 2
            else
                cur = cur .. ch
                i = i + 1
            end
        elseif ch == "|" then
            parts[#parts + 1] = cur
            cur = ""
            i = i + 1
        else
            cur = cur .. ch
            i = i + 1
        end
    end
    parts[#parts + 1] = cur
    return parts
end

--- Parse a markdown table row into trimmed cells.
--- @param line string
--- @return string[]
local function parse_row(line)
    local cells = {}
    local inner = line:match("^%s*|(.*)$")
    if not inner then
        return cells
    end
    local parts = split_on_pipes(inner)
    -- Drop trailing empty part from the closing `|`
    if #parts > 0 and vim.trim(parts[#parts]) == "" then
        parts[#parts] = nil
    end
    for _, cell in ipairs(parts) do
        cells[#cells + 1] = vim.trim(cell)
    end
    return cells
end

--- True if every cell looks like `---`, `:---`, `---:`, or `:---:`.
--- @param cells string[]
--- @return boolean
local function is_separator_row(cells)
    for _, cell in ipairs(cells) do
        if not cell:match("^:?%-+:?$") then
            return false
        end
    end
    return #cells > 0
end

--- Build a separator cell of given width preserving alignment markers.
--- @param original string
--- @param width integer
--- @return string
local function build_separator_cell(original, width)
    local left = original:match("^:") and ":" or ""
    local right = original:match(":$") and ":" or ""
    local dashes = width - #left - #right
    if dashes < 1 then
        dashes = 1
    end
    return left .. string.rep("-", dashes) .. right
end

--- Visual width of a cell after stripping markdown syntax that's concealed
--- at conceallevel >= 1. Tracks nvim-treesitter's bundled markdown_inline
--- highlights query plus the local override in
--- `after/queries/markdown_inline/highlights.scm` — which conceals the four
--- ~ delimiters of true `~~double~~` strikethroughs while leaving stray
--- single tildes visible. Autolink brackets (`<https://...>`) are NOT
--- concealed and intentionally not stripped here. Order matters: code spans
--- first (their content is opaque), then images before links (`![..](..)`
--- would otherwise match `[..](..)`), then longest emphasis delimiters
--- before shorter.
--- @param cell string
--- @return integer
local function cell_visual_width(cell)
    local s = cell
    s = s:gsub("`([^`]+)`", "%1")
    s = s:gsub("!%[([^%]]*)%]%([^)]*%)", "%1") -- image
    s = s:gsub("%[([^%]]*)%]%([^)]*%)", "%1") -- inline link
    s = s:gsub("%[([^%]]*)%]%[[^%]]*%]", "%1") -- full reference link
    s = s:gsub("%*%*%*(.-)%*%*%*", "%1")
    s = s:gsub("___(.-)___", "%1")
    s = s:gsub("%*%*(.-)%*%*", "%1")
    s = s:gsub("__(.-)__", "%1")
    s = s:gsub("%*(.-)%*", "%1")
    s = s:gsub("_(.-)_", "%1")
    s = s:gsub("~~(.-)~~", "%1") -- only ~~double~~ — single tildes stay visible
    return vim.api.nvim_strwidth(s)
end

--- Format a contiguous block of table lines (no validation that they all are).
--- @param table_lines string[]
--- @param opts? utils.markdown_table.Opts
--- @return string[]
function M.format_block(table_lines, opts)
    opts = opts or {}
    local max_width = opts.max_width

    local rows = {}
    for _, line in ipairs(table_lines) do
        rows[#rows + 1] = parse_row(line)
    end

    local num_cols = 0
    for _, row in ipairs(rows) do
        if #row > num_cols then
            num_cols = #row
        end
    end
    if num_cols == 0 then
        return table_lines
    end

    local sep_idx = nil
    for i, row in ipairs(rows) do
        if is_separator_row(row) then
            sep_idx = i
            break
        end
    end

    local col_widths = {}
    for c = 1, num_cols do
        col_widths[c] = 0
    end
    for idx, row in ipairs(rows) do
        if idx ~= sep_idx then
            for c = 1, num_cols do
                local cell = row[c] or ""
                local vw = cell_visual_width(cell)
                if max_width and vw > max_width then
                    vw = max_width
                end
                if vw > col_widths[c] then
                    col_widths[c] = vw
                end
            end
        end
    end

    -- Separator dashes need at least 3 to remain valid markdown.
    for c = 1, num_cols do
        if col_widths[c] < 3 then
            col_widths[c] = 3
        end
    end

    local result = {}
    for _, row in ipairs(rows) do
        local parts = {}
        local is_sep = is_separator_row(row)
        for c = 1, num_cols do
            local cell = row[c] or ""
            if is_sep then
                parts[#parts + 1] = build_separator_cell(cell, col_widths[c])
            else
                local pad = col_widths[c] - cell_visual_width(cell)
                if pad < 0 then
                    pad = 0 -- overflow: cell wider than capped column
                end
                parts[#parts + 1] = cell .. string.rep(" ", pad)
            end
        end
        result[#result + 1] = "| " .. table.concat(parts, " | ") .. " |"
    end

    return result
end

--- Find the contiguous range of table lines covering `row` (0-indexed).
--- Returns nil if `row` is not on a table line.
--- @param bufnr integer
--- @param row integer 0-indexed
--- @return integer? start_row 0-indexed, inclusive
--- @return integer? end_row 0-indexed, inclusive
function M.find_block_range(bufnr, row)
    local line = vim.api.nvim_buf_get_lines(bufnr, row, row + 1, false)[1]
    if not line or not M.is_table_line(line) then
        return nil, nil
    end
    local total = vim.api.nvim_buf_line_count(bufnr)
    local s = row
    while s > 0 do
        local prev = vim.api.nvim_buf_get_lines(bufnr, s - 1, s, false)[1]
        if not prev or not M.is_table_line(prev) then
            break
        end
        s = s - 1
    end
    local e = row
    while e < total - 1 do
        local nxt = vim.api.nvim_buf_get_lines(bufnr, e + 1, e + 2, false)[1]
        if not nxt or not M.is_table_line(nxt) then
            break
        end
        e = e + 1
    end
    return s, e
end

--- Align the table block containing the given row, in place.
--- @param bufnr integer
--- @param row integer 0-indexed
--- @param opts? utils.markdown_table.Opts
--- @return boolean aligned True if a table block was found and replaced.
function M.align_at(bufnr, row, opts)
    local s, e = M.find_block_range(bufnr, row)
    if not s or not e then
        return false
    end
    local lines = vim.api.nvim_buf_get_lines(bufnr, s, e + 1, false)
    local formatted = M.format_block(lines, opts)
    vim.api.nvim_buf_set_lines(bufnr, s, e + 1, false, formatted)
    return true
end

return M
