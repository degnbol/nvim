-- Text objects for mini.ai by id: a `custom_textobjects` spec and a
-- description. An entry without a spec describes a mini.ai built-in.
-- `spec = false` disables mini.ai's object, so Vim's own applies.
local ai = require "mini.ai"
local util = require "utils/init"
local ts = require "utils/treesitter"

---A mini.ai region in the current buffer: 1-indexed lines and byte columns,
---`to` inclusive. An empty region has no `to` (`:h MiniAi-glossary`).
---@alias TextobjectRegion { from: { line: integer, col: integer }, to: { line: integer, col: integer }|nil }

---Region of a range in the current buffer, as `MiniAi.gen_spec.treesitter`
---converts a tree-sitter range. A range that ends at column 0 ends on the
---newline of the line before. An empty range gives an empty region.
---@param range Range4 0-indexed, end exclusive
---@return TextobjectRegion
local function range_region(range)
    local from = { line = range[1] + 1, col = range[2] + 1 }
    if range[1] == range[3] and range[2] == range[4] then return { from = from } end
    local to = { line = range[3] + 1, col = range[4] }
    if to.col == 0 then
        to.line = to.line - 1
        to.col = #util.get_line(to.line - 1) + 1
    end
    return { from = from, to = to }
end

---Range of a region of the current buffer, the inverse of `range_region`. The
---range ends after the last byte of the character at `to`. A region that ends
---on a newline gives a range that ends at column 0 of the next row.
---@param region TextobjectRegion
---@return Range4 0-indexed, end exclusive; empty at `from` for an empty region
local function region_range(region)
    local from, to = region.from, region.to
    if not to then return { from.line - 1, from.col - 1, from.line - 1, from.col - 1 } end
    local line = util.get_line(to.line - 1)
    if to.col > #line then return { from.line - 1, from.col - 1, to.line, 0 } end
    return { from.line - 1, from.col - 1, to.line - 1, to.col + vim.str_utf_end(line, to.col) }
end

---Whether the `textobjects` query of a language tree has a capture.
---@param lang_tree vim.treesitter.LanguageTree
---@param capture string without `@`, e.g. "comment.outer"
---@return boolean
local function has_capture(lang_tree, capture)
    local query = vim.treesitter.query.get(lang_tree:lang(), "textobjects")
    return query ~= nil and vim.list_contains(query.captures, capture)
end

---Rows of mini.ai's search: those of the reference region and `opts.n_lines`
---on each side.
---@param opts { reference_region: TextobjectRegion, n_lines: integer } see `MiniAi.find_textobject`
---@return integer start first row, 0-indexed
---@return integer stop row after the last, 0-indexed
local function search_rows(opts)
    local reference = opts.reference_region
    local to = reference.to or reference.from
    return reference.from.line - 1 - opts.n_lines, to.line + opts.n_lines
end

---One region per full-line comment, from `leader` to the end of the line.
---@param buf integer
---@param leader string what a comment opens with, e.g. "#"
---@return TextobjectRegion[]
local function comment_line_regions(buf, leader)
    local regions = {}
    for row, line in ipairs(vim.api.nvim_buf_get_lines(buf, 0, -1, false)) do
        if util.is_comment_line(line, leader) then
            regions[#regions + 1] = {
                from = { line = row, col = line:find "%S" },
                to = { line = row, col = #line },
            }
        end
    end
    return regions
end

---Position of byte `offset` of `text`, where `text` starts at `from`.
---@param from { line: integer, col: integer }
---@param text string lines joined by "\n"
---@param offset integer 1-indexed
---@return { line: integer, col: integer }
local function offset_position(from, text, offset)
    local before = text:sub(1, offset - 1)
    local _, newlines = before:gsub("\n", "")
    if newlines == 0 then return { line = from.line, col = from.col + offset - 1 } end
    return { line = from.line + newlines, col = offset - before:match ".*\n()" + 1 }
end

---Text of a `region` of the current buffer, lines joined by "\n".
---A newline at its end is left out.
---@param region TextobjectRegion
---@return string "" for an empty region
local function region_text(region)
    if not region.to then return "" end
    local to_line = util.get_line(region.to.line - 1)
    return table.concat(vim.api.nvim_buf_get_text(0,
        region.from.line - 1, region.from.col - 1,
        region.to.line - 1, math.min(region.to.col, #to_line), {}), "\n")
end

---The region of bytes `first` to `last` of `text`.
---@param region TextobjectRegion
---@param text string text of `region`, lines joined by "\n"
---@param first integer 1-indexed
---@param last integer 1-indexed, inclusive
---@return TextobjectRegion empty at `first` when `first > last`
local function subregion(region, text, first, last)
    local from = offset_position(region.from, text, first)
    if first > last then return { from = from } end
    return { from = from, to = offset_position(region.from, text, last) }
end

---`region` of the current buffer without its comment delimiters. The delimiter
---pair used is the one with the longest left delimiter that `region` starts
---with, e.g. `---` over `--`. Its right delimiter is removed only when
---`region` ends with it.
---@param region TextobjectRegion
---@param delimiters { left: string, right: string|nil }[] e.g. `{ left = "/*", right = "*/" }`
---@return TextobjectRegion
local function strip_delimiters(region, delimiters)
    local text = region_text(region)
    local delimiter = { left = "" }
    for _, d in ipairs(delimiters) do
        if #d.left > #delimiter.left and vim.startswith(text, d.left) then delimiter = d end
    end
    local last = #text
    local right = delimiter.right or ""
    if right ~= "" and vim.endswith(text, right) then last = last - #right end
    return subregion(region, text, #delimiter.left + 1, last)
end

---`region` of the current buffer without its leading and trailing whitespace.
---@param region TextobjectRegion
---@return TextobjectRegion empty at its end when it is all whitespace
local function trim_whitespace(region)
    if not region.to then return region end
    local text = region_text(region)
    return subregion(region, text, text:match "^%s*()", text:match "()%s*$" - 1)
end

---Whether position `p` is strictly before position `q`.
---@param p { line: integer, col: integer }
---@param q { line: integer, col: integer }
---@return boolean
local function before(p, q) return util.before(p.line, p.col, q.line, q.col) end

---The part of region `a` that is also in region `b`.
---@param a TextobjectRegion
---@param b TextobjectRegion
---@return TextobjectRegion empty at the later start when they do not overlap
local function intersect(a, b)
    local from = before(a.from, b.from) and b.from or a.from
    local to = a.to and b.to and (before(a.to, b.to) and a.to or b.to)
    if not to or before(to, from) then return { from = from } end
    return { from = from, to = to }
end

---Comment regions in the current buffer, as a mini.ai spec. `a` is
---`@comment.outer` of the language at the cursor, or each full-line comment
---when that language's `textobjects` query has no such capture. `i` is the
---part of `a` inside both the delimiters of the 'comments' and
---'commentstring' at the cursor and the match's `@comment.inner`, if any,
---without leading and trailing whitespace. When `a` of a comment covers the
---reference region, `i` gives only that comment, since mini.ai checks
---coverage against the `i` candidates.
---@param ai_type "a"|"i"
---@param opts { reference_region: TextobjectRegion, n_lines: integer } see `MiniAi.find_textobject`
---@return TextobjectRegion|TextobjectRegion[] candidates
local function comment(ai_type, _, opts)
    local buf = vim.api.nvim_get_current_buf()
    local cursor = vim.api.nvim_win_get_cursor(0)
    local row, col = cursor[1] - 1, cursor[2]
    local reference = opts.reference_region
    local reference_to = reference.to or reference.from
    local commentstring_delimiter = util.commentstring_delimiter(ts.commentstring_at(buf, { row, col }))
    local lang_tree = ts.lang_tree_at(buf, { row, col })
    ---@type { outer: TextobjectRegion, inner: TextobjectRegion|nil }[]
    local comments = {}
    if lang_tree and has_capture(lang_tree, "comment.outer") then
        local captures = { "comment.outer", "comment.inner" }
        local start, stop = search_rows(opts)
        for _, match in ipairs(ts.capture_matches(lang_tree, captures, start, stop)) do
            if match["comment.outer"] then
                local inner = match["comment.inner"]
                comments[#comments + 1] = {
                    outer = range_region(match["comment.outer"]),
                    inner = inner and range_region(inner),
                }
            end
        end
    else
        for _, region in ipairs(comment_line_regions(buf, commentstring_delimiter.left)) do
            comments[#comments + 1] = { outer = region }
        end
    end
    if ai_type == "a" then return vim.tbl_map(function(c) return c.outer end, comments) end
    local delimiters = util.comment_delimiters(ts.option_at(buf, { row, col }, "comments"))
    delimiters[#delimiters + 1] = commentstring_delimiter
    local function inner(c)
        local region = strip_delimiters(c.outer, delimiters)
        if c.inner then region = intersect(region, c.inner) end
        return trim_whitespace(region)
    end
    for _, c in ipairs(comments) do
        if not before(reference.from, c.outer.from) and not before(c.outer.to, reference_to) then
            return inner(c)
        end
    end
    return vim.tbl_map(inner, comments)
end

---Whether a range is empty.
---@param range Range4
---@return boolean
local function is_empty(range) return range[1] == range[3] and range[2] == range[4] end

---A key/value pair of the current buffer.
---@alias KeyValue { key: Range4, value: Range4 } 0-indexed, end exclusive

---Range of the key or the value of a pair. For `ai_type` `i`, the part
---alone. For `a`, the part with the separator and its whitespace: the key
---up to the start of the value, or the value from the end of the key.
---@param pair KeyValue
---@param ai_type "a"|"i"
---@param part "key"|"value"
---@return Range4 0-indexed, end exclusive
local function part_range(pair, ai_type, part)
    local key, value = pair.key, pair.value
    if ai_type == "i" then return pair[part] end
    if part == "key" then return { key[1], key[2], value[1], value[2] } end
    return { key[3], key[4], value[3], value[4] }
end

---Key/value pairs of the `textobjects` query matches of a language tree that
---intersect a range of rows and have both `@assignment.lhs` and
---`@assignment.rhs`. Keys and values are without surrounding whitespace; a
---match whose key or value is all whitespace is left out.
---@param lang_tree vim.treesitter.LanguageTree of the current buffer
---@param start integer first row, 0-indexed
---@param stop integer row after the last, 0-indexed
---@return KeyValue[]
local function query_pairs(lang_tree, start, stop)
    local captures = { "assignment.lhs", "assignment.rhs" }
    local kv_pairs = {}
    for _, match in ipairs(ts.capture_matches(lang_tree, captures, start, stop)) do
        local lhs, rhs = match["assignment.lhs"], match["assignment.rhs"]
        if lhs and rhs then
            local key = region_range(trim_whitespace(range_region(lhs)))
            local value = region_range(trim_whitespace(range_region(rhs)))
            if not is_empty(key) and not is_empty(value) then
                kv_pairs[#kv_pairs + 1] = { key = key, value = value }
            end
        end
    end
    return kv_pairs
end

---Key/value pair in the argument (textobject `i,`) around the cursor, else
---in the cursor line, split by `util.split_pair`.
---@return KeyValue|nil
local function text_pair()
    local cursor = vim.api.nvim_win_get_cursor(0)
    local line = util.get_line(cursor[1] - 1)
    local argument = ai.find_textobject("i", ",", {
        search_method = "cover",
        reference_region = { from = { line = cursor[1], col = cursor[2] + 1 } },
    })
    local regions = { { from = { line = cursor[1], col = 1 }, to = { line = cursor[1], col = math.max(#line, 1) } } }
    if argument then table.insert(regions, 1, argument) end
    for _, region in ipairs(regions) do
        local text = region_text(region)
        local split = util.split_pair(text)
        if split then
            return {
                key = region_range(subregion(region, text, split.key[1], split.key[2])),
                value = region_range(subregion(region, text, split.value[1], split.value[2])),
            }
        end
    end
end

---Key or value regions in the current buffer, as a mini.ai spec. The pairs
---are the `textobjects` query matches of the language at the cursor with
---both `@assignment.lhs` and `@assignment.rhs`. Of the pairs whose part is
---not inside the reference region, the `opts.n_times`-th smallest whose span
---from key to value contains the reference gives the region. When none
---does, every pair's part is a candidate for mini.ai's search. A language
---whose query has no `@assignment.lhs` gets the pair in the argument around
---the cursor, else in the cursor line, and ignores `opts.n_times`.
---@param part "key"|"value"
---@return fun(ai_type: "a"|"i", id: string, opts: { reference_region: TextobjectRegion, n_times: integer, n_lines: integer }): TextobjectRegion|TextobjectRegion[]|nil spec
local function pair(part)
    return function(ai_type, _, opts)
        local cursor = vim.api.nvim_win_get_cursor(0)
        local lang_tree = ts.lang_tree_at(0, { cursor[1] - 1, cursor[2] })
        if not (lang_tree and has_capture(lang_tree, "assignment.lhs")) then
            local kv = text_pair()
            return kv and range_region(part_range(kv, ai_type, part))
        end
        local reference = region_range(opts.reference_region)
        local kv_pairs = query_pairs(lang_tree, search_rows(opts))
        local function part_region(kv) return range_region(part_range(kv, ai_type, part)) end
        local stepped_out = vim.tbl_filter(function(kv)
            return not ts.contains(reference, part_range(kv, ai_type, part))
        end, kv_pairs)
        local spans = vim.tbl_map(function(kv)
            return { kv.key[1], kv.key[2], kv.value[3], kv.value[4] }
        end, stepped_out)
        local chosen = ts.nth_enclosing(spans, reference, opts.n_times)
        if chosen then return part_region(stepped_out[chosen]) end
        return vim.tbl_map(part_region, kv_pairs)
    end
end

local treesitter = ai.gen_spec.treesitter

---@type table<string, { spec: any, desc: string|nil }>
return {
    f = { spec = treesitter { a = "@function.outer", i = "@function.inner" }, desc = "Function" },
    C = { spec = treesitter { a = "@class.outer", i = "@class.inner" }, desc = "Class" },
    a = { spec = treesitter { a = "@parameter.outer", i = "@parameter.inner" }, desc = "Parameter" },
    ["?"] = { spec = treesitter { a = "@conditional.outer", i = "@conditional.inner" }, desc = "Conditional" },
    o = { spec = treesitter { a = "@loop.outer", i = "@loop.inner" }, desc = "Loop" },
    F = { spec = treesitter { a = "@fence.outer", i = "@fence.inner" }, desc = "Code fence" },
    k = { spec = pair "key", desc = "Key" },
    v = { spec = pair "value", desc = "Value" },
    c = { spec = comment, desc = "Comment" },
    [","] = { spec = ai.gen_spec.argument(), desc = "Argument" },
    b = { desc = "Bracket" },
    q = { desc = "Quote" },
    -- mini.ai's `t` pairs an opening tag with the first closing tag of the same
    -- name, so it fails on the outer of nested same-name tags.
    t = { spec = false },
}
