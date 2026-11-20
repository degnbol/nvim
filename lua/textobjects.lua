-- Text objects for mini.ai by id: a `custom_textobjects` spec and a
-- description. An entry without a spec describes a mini.ai built-in.
-- `spec = false` disables mini.ai's object, so Vim's own applies.
local ai = require "mini.ai"
local util = require "utils/init"
local ts = require "utils/treesitter"

---A mini.ai region in the current buffer: 1-indexed lines and byte columns,
---`to` inclusive. An empty region has no `to` (`:h MiniAi-glossary`).
---@alias TextobjectRegion { from: { line: integer, col: integer }, to: { line: integer, col: integer }|nil }

---Region of a tree-sitter range in the current buffer, as
---`MiniAi.gen_spec.treesitter` converts it: a range that ends at column 0
---ends on the newline of the line before.
---@param range Range6 0-indexed, end exclusive
---@return TextobjectRegion
local function range_region(range)
    local to = { line = range[4] + 1, col = range[5] }
    if to.col == 0 then
        to.line = to.line - 1
        to.col = #util.get_line(to.line - 1) + 1
    end
    return { from = { line = range[1] + 1, col = range[2] + 1 }, to = to }
end

---Matches of the `textobjects` query of one language tree that intersect a
---range of rows, each as the regions of its captures, with directives such as
---`#offset!` applied. A quantified capture gives one region from its first
---node to its last. Trees of other languages, parent or injected, are not
---searched.
---@param lang_tree vim.treesitter.LanguageTree of the current buffer
---@param start integer first row, 0-indexed
---@param stop integer row after the last, 0-indexed
---@return table<string, TextobjectRegion>[] matches by capture name without `@`, e.g. "comment.outer"
local function capture_matches(lang_tree, start, stop)
    local query = vim.treesitter.query.get(lang_tree:lang(), "textobjects")
    if not query then return {} end
    local buf = lang_tree:source() --[[@as integer]]
    local matches = {}
    for _, tree in ipairs(lang_tree:trees()) do
        for _, match, metadata in query:iter_matches(tree:root(), buf, start, stop) do
            local regions = {}
            for id, nodes in pairs(match) do
                local first = vim.treesitter.get_range(nodes[1], buf, metadata[id])
                local last = vim.treesitter.get_range(nodes[#nodes], buf, metadata[id])
                regions[query.captures[id]] = range_region {
                    first[1], first[2], first[3], last[4], last[5], last[6],
                }
            end
            matches[#matches + 1] = regions
        end
    end
    return matches
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

---Text of a non-empty `region` of the current buffer, lines joined by "\n".
---A newline at its end is left out.
---@param region TextobjectRegion
---@return string
local function region_text(region)
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
    local query = lang_tree and vim.treesitter.query.get(lang_tree:lang(), "textobjects")
    ---@type { outer: TextobjectRegion, inner: TextobjectRegion|nil }[]
    local comments = {}
    if lang_tree and query and vim.list_contains(query.captures, "comment.outer") then
        local start = reference.from.line - 1 - opts.n_lines
        local stop = reference_to.line + opts.n_lines
        for _, match in ipairs(capture_matches(lang_tree, start, stop)) do
            if match["comment.outer"] then
                comments[#comments + 1] = { outer = match["comment.outer"], inner = match["comment.inner"] }
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

local treesitter = ai.gen_spec.treesitter

---@type table<string, { spec: any, desc: string|nil }>
return {
    f = { spec = treesitter { a = "@function.outer", i = "@function.inner" }, desc = "Function" },
    C = { spec = treesitter { a = "@class.outer", i = "@class.inner" }, desc = "Class" },
    a = { spec = treesitter { a = "@parameter.outer", i = "@parameter.inner" }, desc = "Parameter" },
    ["?"] = { spec = treesitter { a = "@conditional.outer", i = "@conditional.inner" }, desc = "Conditional" },
    o = { spec = treesitter { a = "@loop.outer", i = "@loop.inner" }, desc = "Loop" },
    k = { spec = treesitter { a = "@assignment.lhs", i = "@assignment.lhs" }, desc = "Assignment LHS" },
    v = { spec = treesitter { a = "@assignment.rhs", i = "@assignment.rhs" }, desc = "Assignment RHS" },
    c = { spec = comment, desc = "Comment" },
    [","] = { spec = ai.gen_spec.argument(), desc = "Argument" },
    b = { desc = "Bracket" },
    q = { desc = "Quote" },
    -- mini.ai's `t` pairs an opening tag with the first closing tag of the same
    -- name, so it fails on the outer of nested same-name tags.
    t = { spec = false },
}
