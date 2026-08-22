-- Colours for the @chem.* captures of the smarts queries, and the extmarks that
-- give an atom the colour of the element it spells. Element identity is node
-- text, and the treesitter highlighter derives a group from the capture name
-- alone, so identity cannot be a capture: queries/smarts/atoms.scm finds the
-- atoms and the lookup below says which element each one is.
--
-- A private namespace, so no colorscheme paints the groups by accident. Element
-- colours come from chem.elements: Jmol's CPK table, corrected per background so
-- every one of them stays legible — see elements.lua.py for the rule and the
-- source. The elements no structure can contain are called out instead of
-- coloured. The reaction arrow is inverted rather than coloured: a solid block
-- reads as a divider between the two sides, and no hue can be mistaken for an
-- element. Bonds and recursion anchors set a style and no colour, so they keep
-- whatever the capture underneath them resolves to: stacked extmarks merge
-- attributes, the higher priority only winning the ones it sets.
local M = {}

local hi = require "utils/highlights"
local util = require "utils/init"
local elements = require "chem/elements"

local ns = vim.api.nvim_create_namespace("chem.elements")

-- Above the grammar's own query, which gives every element node @constant at the
-- default 100. Conflicting attributes at equal priority resolve by set order.
local PRIORITY = 101

--- Every spelling of an element identity, mapped to its highlight group.
--- `Fe`, `fe` and `#26` are all iron. Lookup is exact, which is what keeps
--- `[#0]`, `[#06]` and `[#999]` out: all three parse, none is an element.
--- @return table<string, string> node text -> group
local function element_groups()
    local groups = {}
    for _, element in ipairs(elements) do
        local group = "@chem.element." .. element.name
        groups["#" .. element.z] = group
        -- `H` never nodes as an element_symbol: the same `H` is the atom in the
        -- bare spelling and a count of attached hydrogens everywhere else, both
        -- under (hydrogen), and no text tells those apart. The atom is matched
        -- structurally in after/queries/smarts/highlights.scm instead.
        if element.symbol ~= "H" then groups[element.symbol] = group end
        if element.aromatic then groups[element.symbol:lower()] = group end
    end
    for symbol, z in pairs(elements.suspect) do
        groups[symbol] = "@chem.suspect"
        groups["#" .. z] = "@chem.suspect"
    end
    return groups
end

local groups = element_groups()

--- @alias ChemHighlight [integer, integer, integer, integer, string]
--- start row, start column, end row, end column, highlight group

--- True when a position is within a range, both ends included — the same
--- containment `nvim_buf_get_extmarks` applies to a mark's start, so what
--- `clear` removes is exactly what `paint` puts back.
--- @param range Range4
--- @param row integer
--- @param col integer
--- @return boolean
local function within(range, row, col)
    return not util.before(row, col, range[1], range[2])
        and not util.before(range[3], range[4], row, col)
end

--- Where the atoms of a range are, and which element each one spells.
--- @param buf integer buffer the tree was parsed from
--- @param root TSNode root of a smarts tree
--- @param range Range4
--- @return ChemHighlight[]
local function element_highlights(buf, root, range)
    local query = assert(vim.treesitter.query.get("smarts", "atoms"))
    local highlights = {}
    local erow = math.min(range[3] + 1, vim.api.nvim_buf_line_count(buf))
    for _, node in query:iter_captures(root, buf, range[1], erow) do
        local group = groups[vim.treesitter.get_node_text(node, buf)]
        local r1, c1, r2, c2 = node:range()
        if group and within(range, r1, c1) then
            highlights[#highlights + 1] = { r1, c1, r2, c2, group }
        end
    end
    return highlights
end

--- The part of a tree's range one change touched: the change's rows, but the
--- range's own columns, so that a sibling region on the same row keeps its
--- marks. Whole rows because a change is bytes, while the atom a mark colours is
--- the unit being repainted; a region spanning one row — a TSV cell — is
--- therefore repainted whole.
--- @param range Range4
--- @param change Range6
--- @return Range4|nil touched nil when the change is outside the range
local function changed_part(range, change)
    local srow, scol, erow, ecol = range[1], range[2], range[3], range[4]
    if change[4] < srow or erow < change[1] then return end
    if change[1] > srow then srow, scol = change[1], 0 end
    -- The start of the row after the change is the end of the change's last row.
    if change[4] < erow then erow, ecol = change[4] + 1, 0 end
    return { srow, scol, erow, ecol }
end

--- Drop the element marks of a range, and any invalidated mark just past it. By
--- range and not by row, because two structures can share a row — a TSV row with
--- two chemical columns — as separate regions, and a row-wide clear would wipe
--- the other one's marks.
---
--- The invalid ones are the marks of text an edit deleted wholesale:
--- `nvim_buf_set_lines` over a line collapses every mark on it to the start of
--- the next line, outside the columns of every region on the line it came from.
--- That position is indistinguishable from a live mark in the first column of
--- the next row, so 'invalidate' is what tells the two apart.
--- @param buf integer
--- @param range Range4
local function clear(buf, range)
    -- A root tree's range ends at a sentinel row past any buffer, and
    -- nvim_buf_get_extmarks answers a row that large with a single mark instead
    -- of clamping it: 2^32-1 returns one where 2^31-1 returns all of them.
    local erow = math.min(range[3] + 1, vim.api.nvim_buf_line_count(buf))
    for _, mark in ipairs(vim.api.nvim_buf_get_extmarks(
        buf, ns, { range[1], range[2] }, { erow, 0 }, { details = true })) do
        if assert(mark[4]).invalid or within(range, mark[2], mark[3]) then
            vim.api.nvim_buf_del_extmark(buf, ns, mark[1])
        end
    end
end

--- Drop the element marks of every smarts region at or below a language tree.
--- Off the regions and not the trees, because `invalidate(true)` clears the
--- trees of every tree it is about to reparse: a language dropped by that
--- reparse then reports no tree to take the rows from, and its regions are the
--- only record left of where it was.
--- @param buf integer
--- @param ltree vim.treesitter.LanguageTree
local function clear_regions(buf, ltree)
    if ltree:lang() == "smarts" then
        for _, region in ipairs(ltree:included_regions()) do
            for _, range in ipairs(region) do
                clear(buf, { range[1], range[2], range[4], range[5] })
            end
        end
    end
    for _, child in pairs(ltree:children()) do clear_regions(buf, child) end
end

--- Replace the element marks of a range with the atoms one tree holds there.
--- @param buf integer
--- @param root TSNode root of a smarts tree
--- @param range Range4
local function paint(buf, root, range)
    clear(buf, range)
    for _, mark in ipairs(element_highlights(buf, root, range)) do
        vim.api.nvim_buf_set_extmark(buf, ns, mark[1], mark[2], {
            end_row = mark[3], end_col = mark[4],
            hl_group = mark[5], priority = PRIORITY,
            -- Not to hide the mark, which has nothing left to colour anyway, but
            -- so that `clear` can recognise one whose text is gone.
            invalidate = true,
        })
    end
end

--- Repaint the part of one tree's ranges that a parse reported changed.
--- @param buf integer
--- @param tree TSTree the tree as just parsed
--- @param changes Range6[]
local function repaint(buf, tree, changes)
    for _, range in ipairs(tree:included_ranges(false)) do
        for _, change in ipairs(changes) do
            -- A newly created tree reports the whole document, spelled as a
            -- sentinel end row (highlighter.lua reads the same one), so the
            -- tree's own ranges are what bound the work: a sibling region keeps
            -- its marks.
            local changed = changed_part(range, change)
            if changed then paint(buf, tree:root(), changed) end
        end
    end
end

--- Call a function on a language tree and every tree injected below it.
--- @param ltree vim.treesitter.LanguageTree
--- @param fn fun(ltree: vim.treesitter.LanguageTree)
local function walk(ltree, fn)
    fn(ltree)
    for _, child in pairs(ltree:children()) do walk(child, fn) end
end

--- Drop every element mark in a buffer, for a change in what counts as a
--- structure rather than a change in the text. Pair it with a reparse, which
--- repaints the regions that are still chemical: a region that stopped being one
--- is gone from the tree, and nothing left reports its rows.
--- @param buf integer
function M.drop_marks(buf)
    vim.api.nvim_buf_clear_namespace(buf, ns, 0, -1)
end

--- Colour the atoms of a buffer, and recolour them after every later reparse.
--- @param buf integer
function M.attach(buf)
    local parser = assert(vim.treesitter.get_parser(buf))

    --- @param ltree vim.treesitter.LanguageTree
    local function watch(ltree)
        if ltree:lang() ~= "smarts" then return end
        -- Per language tree and not recursive, because on_changedtree says which
        -- rows changed but not which language did, and it fires before the tree
        -- is stored — so the tree it is handed is the only reach to the new
        -- parse. Not on_lines either: that fires before the reparse, when the
        -- tree still describes the old text. And no debounce — on_changedtree
        -- already follows a parse, and parses are async and range-limited.
        ltree:register_cbs {
            on_changedtree = function(changes, tree) repaint(buf, tree, changes) end,
        }
        -- Whatever is already parsed reports no change, and a buffer nvim has
        -- drawn once is parsed before any FileType handler runs.
        for _, tree in pairs(ltree:trees()) do
            for _, range in ipairs(tree:included_ranges(false)) do
                paint(buf, tree:root(), range)
            end
        end
    end

    walk(parser, watch)
    -- A fence's structure is a child added by a later parse, and on_child_added
    -- fires before that child first parses. Recursive, so a fence inside an
    -- injection is covered too.
    parser:register_cbs({
        on_child_added = watch,
        -- The child is unlinked before this fires and nothing will report its
        -- rows again, so its marks have to go here.
        on_child_removed = function(child) clear_regions(buf, child) end,
    }, true)
end

--- Define the groups, and attach to every buffer treesitter highlights.
function M.setup()
    -- Fired by plugin/treesitter.lua once a buffer's parser exists, so attach
    -- inherits which buffers get highlighted at all without restating the policy.
    vim.api.nvim_create_autocmd("User", {
        pattern = "TSHighlightStart",
        group = vim.api.nvim_create_augroup("chem.elements", { clear = true }),
        callback = function(args) M.attach(args.data.buf) end,
    })
    hi.onColorScheme(function()
        hi.set("@chem.reaction", { reverse = true })
        hi.set("@chem.bond", { bold = true })
        hi.set("@chem.anchor", { underline = true })
        hi.link("@chem.suspect", "DiagnosticWarn")
        -- 'background' is set by the colorscheme before ColorScheme fires, and
        -- names the column it wants.
        for _, element in ipairs(elements) do
            hi.setfg("@chem.element." .. element.name, element[vim.o.background])
        end
    end)
end

return M
