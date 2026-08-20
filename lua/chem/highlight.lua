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

--- Where the atoms of a row range are, and which element each one spells.
--- @param buf integer buffer the tree was parsed from
--- @param root TSNode root of a smarts tree
--- @param srow integer first row to search, 0-based
--- @param erow integer row to stop before
--- @return ChemHighlight[]
local function element_highlights(buf, root, srow, erow)
    local query = assert(vim.treesitter.query.get("smarts", "atoms"))
    local highlights = {}
    for _, node in query:iter_captures(root, buf, srow, erow) do
        local group = groups[vim.treesitter.get_node_text(node, buf)]
        if group then
            local r1, c1, r2, c2 = node:range()
            highlights[#highlights + 1] = { r1, c1, r2, c2, group }
        end
    end
    return highlights
end

--- The rows one tree covers, as a range clamped to the buffer.
--- @param buf integer
--- @param tree TSTree
--- @return integer srow first row, 0-based
--- @return integer erow row after the last
local function tree_rows(buf, tree)
    local ranges = tree:included_ranges(false)
    return ranges[1][1],
        math.min(ranges[#ranges][3] + 1, vim.api.nvim_buf_line_count(buf))
end

--- Replace the element marks of a row range with one tree's.
--- @param buf integer
--- @param tree TSTree
--- @param srow integer first row to paint, 0-based
--- @param erow integer row to stop before
local function paint(buf, tree, srow, erow)
    vim.api.nvim_buf_clear_namespace(buf, ns, srow, erow)
    for _, mark in ipairs(element_highlights(buf, tree:root(), srow, erow)) do
        vim.api.nvim_buf_set_extmark(buf, ns, mark[1], mark[2], {
            end_row = mark[3], end_col = mark[4],
            hl_group = mark[5], priority = PRIORITY,
        })
    end
end

--- Repaint the rows a parse reported changed in one tree.
--- @param buf integer
--- @param tree TSTree the tree as just parsed
--- @param changes Range6[]
local function repaint(buf, tree, changes)
    local tstart, tend = tree_rows(buf, tree)
    for _, change in ipairs(changes) do
        -- A newly created tree reports the whole document, spelled as a sentinel
        -- end row (highlighter.lua reads the same one), so the tree's own rows
        -- are what bound the work: a sibling region keeps its marks.
        local srow, erow = math.max(change[1], tstart), math.min(change[4] + 1, tend)
        if srow < erow then paint(buf, tree, srow, erow) end
    end
end

--- Call a function on a language tree and every tree injected below it.
--- @param ltree vim.treesitter.LanguageTree
--- @param fn fun(ltree: vim.treesitter.LanguageTree)
local function walk(ltree, fn)
    fn(ltree)
    for _, child in pairs(ltree:children()) do walk(child, fn) end
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
            paint(buf, tree, tree_rows(buf, tree))
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
        on_child_removed = function(child)
            child:for_each_tree(function(tree, ltree)
                if ltree:lang() ~= "smarts" then return end
                local srow, erow = tree_rows(buf, tree)
                vim.api.nvim_buf_clear_namespace(buf, ns, srow, erow)
            end)
        end,
    }, true)
end

--- Define the groups, and redefine them on every later colorscheme change.
function M.setup()
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
