local Kind = vim.lsp.protocol.CompletionItemKind

local SELECTORS = {
    'name', 'resn', 'resi', 'chain', 'segi', 'alt', 'model',
    'index', 'id', 'rank', 'pepseq', 'label', 'elem', 'flag',
    'ss', 'rep', 'color', 'cartoon_color', 'ribbon_color',
    'numeric_type', 'formal_charge', 'partial_charge', 'state',
    'b', 'q', 'x', 'y', 'z',
}

local SELECTOR_SHORTHANDS = {
    'n.', 'r.', 'i.', 'c.', 's.', 'm.', 'idx.', 'ps.',
    'e.', 'f.', 'nt.', 'pc.', 'fc.', 'v.', 'pr.',
}

local BUILTINS = {
    'all', 'none', 'enabled', 'visible', 'bonded', 'protected',
    'fixed', 'restrained', 'masked', 'organic', 'inorganic',
    'solvent', 'polymer', 'guide', 'hetatm', 'hydrogens',
    'backbone', 'sidechain', 'metals', 'donors', 'acceptors',
    'present', 'center', 'origin',
    'polymer.protein', 'polymer.nucleic',
}

local BUILTIN_SHORTHANDS = {
    'org.', 'ino.', 'sol.', 'pol.', 'h.', 'bb.', 'sc.',
    'don.', 'acc.', 'fxd.', 'rst.', 'msk.',
}

local EXPANSION_KEYWORDS = {
    'byres', 'bychain', 'bymolecule', 'byfragment', 'bysegi',
    'byobject', 'bycalpha', 'byring', 'bycell', 'bound_to',
    'neighbor', 'extend', 'first', 'last', 'in', 'like',
}

local EXPANSION_SHORTHANDS = {
    'br.', 'bc.', 'bm.', 'bf.', 'bs.', 'bca.', 'bto.',
    'nbr.', 'xt.',
}

local PROXIMITY_KEYWORDS = {
    'within', 'around', 'expand', 'gap', 'near_to', 'beyond',
}

local PROXIMITY_SHORTHANDS = {
    'w.', 'a.', 'x.', 'nto.', 'be.',
}

local OPERATORS = { 'and', 'or', 'not', 'of' }

local REPRESENTATIONS = {
    'everything', 'lines', 'sticks', 'spheres', 'dots', 'surface', 'mesh',
    'cartoon', 'ribbon', 'labels', 'nb_spheres', 'nonbonded', 'volume',
    'slice', 'licorice', 'cell', 'extent',
}

local function make_items(words, kind, detail)
    local items = {}
    for _, w in ipairs(words) do
        items[#items + 1] = {
            label = w,
            kind = kind,
            labelDetails = { description = detail },
            source = "pymol_select",
        }
    end
    return items
end

-- Check if the pymol_select injection is active at the cursor position.
-- Piggybacks on the injection scoping (function args, assignments — not docstrings).
local function in_pymol_context()
    local ok, parser = pcall(vim.treesitter.get_parser, 0)
    if not ok or not parser then return false end
    local pos = vim.api.nvim_win_get_cursor(0)
    local range = { pos[1] - 1, pos[2], pos[1] - 1, pos[2] }
    return parser:language_for_range(range):lang() == "pymol_select"
end

local M = {}

function M.new()
    local self = setmetatable({}, { __index = M })
    self.items = {}
    vim.list_extend(self.items, make_items(SELECTORS, Kind.Keyword, "selector"))
    vim.list_extend(self.items, make_items(SELECTOR_SHORTHANDS, Kind.Keyword, "selector"))
    vim.list_extend(self.items, make_items(BUILTINS, Kind.Constant, "builtin"))
    vim.list_extend(self.items, make_items(BUILTIN_SHORTHANDS, Kind.Constant, "builtin"))
    vim.list_extend(self.items, make_items(EXPANSION_KEYWORDS, Kind.Keyword, "expansion"))
    vim.list_extend(self.items, make_items(EXPANSION_SHORTHANDS, Kind.Keyword, "expansion"))
    vim.list_extend(self.items, make_items(PROXIMITY_KEYWORDS, Kind.Keyword, "proximity"))
    vim.list_extend(self.items, make_items(PROXIMITY_SHORTHANDS, Kind.Keyword, "proximity"))
    vim.list_extend(self.items, make_items(OPERATORS, Kind.Operator, "operator"))
    vim.list_extend(self.items, make_items(REPRESENTATIONS, Kind.Enum, "representation"))
    return self
end

function M:get_completions(_, callback)
    if not in_pymol_context() then
        callback({ is_incomplete_forward = true, is_incomplete_backward = true, items = {} })
        return function() end
    end
    callback({
        is_incomplete_forward = false,
        is_incomplete_backward = false,
        items = self.items,
    })
    return function() end
end

return M
