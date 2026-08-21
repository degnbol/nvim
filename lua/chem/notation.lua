-- Which column names hold a chemical line notation, and which notation each of
-- them names. An exact set rather than a pattern: a pattern would catch the
-- spellings not listed for free, but eventually matches a column holding
-- something else, and a column the set misses can be marked by hand
-- (chem/tsv.lua).
local M = {}

--- @alias ChemNotation "smiles"|"smarts"

--- The notation a bare name spells. SMIRKS is a reaction transform written in
--- SMARTS, so the same parser reads it.
--- @type table<string, ChemNotation>
local BARE = {
    smiles = "smiles",
    cxsmiles = "smiles",
    smarts = "smarts",
    smirks = "smarts",
}

--- The roles a structure column takes in a reaction table, as name prefixes.
--- @type string[]
local ROLES = { "", "reaction_", "product_", "reactant_" }

--- @type table<string, ChemNotation>
local NOTATIONS = {}
for name, notation in pairs(BARE) do
    for _, role in ipairs(ROLES) do NOTATIONS[role .. name] = notation end
end

--- The chemical line notation a column name spells, or nil for a name that
--- spells none. Case-insensitive: `SMILES` is as common a heading as `smiles`.
--- @param name string
--- @return ChemNotation|nil
function M.of(name)
    return NOTATIONS[name:lower()]
end

return M
