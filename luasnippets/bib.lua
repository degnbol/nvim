---@diagnostic disable: unused-local
local ls = require "utils/luasnip"
local s, t, i, c, f, d, sn, fmta, conds, rep, ms = ls.s, ls.t, ls.i, ls.c, ls.f, ls.d, ls.sn, ls.fmta, ls.conds, ls.rep, ls.ms

---https://github.com/L3MON4D3/LuaSnip/blob/master/DOC.md#dynamicnode
local function _dUpper(args)
    -- the returned snippetNode doesn't need a position; it's inserted
    -- "inside" the dynamicNode.
    return sn(nil, {
        -- jump-indices are local to each snippetNode, so restart at 1.
        i(1, args[1][1]:upper())
    })
end
local function dUpper(jump_index, node_reference)
    return d(jump_index, _dUpper, {node_reference})
end

return {}, {

-- snippets for writing glossary.bib

s({trig="ac", dscr="Acronym", condition=conds.line_begin},
fmta([[@acronym{<>,
    short={<>},
    long={<>}
}
]], {i(1, "LABEL"), dUpper(2, 1), dUpper(3, 1)})),

s({trig="en", dscr="Entry", condition=conds.line_begin},
fmta([[@entry{<>,
    name={<>},
    description={<>}
}
]], {i(1, "LABEL"), dUpper(2, 1), dUpper(3, 1)})),

}
