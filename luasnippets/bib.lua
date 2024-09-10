local ls = require "luasnip"
local sn = ls.sn
local s = ls.s
local t = ls.t
local i = ls.i
local d = ls.d
local fmta = ls.fmta

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
