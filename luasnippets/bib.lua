---@diagnostic disable: unused-local
local ls = require "utils/luasnip"
local s, t, i, c, f, d, sn, fmta, conds, rep, ms = ls.s, ls.t, ls.i, ls.c, ls.f, ls.d, ls.sn, ls.fmta, ls.conds, ls.rep, ls.ms

return {}, {

-- snippets for writing glossary.bib

s({trig="ac", dscr="Acronym", condition=conds.line_begin},
fmta([[@acronym{<>,
    short={<>},
    long={<>}
}
]], {i(1, "LABEL"), ls.upper_jump(2, 1), ls.upper_jump(3, 1)})),

s({trig="en", dscr="Entry", condition=conds.line_begin},
fmta([[@entry{<>,
    name={<>},
    description={<>}
}
]], {i(1, "LABEL"), ls.upper_jump(2, 1), ls.upper_jump(3, 1)})),

}
