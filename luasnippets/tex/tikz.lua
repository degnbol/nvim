---@diagnostic disable: unused-local
local ls = require "utils/luasnip"
local s, t, i, c, f, d, sn, fmta, conds, rep, ms = ls.s, ls.t, ls.i, ls.c, ls.f, ls.d, ls.sn, ls.fmta, ls.conds, ls.rep, ls.ms
-- TODO: only complete in tikz env
return {
s({trig="node", dscr="tikz node"},
fmta([[\node<><> at (<>) {<>};
]], {
    i(1, " (NEW ID)"),
    i(2, " [blue,fill=red]"),
    i(3, "0,0"),
    i(4, "TEXT"),
})),
}
