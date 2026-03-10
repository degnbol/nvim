---@diagnostic disable: unused-local
local lsu = require "utils/luasnip"
local s, t, i, c, f, d, sn, fmta, conds, rep, ms = lsu.s, lsu.t, lsu.i, lsu.c, lsu.f, lsu.d, lsu.sn, lsu.fmta, lsu.conds, lsu.rep, lsu.ms
local re = lsu.re

return {
    ms({
        "^([\t ]*)(%w+) ?=",
        "^([\t ]*)set (%w+) ?=",
        common={
        trigEngine="pattern",
        snippetType='autosnippet',
        dscr="Variable assignment",
    }},
    {re(1), t"set ", re(2), t" to"}),
}
