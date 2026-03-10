---@diagnostic disable: unused-local
local lsu = require "utils/luasnip"
local s, t, i, c, f, d, sn, fmta, conds, rep, ms = lsu.s, lsu.t, lsu.i, lsu.c, lsu.f, lsu.d, lsu.sn, lsu.fmta, lsu.conds, lsu.rep, lsu.ms
local re = lsu.re

return {
}, {
-- autosnippets

s({trig="(%s*)([%w-]+):", dscr="Auto quote unquoted key.",
trigEngine="pattern"},
{re(1), t'"', re(2), t'": '}),

s({trig=": ?{{", dscr="Opening bracket after a key.",
trigEngine="pattern"},
-- two spaces
{t{": {", "  "}, i(1), t{"", "}"}}),

}
