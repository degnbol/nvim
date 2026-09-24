---@diagnostic disable: unused-local
local lsu = require "utils/luasnip"
local s, t, i, c, f, d, sn, fmta, conds, rep, ms = lsu.s, lsu.t, lsu.i, lsu.c, lsu.f, lsu.d, lsu.sn, lsu.fmta, lsu.conds, lsu.rep, lsu.ms
local re = lsu.re

return {
s({trig="glossy", dscr="New glossy entry", condition=conds.line_begin},
fmta([[<>:
  short: <>
  long: <>
  <>
]], {i(1, "key"), lsu.upper_jump(2, 1), i(3, ""), c(4, {t"", t"description: "})})),
}

