---@diagnostic disable: unused-local
local lsu = require "utils/luasnip"
local s, t, i, c, f, d, sn, fmta, conds, rep, ms = lsu.s, lsu.t, lsu.i, lsu.c, lsu.f, lsu.d, lsu.sn, lsu.fmta, lsu.conds, lsu.rep, lsu.ms


return {
s("template", lsu.putfilenode({"preamble"}, "typst", "typ"), {show_condition = conds.line_end}),
}
