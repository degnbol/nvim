---@diagnostic disable: unused-local
local ls = require "utils/luasnip"
local s, t, i, c, f, d, sn, fmta, conds, rep, ms = ls.s, ls.t, ls.i, ls.c, ls.f, ls.d, ls.sn, ls.fmta, ls.conds, ls.rep, ls.ms
return {
-- not auto

}, {
-- auto

s({trig="?", dscr="Debug print", condition=conds.line_begin},
{t'println!("{:?}", ', i(1, "VALUE"), t');'}),

}
