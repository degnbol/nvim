---@diagnostic disable: unused-local
local ls = require "utils/luasnip"
local s, t, i, c, f, d, sn, fmta, conds, rep, ms = ls.s, ls.t, ls.i, ls.c, ls.f, ls.d, ls.sn, ls.fmta, ls.conds, ls.rep, ls.ms
return {
--

-- NOTE: requires concealcursor-=v if conceallevel>0
s({
    trig='````',
    dscr="Tripple backticks",
    snippetType='autosnippet',
    condition=conds.line_begin,
},
{t"```", i(1, "LANGUAGE"), t{"", "```"}}),


}
