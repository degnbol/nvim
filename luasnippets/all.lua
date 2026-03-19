---@diagnostic disable: unused-local
local ls = require "utils/luasnip"
local s, t, i, c, f, d, sn, fmta, conds, rep, ms = ls.s, ls.t, ls.i, ls.c, ls.f, ls.d, ls.sn, ls.fmta, ls.conds, ls.rep, ls.ms
return {
-- for all filetypes

-- typos that can't be corrected with abb
s({trig="i..e", dscr="I.e. typo", snippetType='autosnippet'}, {t"i.e."}),
s({trig="e..g", dscr="E.g. typo", snippetType='autosnippet'}, {t"e.g."}),
s({trig="i.e..", dscr="I.e. typo", snippetType='autosnippet'}, {t"i.e."}),
s({trig="e.g..", dscr="E.g. typo", snippetType='autosnippet'}, {t"e.g."}),
s({trig="e.g ", dscr="E.g. typo", snippetType='autosnippet'}, {t"e.g. "}),
s({trig="ya'll", dscr="You all", snippetType='autosnippet'}, {t"y'all"}),

}
