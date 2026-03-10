---@diagnostic disable: unused-local
local ls = require "utils/luasnip"
local s, t, i, c, f, d, sn, fmta, conds, rep, ms = ls.s, ls.t, ls.i, ls.c, ls.f, ls.d, ls.sn, ls.fmta, ls.conds, ls.rep, ls.ms
return {
--
s({trig="priority", dscr="Set priority."},
{t'(#set! "priority" ', i(1, "200"), t')'}),


s({
    -- you can also write "extends", it will be close enough.
    trig=";extends", 
    dscr="File extends previous definitions.",
    show_condition=conds.line_end,
},
t";extends"),


}
