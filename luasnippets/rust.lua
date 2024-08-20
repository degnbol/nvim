-- should be unnessary but s keeps being understood as a bool set somewhere
local ls = require "luasnip"
local s = ls.s
local t = ls.t
local i = ls.i

return {
-- not auto

}, {
-- auto

s({trig="?", dscr="Debug print", condition=conds.line_begin},
{t'println!("{:?}", ', i(1, "VALUE"), t');'}),

}
