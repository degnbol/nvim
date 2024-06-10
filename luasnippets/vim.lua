local ls = require "luasnip"
local s = ls.s

return {
-- not auto

s({trig="esc", dscr="Temp ESC"},
{t[[<C-\><C-O>]]}),

}, {
-- auto

-- Since it is a wordTrig it will not get corrected in valid use, such as pathogen#infect()
-- TODO: make a condition for being inside a comment and don't exclude that with condition=-in_comment
s({trig="#", dscr="Correct wrong comment char."},
{t'"'}),

}
