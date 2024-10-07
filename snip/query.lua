local ls = require "luasnip"
local s = ls.s

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
