local util = require "utils/luasnip"

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
