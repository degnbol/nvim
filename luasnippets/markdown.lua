local ls = require "luasnip"
local s = ls.s
local t = ls.t
local i = ls.i
local conds = require("luasnip.extras.expand_conditions")
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
