local util = require "utils/luasnip"

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
