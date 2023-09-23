local util = require "utils/luasnip"

return {
--

s({
    trig='```',
    dscr="Tripple backticks",
    snippetType='autosnippet',
},
fmta([[```<>
```
]], i(1))),
s({
    trig='```',
    dscr="Tripple quotes when the second quote is written last. This is relevant when something autopairs an opening quote.",
    snippetType='autosnippet',
    trigEngine=function (trigger) return util.match_ahead(1) end,
},
-- NOTE: there are only two final backticks since one is left over from end of trigger.
fmta([[```<>
``
]], i(1))),

}
