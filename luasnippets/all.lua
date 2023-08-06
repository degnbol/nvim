local util = require "utils/luasnip"

return {
-- for all filetypes

-- alt to autopairs
s({trig='""', dscr="autopairs", snippetType='autosnippet'},
{t'"', i(1), t'"'}),
s({trig="''", dscr="autopairs", snippetType='autosnippet'},
{t"'", i(1), t"'"}),
s({trig="``", dscr="autopairs", snippetType='autosnippet'},
{t"`", i(1), t"`"}),
s({trig="((", dscr="autopairs", snippetType='autosnippet', wordTrig=false},
{t"(", i(1), t")"}),
s({trig="{{", dscr="autopairs", snippetType='autosnippet', wordTrig=false, priority=100},
{t"{", i(1), t"}"}),
s({
    trig='{{}',
    dscr="Tripple press open bracket. Second bracket is typed last. The last one turns the other way due to the snippet above that turns {{ into {}",
    snippetType='autosnippet',
    trigEngine=function (trigger) return util.match_oneAhead end,
},
-- NOTE: there is no final bracket in replacement since one is left over from end of trigger.
fmta([[{
    <>

]], i(1))),

}
