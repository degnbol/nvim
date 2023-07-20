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
s({trig="{{", dscr="autopairs", snippetType='autosnippet', wordTrig=false},
{t"{", i(1), t"}"}),



}
