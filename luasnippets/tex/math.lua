local lsu = require"luasnip_util"
local in_math = lsu.in_math
local re = lsu.re
local get_visual = lsu.get_visual

return {
-- all snippets here are expected to be called within a math env.
-- Any snippet that may have parts replaced with unicode are in unicode.lua and 
-- duplicated in nouniode.lua

s({trig="case", dscr="cases", snippetType="autosnippet"},
  fmta(
[[
\begin{cases}
	<>
\end{cases}
]],
    { i(1), }
  ),
  {condition=in_math}
),

s(
{trig="(%a)(%d)", dscr="auto subscript", regTrig=true, snippetType="autosnippet"},
{re(1), t("_"), re(2)},
{condition=in_math}
),
s(
{trig="(%a)_(%d%d)", dscr="auto subscript 2", regTrig=true, snippetType="autosnippet"},
{re(1), t("_{"), re(2), t("}") },
{condition=in_math}
),

s(
{trig="//", dscr="fraction", snippetType="autosnippet"},
fmta("\\frac{<>}{<>}", {d(1, get_visual), i(2)}),
{condition=in_math}
),

s(
{trig='(%b())/', dscr="fraction after paren", regTrig=true, snippetType="autosnippet"},
{t("\\frac{"),
f(function(_, snip)
    s = snip.captures[1]
    return string.sub(s, 2, #s-1)
end),
t("}{"), i(1), t("}")},
{condition=in_math}
),

s(
{trig='\\left(%b())/', dscr="fraction after paren 2", regTrig=true, snippetType="autosnippet"},
{t("\\frac{"),
f(function(_, snip)
    s = snip.captures[1]
    -- 7 for length of "\right)"
    return string.sub(s, 2, #s-7)
end),
t("}{"), i(1), t("}")},
{condition=in_math}
),

-- [^%s(){}]+ instead of %w+ since the latter doesn't capture unicode.
-- I'm not excluding + and * etc. on purpose, so I can type "a+b/" to get a 
-- division (a+b)/ and type "a + b/" to get a division only on b
s(
{trig="([^%s(){}]+)/", dscr="fraction after word", regTrig=true, snippetType="autosnippet"},
{t"\\frac{", re(1), t"}{", i(1), t"}"},
{condition=in_math}
),

s({trig = '([^%a])ee', regTrig = true, wordTrig=false, snippetType="autosnippet"},
  fmta(
    "<>e^{<>}",
    { re(1), d(1, get_visual) }
  ),
  {condition=in_math}
),
s({trig = 'invs', dscr="inverse", snippetType="autosnippet"},
t"^{-1}",
{condition=in_math}
),

-- it's ok that they are wordTrig, e.g. {} are considered word separators, so they autocomplete immediately after these.
s({trig="(%a)hat", dscr="hat postfix", regTrig=true, snippetType="autosnippet"},
{t"\\hat{", re(1), t"}"},
{condition=in_math}),
s({trig="hat", dscr="hat", snippetType="autosnippet"},
{t"\\hat{", i(1), t"}"},
{condition=in_math}),

s({trig="(%a)bar", dscr="bar postfix", regTrig=true, snippetType="autosnippet"},
{t"\\overline{", re(1), t"}"},
{condition=in_math}),
s({trig="bar", dscr="bar", snippetType="autosnippet"},
{t"\\overline{", i(1), t"}"},
{condition=in_math}),

s({trig="()", dscr="left( right)", wordTrig=false, snippetType="autosnippet"},
{t"\\left(", d(1, get_visual), t"\\right)"},
{condition=in_math}),
s({trig="lr", dscr="left( right)", wordTrig=false},
{t"\\left(", d(1, get_visual), t"\\right)"},
{condition=in_math}),
s({trig="lr(", dscr="left( right)", wordTrig=false},
{t"\\left(", d(1, get_visual), t"\\right)"},
{condition=in_math}),
s({trig="lr|", dscr="left| right|", wordTrig=false},
{t"\\left|", d(1, get_visual), t"\\right|"},
{condition=in_math}),
s({trig="lr{", dscr="left{ right}", wordTrig=false},
{t"\\left{", d(1, get_visual), t"\\right}"},
{condition=in_math}),
s({trig="lrb", dscr="left{ right}", wordTrig=false},
{t"\\left{", d(1, get_visual), t"\\right}"},
{condition=in_math}),
s({trig="lr[", dscr="left[ right]", wordTrig=false},
{t"\\left[", d(1, get_visual), t"\\right]"},
{condition=in_math}),
s({trig="lr<", dscr="left< right>", wordTrig=false},
{t"\\left<", d(1, get_visual), t"\\right>"},
{condition=in_math}),
s({trig="lra", dscr="left< right>", wordTrig=false},
{t"\\left<", d(1, get_visual), t"\\right>"},
{condition=in_math}),

s({trig="norm", wordTrig=false, snippetType="autosnippet"},
{t"\\|", i(1), t"\\|"},
{condition=in_math}),

s({trig="sr", dscr="sqrt", wordTrig=false, snippetType="autosnippet"},
{t"\\sqrt{", d(1, get_visual), t"}"},
{condition=in_math}
),
s({trig="sq", descr="squared", wordTrig=false, snippetType="autosnippet"},
t"^2",
{condition=in_math}),

s({trig="cb", descr="cubed", wordTrig=false, snippetType="autosnippet"},
t"^3",
{condition=in_math}),

s({trig="td", descr="to the power", wordTrig=false, snippetType="autosnippet"},
{t"^{", i(1), t"}"},
{condition=in_math}),
s({trig="^^", descr="superscript", wordTrig=false, snippetType="autosnippet"},
{t"^{", i(1), t"}"},
{condition=in_math}),
s({trig="tD", descr="to the (power)", wordTrig=false, snippetType="autosnippet"},
{t"^{(", i(1), t")}"},
{condition=in_math}),
s({trig="__", descr="subscript", wordTrig=false, snippetType="autosnippet"},
{t"_{", i(1), t"}"},
{condition=in_math}),

s({trig="=>", descr="implies", wordTrig=false, snippetType="autosnippet"},
{t"\\implies"},
{condition=in_math}),
s({trig="=<", descr="implied by", wordTrig=false, snippetType="autosnippet"},
{t"\\impliedby"},
{condition=in_math}),

s({trig="==", descr="equal", wordTrig=false, snippetType="autosnippet"},
{t"&= ", i(1)},
{condition=in_math}),

s({trig=[[\\\]], descr="set minus", wordTrig=false, snippetType="autosnippet"},
t"\\setminus",
{condition=in_math}),

s({trig="RR", descr="real", wordTrig=false, snippetType="autosnippet"},
t"\\R",
{condition=in_math}),

s({trig="ZZ", descr="real", wordTrig=false, snippetType="autosnippet"},
t"\\Z",
{condition=in_math}),

s({trig="tt", descr="text", wordTrig=false, snippetType="autosnippet"},
{t"\\text{",i(1),t"}"},
{condition=in_math}),

s({trig="SI", descr="SI", wordTrig=false, snippetType="autosnippet"},
{t"\\SI{",i(1),t"}{", i(2), t"}"},
{condition=in_math}),

}
