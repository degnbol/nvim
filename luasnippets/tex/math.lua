local lsu = require"luasnip_util"
local in_math = lsu.in_math
local re = lsu.re
local get_visual = lsu.get_visual

local optm = {condition=in_math}

-- all snippets here are expected to be called within a math env.
-- Any snippet that may have parts replaced with unicode are in unicode.lua and 
-- duplicated in nouniode.lua
return {
-- not autosnippet

s({trig="lr", dscr="left/right", snippetType='autosnippet'},
{t"\\left", i(1, '('), d(2, get_visual), t"\\right",
f(function(val)
    -- first char only since the node technically contains all the remaining 
    -- text inside the \left and \right
    val = val[1][1]:sub(1,1)
    if val == "(" then
        return ")"
    elseif val == "|" then
        return "|"
    elseif val == "[" then
        return "]"
    elseif val == "<" then
        return ">"
    end
end, {1}),
}),
-- doesn't work if autopairing
s({trig="()", dscr="left( right)", wordTrig=false, snippetType="autosnippet"},
{t"\\left(", d(1, get_visual), t"\\right)"}, optm),

},
{
-- autosnippet

s({trig="case", dscr="cases"},
fmta([[
\begin{cases}
	<>
\end{cases}
]], {i(1)}), optm),

s({trig="(%a)(%d)", dscr="auto subscript", regTrig=true},
{re(1), t("_"), re(2)}, optm),
s({trig="(%a)_(%d%d)", dscr="auto subscript 2", regTrig=true},
{re(1), t("_{"), re(2), t("}") }, optm),

s({trig="//", dscr="fraction"},
fmta("\\frac{<>}{<>}", {d(1, get_visual), i(2)}), optm),

s({trig='(%b())/', dscr="fraction after paren", regTrig=true},
{t"\\frac{",
f(function(_, snip)
    s = snip.captures[1]
    return string.sub(s, 2, #s-1)
end),
t"}{", i(1), t"}"}, optm),

s({trig='\\left(%b())/', dscr="fraction after paren 2", regTrig=true},
{t"\\frac{",
f(function(_, snip)
    s = snip.captures[1]
    -- 7 for length of "\right)"
    return string.sub(s, 2, #s-7)
end),
t"}{", i(1), t"}"}, optm),

-- [^%s(){}]+ instead of %w+ since the latter doesn't capture unicode.
-- I'm not excluding + and * etc. on purpose, so I can type "a+b/" to get a 
-- division (a+b)/ and type "a + b/" to get a division only on b
s({trig="([^%s(){}]+)/", dscr="fraction after word", regTrig=true},
{t"\\frac{", re(1), t"}{", i(1), t"}"}, optm),

s({trig='([^%a])ee', regTrig=true, wordTrig=false},
fmta("<>e^{<>}", { re(1), d(1, get_visual) }), optm),
s({trig = 'invs', dscr="inverse"},
t"^{-1}", optm),

-- it's ok that they are wordTrig, e.g. {} are considered word separators, so they autocomplete immediately after these.
-- [^%s(){}0-9]+ instead of %a+ since the latter doesn't capture unicode.
s({trig="([^%s(){}0-9]+)hat", dscr="hat postfix", regTrig=true},
{t"\\hat{", re(1), t"}"}, optm),
s({trig="hat", dscr="hat"},
{t"\\hat{", i(1), t"}"}, optm),
-- [^%s(){}0-9]+ instead of %a+ since the latter doesn't capture unicode.
s({trig="([^%s(){}0-9]+)bar", dscr="bar postfix", regTrig=true},
{t"\\overline{", re(1), t"}"}, optm),
s({trig="bar", dscr="bar"},
{t"\\overline{", i(1), t"}"}, optm),
-- [^%s(){}0-9]+ instead of %a+ since the latter doesn't capture unicode.
-- There's also \widetilde that becomes longer for multiple letters or wide letters like M.
s({trig="([^%s(){}0-9]+)~", dscr="tilde postfix", regTrig=true},
{t"\\tilde{", re(1), t"}"}, optm),

s({trig="norm", wordTrig=false},
{t"\\|", i(1), t"\\|"}, optm),

s({trig="sr", dscr="sqrt", wordTrig=false},
{t"\\sqrt{", d(1, get_visual), t"}"}, optm),
s({trig="sq", descr="squared", wordTrig=false},
t"^2", optm),
s({trig="cb", descr="cubed", wordTrig=false},
t"^3", optm),

s({trig="td", descr="to the power", wordTrig=false},
{t"^{", i(1), t"}"}, optm),
s({trig="^^", descr="superscript", wordTrig=false},
{t"^{", i(1), t"}"}, optm),
s({trig="tD", descr="to the (power)", wordTrig=false},
{t"^{(", i(1), t")}"}, optm),
s({trig="__", descr="subscript", wordTrig=false},
{t"_{", i(1), t"}"}, optm),

s({trig="=>", descr="implies", wordTrig=false},
{t"\\implies"}, optm),
s({trig="=<", descr="implied by", wordTrig=false},
{t"\\impliedby"}, optm),

s({trig="==", descr="equal", wordTrig=false},
t"&= ", optm),

s({trig=[[\\\]], descr="set minus", wordTrig=false},
t"\\setminus", optm),

s({trig="RR", descr="real", wordTrig=false},
t"\\R", optm),
s({trig="ZZ", descr="real", wordTrig=false},
t"\\Z", optm),

s({trig="tt", descr="text", wordTrig=false},
{t"\\text{",i(1),t"}"}, optm),

s({trig="\\?SI", descr="SI", regTrig=true},
{t"\\SI{",i(1),t"}{", i(2), t"}"}, optm),

}
