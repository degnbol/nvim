local lsu = require"luasnip_util"
local in_math = lsu.in_math
local re = lsu.re
local get_visual = lsu.get_visual

-- conditions can be set from either first arg of s or third, but shouldn't be set from both.
-- Source: :h luasnip-snippets
optm = {condition=in_math}

return {},{
-- all snippets here are expected to be called within a math env

s({trig=";q", wordTrig=false}, { t"θ" }),
s({trig=";w", wordTrig=false}, { t"ω" }),
s({trig=";e", wordTrig=false}, { t"ε" }), -- using \varepsilon over ϵ
s({trig=";r", wordTrig=false}, { t"ρ" }),
s({trig=";t", wordTrig=false}, { t"τ" }),
s({trig=";y", wordTrig=false}, { t"ψ" }),
s({trig=";u", wordTrig=false}, { t"υ" }),
s({trig=";i", wordTrig=false}, { t"ι" }),
s({trig=";p", wordTrig=false}, { t"π" }),
s({trig=";a", wordTrig=false}, { t"α" }),
s({trig=";s", wordTrig=false}, { t"σ" }),
s({trig=";d", wordTrig=false}, { t"δ" }),
s({trig=";f", wordTrig=false}, { t"ϕ" }),
s({trig=";g", wordTrig=false}, { t"γ" }),
s({trig=";h", wordTrig=false}, { t"η" }),
s({trig=";k", wordTrig=false}, { t"κ" }),
s({trig=";l", wordTrig=false}, { t"λ" }),
s({trig=";z", wordTrig=false}, { t"ζ" }),
s({trig=";x", wordTrig=false}, { t"ξ" }),
s({trig=";c", wordTrig=false}, { t"χ" }),
s({trig=";v", wordTrig=false}, { t"\\sqrt{",i(1), t"}" }),
s({trig=";b", wordTrig=false}, { t"β" }),
s({trig=";n", wordTrig=false}, { t"ν" }),
s({trig=";m", wordTrig=false}, { t"μ" }),

s({trig=";Q", wordTrig=false}, { t"Θ", }),
s({trig=";W", wordTrig=false}, { t"Ω", }),
s({trig=";Y", wordTrig=false}, { t"Ψ", }),
s({trig=";U", wordTrig=false}, { t"Υ", }),
s({trig=";P", wordTrig=false}, { t"Π", }), -- product doesn't have a unicode equivalent
s({trig=";A", wordTrig=false}, { t"Å", }), -- using Ångström for A since greek A is identical to latin A. \AA is also Ångström in text mode.
s({trig=";S", wordTrig=false}, { t"Σ", }), -- note: not the sum sigma. Sigma symbol is sometimes used as a variable.
s({trig=";D", wordTrig=false}, { t"Δ", }),
s({trig=";F", wordTrig=false}, { t"Φ", }),
s({trig=";G", wordTrig=false}, { t"Γ", }),
s({trig=";L", wordTrig=false}, { t"Λ", }),
s({trig=";X", wordTrig=false}, { t"Ξ", }),

s({trig="\\nabla", wordTrig=false}, t"∇", optm),

s({trig="*", dscr="cdot", wordTrig=false}, t"⋅", optm),
s({trig="xx", dscr="cross", wordTrig=false}, t"×", optm),
s({trig="\\pm", dscr="plus/minus", wordTrig=false}, t"±", optm),
s({trig="\\mp", dscr="plus/minus", wordTrig=false}, t"∓", optm),
s({trig="+-", dscr="plus/minus", wordTrig=false}, t"±", optm),
s({trig="-+", dscr="plus/minus", wordTrig=false}, t"∓", optm),

s({trig="\\dots", descr="ellipses", wordTrig=false}, t"…", optm),
-- \dots, not \ldots is technically the equivalent of …
s({trig="\\ldots", descr="ellipses", wordTrig=false}, t"…", optm),
s({trig="\\cdots", descr="center dots", wordTrig=false}, t"⋯", optm),
s({trig="\\vdots", descr="vertical dots", wordTrig=false}, t"⋮", optm),
-- lower priority so c... and v... takes precedence
s({trig="...", descr="ellipses", wordTrig=false, priority=100}, t"…", optm),
s({trig="c...", descr="center dots", wordTrig=true}, t"⋯", optm),
s({trig="v...", descr="vertical dots", wordTrig=true}, t"⋮", optm),

s({trig="\\circ", dscr="Circle/degrees"}, t"∘", optm),
s({trig="deg", dscr="degrees", snippetType="snippet"}, t"^∘\\text{C}", optm),
s({trig="degc", dscr="degrees", wordTrig=false}, t"^∘\\text{C}", optm),

-- √ unfortunately doesn't work like \sqrt so we have to replace it.
-- The reason is there is redundancy, since \surd is a separate command for a √ 
-- by itself.
-- Still a useful trigger, since it is bound to alt+v.
s( {trig="√", dscr="invalid square root"}, {t"\\sqrt{", d(1, get_visual), t"}"}, optm),

-- lower priority so we can replace \oo first
s({trig="oo", descr="infinity", priority=100}, t"∞", optm),
s({trig="\\oo", descr="infinity", wordTrig=false}, t"∞", optm),
s({trig="\\infty", descr="infinity"}, t"∞", optm),

s({trig="([^\\])lim", descr="limit", regTrig=true},
-- {t"\\lim_{", i(1, "n"), t"\\to\\infty}"},
-- space after \to since a vimtex "to" snippet completes on ls.expand_or_jump
-- other solutions than leaving the space would be removing the default "to" snippet somehow,
-- or having separate keys for expand and jump
{re(1), t"\\lim_{", i(1, "n"), t"→", i(2, "∞"), t"}"}, optm),

s({trig="sum", descr="sum"},
{t"∑_{", i(1, "n=1"), t"}^{", i(2,"∞"), t"}"}, optm),

s({trig="part", dscr="d/dx", snippetType="snippet"},
{t"\\frac{∂ ", i(1,"V"), t"}{∂ ", i(2,"x"), t"}"}, optm),
s({trig="partial"}, t"∂", optm),

-- \to and \rightarrow are synonymous
s({trig="->", descr="to/right arrow"}, t"→", optm),
s({trig="<-", descr="left arrow"}, t"←", optm),
-- Currently using \implies instead which is longer with more whitespace
-- s({trig="=>", descr="Double right arrow"}, t"⇒", optm),
-- Currently using \impliedby instead which is longer with more whitespace
-- s({trig="=<", descr="Double left arrow"}, t"⇐", optm),
s({trig="!>", descr="maps to"}, t"↦", optm),
s({trig="\\mapsto", descr="maps to"}, t"↦", optm),

s({trig="!=", descr="not equal", wordTrig=false}, t"≠", optm),
-- can't be used when we have '==' -> '&= ' which is more important
-- s({trig="===", dscr="equivalent"}, t"≡", optm),
s({trig="\\?equiv", dscr="equivalent", regTrig=true}, t"≡", optm),
s({trig="<=", descr="less than or equal", wordTrig=false}, t"≤"),
s({trig=">=", descr="greater than or equal", wordTrig=false}, t"≥"),
s({trig="\\propto", descr="proportional to"}, t"∝", optm),
-- lower priority than =~ (congruent) and e.g. p~ -> ~ on top of p (\tilde{p}).
s({trig="~", descr="tilde", priority=100}, t"\\sim", optm),
s({trig="\\sim~", descr="approximately"}, t"≈", optm),
s({trig="\\sim=", descr="asymptotically equal to"}, t"\\simeq", optm),
s({trig="=~", descr="congruent"}, t"\\cong", optm),

-- logic
s({trig="EE", descr="E (set)", wordTrig=false}, t"∃", optm),
-- \AA is Ångström
s({trig="([^\\])AA", descr="A (set)", regTrig=true, wordTrig=false},
{re(1), t"∀"}, optm),
s({trig="&&", descr="logical and"}, t"∧", optm),
s({trig="||", descr="logical and"}, t"∨", optm),

-- set theory
s({trig="\\in ", descr="in"}, t"∈ ", optm),
s({trig="inn", descr="in"}, t"∈", optm),
s({trig="([^\\])notin", regTrig=true, wordTrig=false, descr="not in"}, {re(1), t"∉"}, optm),
s({trig="\\neg", descr="negate"}, t"¬", optm),
s({trig="OO", descr="empty set", wordTrig=false}, t"∅", optm),
s({trig="\\?empty", descr="empty set", regTrig=true}, t"∅", optm),
s({trig="cc", descr="subset", wordTrig=true}, t"⊂", optm),
s({trig="⊂=", descr="subset or equal"}, t"⊆", optm),
s({trig="c=", descr="subset or equal"}, t"⊆", optm),
s({trig="\\cup", descr="union"}, t"∪", optm),
s({trig="\\cap", descr="intersection"}, t"∩", optm),
s({trig="union", descr="union"}, t"∪", optm),
s({trig="intersect", descr="intersection"}, t"∩", optm),
s({trig="UU", descr="union", wordTrig=false}, t"∪", optm),
s({trig="NN", descr="intersection", wordTrig=false}, t"∩", optm),

s({trig="dint", descr="integral", priority=300},
{t"∫_{", i(1, "-∞"), t"}^{", i(2,"∞"), t"}", d(3, get_visual)}, optm),

s({trig="\\int", descr="integral"}, t"∫", optm),


}
