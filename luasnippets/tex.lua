local in_math = function()
  -- requires the VimTeX plugin
  return vim.fn['vimtex#syntax#in_mathzone']() == 1
end
local in_text = function() return not in_math() end
local in_comment = function()
  return vim.fn['vimtex#syntax#in_comment']() == 1
end
local in_env = function(name)
    local is_inside = vim.fn['vimtex#env#is_inside'](name)
    return (is_inside[1] > 0 and is_inside[2] > 0)
end
local in_itemize = function() return in_env('itemize') end

-- Summary: When `SELECT_RAW` is populated with a visual selection, the function
-- returns an insert node whose initial text is set to the visual selection.
-- When `SELECT_RAW` is empty, the function simply returns an empty insert node.
local get_visual = function(args, parent)
  if (#parent.snippet.env.SELECT_RAW > 0) then
    return sn(nil, i(1, parent.snippet.env.SELECT_RAW))
  else  -- If SELECT_RAW is empty, return a blank insert node
    return sn(nil, i(1))
  end
end

local regGroup = function(_, snip, i) return snip.captures[i] end
local re = function(i) return f(regGroup, nil, {user_args={i}}) end

return {
s({trig=";q", snippetType="autosnippet", wordTrig=false}, { t("\\theta"), }),
s({trig=";w", snippetType="autosnippet", wordTrig=false}, { t("\\omega"), }),
s({trig=";e", snippetType="autosnippet", wordTrig=false}, { t("\\epsilon"), }),
s({trig=";r", snippetType="autosnippet", wordTrig=false}, { t("\\rho"), }),
s({trig=";t", snippetType="autosnippet", wordTrig=false}, { t("\\tau"), }),
s({trig=";y", snippetType="autosnippet", wordTrig=false}, { t("\\psi"), }),
s({trig=";u", snippetType="autosnippet", wordTrig=false}, { t("\\upsilon"), }),
s({trig=";i", snippetType="autosnippet", wordTrig=false}, { t("\\iota"), }),
s({trig=";o", snippetType="autosnippet", wordTrig=false}, { t("\\omicron"), }),
s({trig=";p", snippetType="autosnippet", wordTrig=false}, { t("\\pi"), }),
s({trig=";a", snippetType="autosnippet", wordTrig=false}, { t("\\alpha"), }),
s({trig=";s", snippetType="autosnippet", wordTrig=false}, { t("\\sigma"), }),
s({trig=";d", snippetType="autosnippet", wordTrig=false}, { t("\\delta"), }),
s({trig=";f", snippetType="autosnippet", wordTrig=false}, { t("\\phi"), }),
s({trig=";g", snippetType="autosnippet", wordTrig=false}, { t("\\gamma"), }),
s({trig=";h", snippetType="autosnippet", wordTrig=false}, { t("\\eta"), }),
s({trig=";k", snippetType="autosnippet", wordTrig=false}, { t("\\kappa"), }),
s({trig=";l", snippetType="autosnippet", wordTrig=false}, { t("\\lambda"), }),
s({trig=";z", snippetType="autosnippet", wordTrig=false}, { t("\\zeta"), }),
s({trig=";x", snippetType="autosnippet", wordTrig=false}, { t("\\xi"), }),
s({trig=";c", snippetType="autosnippet", wordTrig=false}, { t("\\chi"), }),
s({trig=";v", snippetType="autosnippet", wordTrig=false}, { t("\\sqrt{"),i(1), t("}") }),
s({trig=";b", snippetType="autosnippet", wordTrig=false}, { t("\\beta"), }),
s({trig=";n", snippetType="autosnippet", wordTrig=false}, { t("\\nu"), }),
s({trig=";m", snippetType="autosnippet", wordTrig=false}, { t("\\mu"), }),

s({trig="beg", snippetType="autosnippet"},
  fmta(
    [[
      \begin{<>}
          <>
      \end{<>}
    ]],
    {
      i(1),
      i(2),
      rep(1),  -- this node repeats insert node i(1)
    }
  ),
  { condition = conds.line_begin, }
),

s({trig="h1", dscr="Top-level section", snippetType="autosnippet"},
  fmta(
    [[\section{<>}]],
    { i(1) }
  ), 
  {condition = conds.line_begin}
),
s({trig="h2", dscr="Sub-section", snippetType="autosnippet"},
  fmta(
    [[\subsection{<>}]],
    { i(1) }
  ), 
  {condition = conds.line_begin}
),
s({trig="h3", dscr="Sub-sub-section", snippetType="autosnippet"},
  fmta(
    [[\subsubsection{<>}]],
    { i(1) }
  ), 
  {condition = conds.line_begin}
),


s({trig="href", dscr="The hyperref package's href{}{} command (for url links)"},
  fmta(
    [[\href{<>}{<>}]],
    {
      i(1, "url"),
      i(2, "display name"),
    }
  )
),


s({trig = "tii", dscr = "Expands 'tii' into LaTeX's textit{} command."},
  fmta("\\textit{<>}",
    {
      d(1, get_visual),
    }
  )
),
s({trig = "tbb", dscr = "Expands 'tii' into LaTeX's textit{} command."},
  fmta("\\textbf{<>}",
    {
      d(1, get_visual),
    }
  )
),


-- mm regex below is cooler
-- s(
--     {trig="mk", dscr="inline math", snippetType="autosnippet"},
--     {
--         t("$"), i(1), t("$")
--     }
-- ),

s({trig = "([^%a])mm", wordTrig = false, regTrig = true, snippetType="autosnippet"},
  fmta(
    "<>$<>$",
    {
      f( function(_, snip) return snip.captures[1] end ),
      d(1, get_visual),
    }
  )
),

s({trig="dm", snippetType="autosnippet"},
  fmta(
    [[
    \[
        <>
    \]

    ]],
    { i(1), }
)),

-- from https://castel.dev/post/lecture-notes-1/
-- and https://www.ejmastnak.com/tutorials/vim-latex/luasnip

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

s(
{trig="([%w\\]+)/", dscr="fraction after word", regTrig=true, snippetType="autosnippet"},
{t"\\frac{", re(1), t"}{", i(1), t"}"}
),

s(
{trig=".", dscr="cdot", snippetType="autosnippet"},
t("\\cdot"),
{condition=in_math}
),
s(
{trig="*", dscr="cdot", wordTrig=false, snippetType="autosnippet"},
t("\\cdot"),
{condition=in_math}
),

s({trig = '([^%a])ee', regTrig = true, wordTrig=false, snippetType="autosnippet"},
  fmta(
    "<>e^{<>}",
    { re(1), d(1, get_visual) }
  ),
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

-- lower priority so we can replace \oo first
s({trig="oo", descr="infinity", snippetType="autosnippet", priority=100},
t"\\infty",
{condition=in_math}),
s({trig="\\oo", descr="infinity", wordTrig=false, snippetType="autosnippet"},
t"\\infty",
{condition=in_math}),

s({trig="lim", descr="limit", snippetType="autosnippet"},
-- {t"\\lim_{", i(1, "n"), t"\\to\\infty}"},
-- space after \to since a vimtex "to" snippet completes on ls.expand_or_jump
-- other solutions than leaving the space would be removing the default "to" snippet somehow,
-- or having separate keys for expand and jump
{t"\\lim_{", i(1, "n"), t"\\to ", i(2, "\\infty"), t"}"},
{condition=in_math}),

s({trig="sum", descr="sum", snippetType="autosnippet"},
{t"\\sum_{", i(1, "n=1"), t"}^{", i(2,"\\infty"), t"}"},
{condition=in_math}),

s({trig="->", descr="to", snippetType="autosnippet"},
t"\\to",
{condition=in_math}),

s({trig="!>", descr="maps to", snippetType="autosnippet"},
t"\\mapsto",
{condition=in_math}),

s({trig="sr", descr="squared", wordTrig=false, snippetType="autosnippet"},
t"^2",
{condition=in_math}),

s({trig="cb", descr="cubed", wordTrig=false, snippetType="autosnippet"},
t"^3",
{condition=in_math}),

s({trig="td", descr="to the power", wordTrig=false, snippetType="autosnippet"},
{t"^{", i(1), t"}"},
{condition=in_math}),



}
