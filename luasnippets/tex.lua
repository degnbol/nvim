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

s({trig=";Q", snippetType="autosnippet", wordTrig=false}, { t("\\Theta"), }),
s({trig=";W", snippetType="autosnippet", wordTrig=false}, { t("\\Omega"), }),
s({trig=";E", snippetType="autosnippet", wordTrig=false}, { t("\\Epsilon"), }),
s({trig=";R", snippetType="autosnippet", wordTrig=false}, { t("\\Rho"), }),
s({trig=";T", snippetType="autosnippet", wordTrig=false}, { t("\\Tau"), }),
s({trig=";Y", snippetType="autosnippet", wordTrig=false}, { t("\\Psi"), }),
s({trig=";U", snippetType="autosnippet", wordTrig=false}, { t("\\Upsilon"), }),
s({trig=";I", snippetType="autosnippet", wordTrig=false}, { t("\\Iota"), }),
s({trig=";O", snippetType="autosnippet", wordTrig=false}, { t("\\Omicron"), }),
s({trig=";P", snippetType="autosnippet", wordTrig=false}, { t("\\Pi"), }),
s({trig=";A", snippetType="autosnippet", wordTrig=false}, { t("\\Alpha"), }),
s({trig=";S", snippetType="autosnippet", wordTrig=false}, { t("\\Sigma"), }),
s({trig=";D", snippetType="autosnippet", wordTrig=false}, { t("\\Delta"), }),
s({trig=";F", snippetType="autosnippet", wordTrig=false}, { t("\\Phi"), }),
s({trig=";G", snippetType="autosnippet", wordTrig=false}, { t("\\Gamma"), }),
s({trig=";H", snippetType="autosnippet", wordTrig=false}, { t("\\Eta"), }),
s({trig=";K", snippetType="autosnippet", wordTrig=false}, { t("\\Kappa"), }),
s({trig=";L", snippetType="autosnippet", wordTrig=false}, { t("\\Lambda"), }),
s({trig=";Z", snippetType="autosnippet", wordTrig=false}, { t("\\Zeta"), }),
s({trig=";X", snippetType="autosnippet", wordTrig=false}, { t("\\Xi"), }),
s({trig=";C", snippetType="autosnippet", wordTrig=false}, { t("\\Chi"), }),
s({trig=";B", snippetType="autosnippet", wordTrig=false}, { t("\\Beta"), }),
s({trig=";N", snippetType="autosnippet", wordTrig=false}, { t("\\Nu"), }),
s({trig=";M", snippetType="autosnippet", wordTrig=false}, { t("\\Mu"), }),

-- TODO: maybe add toggling between different templates.
-- https://tex.stackexchange.com/questions/553/what-packages-do-people-load-by-default-in-latex
s("template",
-- < and > chars are escaped in fmta call by typing << and >>
fmta([[
% !TEX program = LuaLaTeX
\documentclass[a4paper,10pt]{article}
\usepackage[margin=2cm, top=0.5in]{geometry}
\usepackage{xspace} % \xspace at end of newcommand allows "\CUSTOM " instead of "\CUSTOM\ "

\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage[english]{babel}
\usepackage{textcomp}
\usepackage{lmodern} % improved default font

% Improve default latex packages.
% allowing (2%) letter stretch
\usepackage{microtype}
% \raggedright ->> \RaggedRight, \flushleft env ->> \FlushLeft, \center ->> \Center
% https://www.overleaf.com/learn/latex/Text_alignment
\usepackage{ragged2e}
% \usepackage{flushend} % enable for two columns, to get equal lengths on last page.

% \begin{itemize}[label={--}] instead of repeating /item[--]
\usepackage{enumitem}
\usepackage{fancyhdr} % header/footer control
% paragraph space instead of indent
\usepackage[parfill]{parskip} 

% math
\usepackage{mathtools} % loads amsmath, plus e.g. \coloneqq,\mathclap,\substack
\usepackage{amssymb}
\usepackage{siunitx} % provides \SI and column type S (align decimal)

\usepackage{graphicx}
% Subfigures. "skip" == spacing between subfigures.
\usepackage[skip=0pt]{subcaption}
\usepackage{float}
% makes figures and tables stay in their section
\usepackage[section]{placeins}

% Tables.
% https://tex.stackexchange.com/questions/12672/which-tabular-packages-do-which-tasks-and-which-packages-conflict
\usepackage{array} % flexible column formatting
\usepackage{tabularx} % Column type X for width filling.
\usepackage{tabulary} % Column type L,C,R,J for balanced width versions of l,c,r,j.
\usepackage{booktabs} % Better vertical spacing. Midrule etc with varying thickness instead of \hline.
% creates missing tabularx column type named "R" to specify right adjustment
\newcolumntype{R}{>>{\raggedleft\arraybackslash}X}
\usepackage{multicol, multirow} % Also see pbox.
% \usepackage{longtable} % Multipage table. Might not support column X. Alts: xltabular, ltxtable

\usepackage{xcolor} % define colors
\usepackage{hyperref}
\usepackage[capitalise]{cleverref} % \cref which auto adds e.g. "Table " to \ref
% Auto define acronyms on first use. https://www.overleaf.com/learn/latex/Glossaries
\usepackage[acronym]{glossaries-extra}
\setabbreviationstyle[acronym]{long-short}
% bibliography
\usepackage{biblatex}
% appendix
\usepackage[toc,page]{appendix}
% Code blocks.
% https://www.overleaf.com/learn/latex/Code_listing
% https://www.overleaf.com/learn/latex/Code_Highlighting_with_minted
% \usepackage{listings}
% \usepackage{minted} % supports julia

\begin{document}
	<>
\end{document}
]], {i(1),}),
{condition=conds.line_begin}),


s({trig="beg", snippetType="autosnippet"},
  fmta(
    [[
      \begin{<>}
          <>
      \end{<>}
    ]],
    -- rep node repeats insert node i(1)
    { i(1), i(2), rep(1), }
  ),
  { condition = conds.line_begin }
),

s({trig="pac", dscr="package", snippetType="autosnippet"},
{t"\\usepackage[", i(1, "options"), t"]{", i(2, "package"), t"}"},
{condition=conds.line_begin}),

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
    { i(1, "url"), i(2, "display name"), }
  )
),


s({trig = "tii", dscr = "Expands 'tii' into LaTeX's textit{} command."},
  fmta("\\textit{<>}",
    { d(1, get_visual), }
  )
),
s({trig="tbb", dscr="Expands 'tii' into LaTeX's textit{} command."},
  fmta("\\textbf{<>}",
    { d(1, get_visual), }
  )
),

s({trig="enum", dscr="enumerate", snippetType="autosnippet"},
  fmta(
[[
\begin{enumerate}
	\item <>
\end{enumerate}
]],
{ i(0), }),
{condition=conds.line_begin}
),

-- TODO: add toggle for bullet or dash
s({trig="item", dscr="itemize", snippetType="autosnippet"},
  fmta(
[[
\begin{itemize}
	\item <>
\end{itemize}
]],
{ i(0), }),
{condition=conds.line_begin}
),

s({trig="desc", dscr="description", snippetType="autosnippet"},
  fmta(
[[
\begin{description}
	\item[<>] <>
\end{description}
]],
{ i(1), i(0), }),
{condition=conds.line_begin}
),

s({trig="ali", dscr="align", snippetType="autosnippet"},
  fmta(
[[
\begin{align*}
	<>
\end{align*}
]],
    { d(1, get_visual), }
  )
),

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

-- sources:
-- https://castel.dev/post/lecture-notes-1/
-- https://www.ejmastnak.com/tutorials/vim-latex/luasnip
-- https://github.com/gillescastel/latex-snippets/blob/master/tex.snippets

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
{t"\\frac{", re(1), t"}{", i(1), t"}"},
{condition=in_math}
),

s(
{trig="*", dscr="cdot", wordTrig=false, snippetType="autosnippet"},
t("\\cdot"),
{condition=in_math}
),
s({trig="...", descr="ellipses", wordTrig=false, snippetType="autosnippet"},
t"\\ldots",
{condition=in_math}),

s(
{trig="xx", dscr="cross", wordTrig=false, snippetType="autosnippet"},
t("\\times"),
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

-- lower priority so we can replace \oo first
s({trig="oo", descr="infinity", snippetType="autosnippet", priority=100},
t"\\infty",
{condition=in_math}),
s({trig="\\oo", descr="infinity", wordTrig=false, snippetType="autosnippet"},
t"\\infty",
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

s({trig="part", dscr="d/dx"},
{t"\\frac{\\partial ", i(1,"V"), t"}{\\partial ", i(2,"x"), t"}"},
{condition=in_math}),

s({trig="->", descr="to", snippetType="autosnippet"},
t"\\to",
{condition=in_math}),

s({trig="!>", descr="maps to", snippetType="autosnippet"},
t"\\mapsto",
{condition=in_math}),

s({trig="sq", dscr="sqrt", wordTrig=false, snippetType="autosnippet"},
{t"\\sqrt{", d(1, get_visual), t"}"},
{condition=in_math}
),

s({trig="sr", descr="squared", wordTrig=false, snippetType="autosnippet"},
t"^2",
{condition=in_math}),

s({trig="cb", descr="cubed", wordTrig=false, snippetType="autosnippet"},
t"^3",
{condition=in_math}),

s({trig="td", descr="to the power", wordTrig=false, snippetType="autosnippet"},
{t"^{", i(1), t"}"},
{condition=in_math}),
s({trig="tD", descr="to the (power)", wordTrig=false, snippetType="autosnippet"},
{t"^{(", i(1), t")}"},
{condition=in_math}),

s({trig="__", descr="subscript", wordTrig=false, snippetType="autosnippet"},
{t"_{", i(1), t"}"},
{condition=in_math}),

s({trig="<=", descr="less than or equal", wordTrig=false, snippetType="autosnippet"},
{t"\\le"}),
s({trig=">=", descr="greater than or equal", wordTrig=false, snippetType="autosnippet"},
{t"\\ge"}),

s({trig="=>", descr="implies", wordTrig=false, snippetType="autosnippet"},
{t"\\implies"},
{condition=in_math}),

s({trig="==", descr="equal", wordTrig=false, snippetType="autosnippet"},
{t"&= ", i(1)},
{condition=in_math}),
s({trig="!=", descr="not equal", wordTrig=false, snippetType="autosnippet"},
t"\\neq",
{condition=in_math}),

s({trig="notin", descr="not in", wordTrig=false, snippetType="autosnippet"},
t"\\not\\in",
{condition=in_math}),

s({trig=[[\\\]], descr="set minus", wordTrig=false, snippetType="autosnippet"},
t"\\setminus",
{condition=in_math}),

s({trig="EE", descr="E (set)", wordTrig=false, snippetType="autosnippet"},
t"\\exists",
{condition=in_math}),
s({trig="AA", descr="A (set)", wordTrig=false, snippetType="autosnippet"},
t"\\forall",
{condition=in_math}),

s({trig="UU", descr="union", wordTrig=false, snippetType="autosnippet"},
t"\\cup",
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

s({trig="dint", descr="integral", snippetType="autosnippet", priority=300},
{t"\\int_{", i(1, "-\\infty"), t"}^{", i(2,"\\infty"), t"}", d(3, get_visual)},
{condition=in_math}),

}
