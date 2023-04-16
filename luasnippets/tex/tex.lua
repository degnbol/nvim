-- sources:
-- https://castel.dev/post/lecture-notes-1/
-- https://www.ejmastnak.com/tutorials/vim-latex/luasnip
-- https://github.com/gillescastel/latex-snippets/blob/master/tex.snippets

local lsu = require"luasnip_util"
local get_visual = lsu.get_visual

return {
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
\usepackage{unicode-math} % experimental support for unicode in math

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
% \usepackage[capitalise]{cleverref} % \cref which auto adds e.g. "Table " to \ref
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

s({trig="pac", dscr="package", snippetType="autosnippet"},
{t"\\usepackage[", i(1, "options"), t"]{", i(2, "package"), t"}"},
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

}
