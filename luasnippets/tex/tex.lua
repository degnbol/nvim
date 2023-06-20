-- sources:
-- https://castel.dev/post/lecture-notes-1/
-- https://www.ejmastnak.com/tutorials/vim-latex/luasnip
-- https://github.com/gillescastel/latex-snippets/blob/master/tex.snippets

local lsu = require"luasnip_util"
local get_visual = lsu.get_visual
local virt = lsu.virt
local re = lsu.re
local vtu = require"vimtex_util"
local in_text = vtu.in_text
local in_itemize = vtu.in_itemize
local in_description = vtu.in_description

return {
-- TODO: maybe add toggling between different templates.
-- E.g. for PdfLaTeX relevant to journal old-fashioned requirements maybe use package textgreek
-- https://tex.stackexchange.com/questions/553/what-packages-do-people-load-by-default-in-latex
s("template",
-- < and > chars are escaped in fmta call by typing << and >>
-- \usepackage[utf8]{inputenc} is no longer required since 2018 https://www.overleaf.com/learn/latex/Greek
-- \usepackage[T1]{fontenc} specifies output encoding and should also no longer 
-- be needed, especially with modern lualatex or xelatex.
-- \usepackage{textcomp} also seems redundant now.
fmta([[
% !TEX program = LuaLaTeX
\documentclass[a4paper,10pt]{article}
\usepackage[margin=2cm, top=0.5in]{geometry}
\usepackage{xspace} % \xspace at end of newcommand allows "\CUSTOM " instead of "\CUSTOM\ "

\usepackage{csquotes}
\usepackage{fontspec}
% don't \usepackage{lmodern}, see latex/fonts.tex for details.
\setmainfont{New Computer Modern 10}

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

s({trig="pdflatex", dscr="Use PdfLaTeX", condition=conds.line_begin},
t{"% !TEX program = PdfLaTeX", ""}),

s({trig="pac", dscr="package", snippetType="autosnippet"},
{t"\\usepackage", c(2, {t"", {t"[", i(1,"options"), t"]"}}), t"{", i(1, "package"), t"}"},
{condition=conds.line_begin}),

s({trig="beg", snippetType="autosnippet"},
fmta([[\begin{<>}
	<>
\end{<>}

]],
-- rep node repeats insert node i(1)
{ i(1), i(2), rep(1) }),
{ condition = conds.line_begin }
),

s({trig="h1", dscr="Top-level section", snippetType="autosnippet"},
  fmta([[\section{<>}]], { i(1) }), 
  {condition = conds.line_begin}
),
s({trig="h2", dscr="Sub-section", snippetType="autosnippet"},
  fmta([[\subsection{<>}]], { i(1) }), 
  {condition = conds.line_begin}
),
s({trig="h3", dscr="Sub-sub-section", snippetType="autosnippet"},
  fmta([[\subsubsection{<>}]], { i(1) }), 
  {condition = conds.line_begin}
),

-- single _ should never be found in regular text.
-- TODO: narrow the condition by considering things like verbatim text blocks.
s({trig="[^_]_", dscr="single _ outside math", condition=in_text, regTrig=true, wordTrig=false, snippetType='autosnippet'},
{t"\\_"}),

s({trig="\\?__", descr="subscript",   condition=in_text, regTrig=true, wordTrig=false, snippetType="autosnippet"},
{t"\\textsubscript{", i(1), t"}"}),
s({trig="^^", descr="superscript", condition=in_text, snippetType="autosnippet", wordTrig=false},
{t"\\textsuperscript{", i(1), t"}"}),


s({trig="href", dscr="The hyperref package's href{}{} command (for url links)"},
  fmta(
    [[\href{<>}{<>}]],
    { i(1, "url"), i(2, "display name"), }
  )
),

s({trig = "tii", dscr = "Italic", wordTrig=false, snippetType="autosnippet"},
  fmta("\\textit{<>}",
    { d(1, get_visual) }
  )
),
s({trig="tbb", dscr="Bold", wordTrig=false, snippetType="autosnippet"},
  fmta("\\textbf{<>}",
    { d(1, get_visual) }
  )
),
s({trig="ttt", dscr="typewriter", wordTrig=false, snippetType="autosnippet"},
  fmta("\\texttt{<>}",
    { d(1, get_visual) }
  )
),

s({trig="enum", dscr="enumerate", snippetType="autosnippet"},
fmta([[\begin{enumerate}
	\item <>
\end{enumerate}
]], { i(1) }),
{condition=conds.line_begin}
),

s({trig="item", dscr="itemize", snippetType="autosnippet"},
fmta([[\begin{itemize}<>
	\item <>
\end{itemize}

]], { c(1, {t("", virt("^l -> dash instead of bullet")), t"[label={--}]"}), i(2) }),
{condition=conds.line_begin}),

s({trig="desc", dscr="description", snippetType="autosnippet"},
fmta([[\begin{description}
	\item[<>] <>
\end{description}

]], { i(1), i(2) }), {condition=conds.line_begin}),

-- lower priority so the item snippets right below gets called
s({trig="- ", dscr="itemize dashed", snippetType="autosnippet", priority=100},
fmta([[\begin{itemize}[label={--}]
	\item <>
\end{itemize}

]], i(1)), {condition=conds.line_begin}),
-- in itemize
s({trig="- ", dscr="item", snippetType="autosnippet"},
t"\\item ", {condition=conds.line_begin and in_itemize}),
-- in description
s({trig="- ", dscr="item", snippetType="autosnippet"},
{t"\\item[", i(1), t"]"}, {condition=conds.line_begin and in_description}),


s({trig="mm", dscr="inline math", snippetType="autosnippet"},
{t"$", d(1, get_visual), t"$" },
{condition=conds.line_begin}),
s({trig="([^%a])mm", dscr="inline math", wordTrig=false, regTrig=true, snippetType="autosnippet"},
{re(1), t"$", d(1, get_visual), t"$" }),

s({trig="dm", snippetType="autosnippet"},
fmta([[\[
	<>
\]

]], i(1))),

s({trig="ali", dscr="align", snippetType="autosnippet"},
fmta([[\begin{align*}
	<>
\end{align*}
]], d(1, get_visual))),

s({trig="(%a)-", dscr="p-value, n-dimensional, ...", regTrig=true, snippetType="autosnippet"},
{t"$", re(1), t"$-"}),

s({trig="fig", dscr="fig", condition=conds.line_begin, snippetType="autosnippet"},
fmta([[\begin{figure}[ht]
	\centering
	\includegraphics[width=0.95\textwidth]{figures/<>}
	\caption{<>}
	\label{fig:<>}
\end{figure}

]], {i(1, "FILENAME"), i(2, "\\textbf{Title.} Caption."), i(3)})),

s({trig="subfig", dscr="subfig", condition=conds.line_begin, snippetType="autosnippet"},
fmta([[\begin{figure}[ht]
	\centering
	\begin{subfigure}[t]{0.49\textwidth}
		\centering
		\caption{}
		\includegraphics[width=0.95\textwidth]{figures/<>}
	\end{subfigure}
	\hfill
	\begin{subfigure}[t]{0.49\textwidth}
		\centering
		\caption{}
		\includegraphics[width=0.95\textwidth]{figures/<>}
	\end{subfigure}
	\caption{<>}
	\label{fig:<>}
\end{figure}

]], {i(1), i(2), i(3, "\\textbf{Title.} Caption."), i(4)})),


s({trig="tabx", dscr="tabularx", condition=conds.line_begin, snippetType="autosnippet"},
-- @{} suppresses space between columns. @{.} would use "." as column separator.
fmta([[\begin{table}[ht]
\caption{<>}
\begin{tabularx}{\textwidth}{@{}<>@{}}
	\toprule
	<> \\
	\midrule
	<> \\
	\bottomrule
\end{tabularx}
\end{table}

]], {
    i(4, [[\textbf{Title.} Caption.]]),
    i(1, [[lcX]]),
    i(2, [[left & right & filling]]),
    i(3, [[left & right & filling]])
})),

-- \cmidrule{2-4} makes a thin line covering column 2 to 4.
s({trig="\\cmidrule", dscr="rule spanning subset of columns", snippetType="autosnippet"},
fmta("\\cmidrule{<>-<>}", {i(1, "FROM"), i(2, "TO")})),

s({trig="\\multicolumn", dscr="cell spanning multiple columns", snippetType="autosnippet"},
fmta("\\multicolumn{<>}{<>}{<>}", {i(1, "3"), i(2, "c"), i(3, "TEXT")})),

}
