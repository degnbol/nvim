-- sources:
-- https://castel.dev/post/lecture-notes-1/
-- https://www.ejmastnak.com/tutorials/vim-latex/luasnip
-- https://github.com/gillescastel/latex-snippets/blob/master/tex.snippets

local lsu = require "utils/luasnip"
local vtu = require "utils/vimtex"

local get_visual = lsu.get_visual
local virt = lsu.virt
local re = lsu.re
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
% Use uppercase instead of lowercase letters for sub-figure numbering
% (Per default \thesubfigure is defined as \alph{subfigure}, i.e. lowercase letters)
\renewcommand\thesubfigure{\Alph{subfigure}}
\renewcommand\thesubtable{\Alph{subtable}}
% singlelinecheck=on means if only single line caption, then ignore raggedright (center).
\captionsetup[subfigure]{singlelinecheck=off,justification=raggedright}
% remove () around subfig numbering
\captionsetup[sub]{format=plain,position=top,font+=Large,labelformat=simple}
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
{condition=conds.line_begin+conds.line_end}),


}
