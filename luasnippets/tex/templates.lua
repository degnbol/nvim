-- sources:
-- https://castel.dev/post/lecture-notes-1/
-- https://www.ejmastnak.com/tutorials/vim-latex/luasnip
-- https://github.com/gillescastel/latex-snippets/blob/master/tex.snippets
local util = require "utils/init"
local lsu = require "utils/luasnip"
local vtu = require "utils/vimtex"

local get_visual = lsu.get_visual
local virt = lsu.virt
local re = lsu.re
local in_text = vtu.in_text
local in_itemize = vtu.in_itemize
local in_description = vtu.in_description

local rtp = vim.opt.runtimepath:get()[1]

--- Functionnode function to put text from file(s) at a node location
---:h luasnip-functionnode
---@user_args table fnames
---@return text contents from file
local function putfile(args, parent, user_args)
    local texts = {}
    for _, fname in ipairs(user_args) do
        table.insert(texts, util.readtext(rtp .. "/luasnippets/tex/templates/" .. fname .. ".tex"))
    end
    local text = table.concat(texts, '\n')
    return vim.split(text, '\n')
end
---@return functionnode
local function putfilenode(fnames)
    return f(putfile, {}, {user_args={fnames}})
end


return {
-- TODO: maybe add toggling between different templates.
-- E.g. for PdfLaTeX relevant to journal old-fashioned requirements maybe use package textgreek
-- https://tex.stackexchange.com/questions/553/what-packages-do-people-load-by-default-in-latex
s("template",
-- NOTE: < and > chars are escaped in fmta call by typing << and >>
-- \usepackage[utf8]{inputenc} is no longer required since 2018 https://www.overleaf.com/learn/latex/Greek
-- \usepackage[T1]{fontenc} specifies output encoding and should also no longer 
-- be needed, especially with modern lualatex or xelatex.
-- \usepackage{textcomp} also seems redundant now.
fmta([[
% !TEX program = LuaLaTeX
\documentclass[a4paper,10pt]{article}
\usepackage[margin=2cm, top=0.5in]{geometry}

<>

\title{TITLE}
\author{Christian Degnbol Madsen}

\begin{document}
\maketitle
% \tableofcontents

<>

% \clearpage

<>

% \printbibliography

\end{document}
]], {
    putfilenode {
        "core",
        "fonts",
        "math",
        "subfig",
        "tables",
        "cleveref",
        "acro", -- or "glossaries",
        "verbatim",
        "biblatex",
    },
    i(2),
    putfilenode {
        "acro-doc",
    },
}),
{show_condition=conds.line_end}),

s("beamer",
-- < and > chars are escaped in fmta call by typing << and >>
-- \usepackage[utf8]{inputenc} is no longer required since 2018 https://www.overleaf.com/learn/latex/Greek
-- \usepackage[T1]{fontenc} specifies output encoding and should also no longer 
-- be needed, especially with modern lualatex or xelatex.
-- \usepackage{textcomp} also seems redundant now.
fmta([[
% !TEX program = LuaLaTeX
\documentclass[8pt,aspectratio=169]{beamer}
\usetheme[progressbar=frametitle]{metropolis}
\beamertemplatenavigationsymbolsempty
\usepackage{appendixnumberbeamer} % reset frame counter when calling \appendix

<>

% set larger default separation between list items with pac enumitem in core.tex
\setlist[1]{itemsep=1em} % only increase sep for highest level (1)

\begin{document}

\begin{frame}
	\thispagestyle{empty}\setcounter{framenumber}{0}
	\vspace{3cm}
	\resizebox{\textwidth}{!}{\huge <>}

	\resizebox{\textwidth}{!}{
	\textcolor{darkgray}{\small
		\textbf{Christian D. Madsen}\textsuperscript{1,2}<>
	}}

	\vspace{2cm}
	\textcolor{darkgray}{\tiny
		\textsuperscript{1}School of Mathematics and Statistics, University of Melbourne \\
		\textsuperscript{2}Melbourne Integrative Genomics, University of Melbourne \\
	}

\end{frame}

\begin{frame}
    <>
\end{frame}

\end{document}
]], {
    putfilenode {
        "core",
        "fonts",
        "math",
        "subfig",
        "tables",
        "cleveref",
        "verbatim",
    },
    i(1, "TITLE"),
    i(2),
    i(3),
}),
{show_condition=conds.line_end}),


}
