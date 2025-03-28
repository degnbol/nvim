-- sources:
-- https://castel.dev/post/lecture-notes-1/
-- https://www.ejmastnak.com/tutorials/vim-latex/luasnip
-- https://github.com/gillescastel/latex-snippets/blob/master/tex.snippets

local ls = require "luasnip"
local lsu = require "utils/luasnip"
local vtu = require "utils/vimtex"
local sn = ls.sn
local s = ls.s
local t = ls.t
local i = ls.i
local d = ls.d
local c = ls.c
local get_visual = lsu.get_visual
local virt = lsu.virt
local re = lsu.re
local in_text = vtu.in_text
local cond_itemize = vtu.cond_itemize
local cond_description = vtu.cond_description
local fmta = require("luasnip.extras.fmt").fmta

---For use with f(acro, {})
---@param _ table args, unused
---@param snip table
---@return string replacement
local function acro(_, snip)
    local text = snip.captures[1]
    local cmd
    -- dsb -> \ac
    -- Dsb -> \Ac
    if text:sub(1,1):match("%u") and text:match("%l") then
        cmd='\\Ac'
        text=text:sub(1,1):lower() .. text:sub(2)
    else
        cmd='\\ac'
    end
    -- dsbs -> \acp (if you need s at end, make without postfix, or postfix, then add s)
    if text:sub(-1):match("s") then
        cmd=cmd..'p'
        text=text:sub(1,-2)
    end
    return cmd .. '{' .. text .. '}'
end

return {
--

s({trig="\\newcommand", dscr="New command", show_condition=lsu.line_end},
{t"\\newcommand{\\", i(1), t"}{", i(2), t"}"}),


s({trig="%!", dscr="TEX directive", condition=conds.line_begin, snippetType="autosnippet"},
{t"% !TEX ", c(1, {t"TS-program", t"encoding", t"root", t"spellcheck"}), t" = ", i(2, "VALUE")}),

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

-- single _ should never be found in regular text, however it is fine in e.g. \includegraphics filepaths so we need to condition better.
-- TODO: narrow the condition by considering things like comments, verbatim text blocks, and filepaths.
-- Similarly % isn't allowed (since it's a comment), and it sometimes appears when pasting a link into \href{url}{display} so it could be nice to either highlight that as error or auto-replace
-- s({trig="([^_])_", dscr="single _ outside math", condition=in_text, regTrig=true, wordTrig=false, snippetType='autosnippet'},
-- {re(1), t"\\_"}),

s({trig="\\?__", descr="subscript",   condition=in_text, regTrig=true, wordTrig=false, snippetType="autosnippet"},
{t"\\textsubscript{", i(1), t"}"}),
s({trig="^^", descr="superscript", condition=in_text, snippetType="autosnippet", wordTrig=false},
{t"\\textsuperscript{", i(1), t"}"}),

s({trig="sansserif", dscr="Sans serif"},
{t"{\\sffamily ", i(1), t"}"}),

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

]], { c(1, {t("", virt("^/ -> dash instead of bullet")), t"[label={--}]"}), i(2) }),
{condition=conds.line_begin}),

s({trig="description", dscr="description"},
fmta([[\begin{description}
	\item[<>] <>
\end{description}

]], { i(1), i(2) }), {condition=conds.line_begin}),

-- lower priority so the item snippets right below gets called
s({
    trig="^- ",
    regTrig=true,
    dscr="itemize dashed",
    snippetType="autosnippet",
    priority=100,
    condition=conds.line_begin*lsu.line_end
},
fmta([[\begin{itemize}[label={--}]
	\item <>
\end{itemize}

]], i(1))),

-- in itemize
s({
    trig="- ",
    dscr="item",
    snippetType="autosnippet",
    condition=conds.line_begin*cond_itemize,
},
t"\\item "),
-- in description
s({trig="- ", dscr="item", snippetType="autosnippet"},
{t"\\item[", i(1), t"] "}, {condition=conds.line_begin*cond_description}),


-- Can't use this since we have e.g. mm10 mouse reference genome.
-- s({trig="mm", dscr="inline math", snippetType="autosnippet"},
-- {t"$", d(1, get_visual), t"$" },
-- {condition=conds.line_begin}),
-- Making sure to not allow alphanumeric prefix, since it will ruin words that contain mm,
-- and not numbers either since it will ruin e.g. 4mm when writing units for tikz.
-- s({trig="([^%w])mm", dscr="inline math", wordTrig=false, regTrig=true, snippetType="autosnippet"},
-- {re(1), t"$", d(1, get_visual), t"$" }),

s({trig="dm", snippetType="autosnippet", condition=conds.line_begin*lsu.line_end},
fmta([[\[
	<>
\]

]], i(1))),

s({trig="align", dscr="align", snippetType="autosnippet", condition=conds.line_begin*lsu.line_end},
fmta([[\begin{align*}
	<>
\end{align*}
]], d(1, get_visual))),

-- causing problems for e.g. writing fig. A-B
-- s({trig="(%a)-", dscr="p-value, n-dimensional, ...", regTrig=true, snippetType="autosnippet"},
-- {t"$", re(1), t"$-"}),

s({trig="\\includegraphics", dscr="includegraphics",
        -- condition=conds.line_begin*lsu.line_end
    },
{t"\\includegraphics[width=\\textwidth,height=\\textheight]{./figures", i(1), t"}"}),

s({trig="fig", dscr="fig", condition=conds.line_begin*lsu.line_end, snippetType="autosnippet"},
fmta([[\begin{figure}[htb]
	\centering
	\includegraphics[width=\textwidth,height=\textheight]{./figures<>}
	\caption{
	    \textbf{<>.}
	    <>.
	}
	\label{fig:<>}
\end{figure}

]], {i(1, "FILENAME"), i(2, "TITLE"), i(3, "CAPTION"), i(3)})),

s({trig="subfig", dscr="subfig", condition=conds.line_begin, snippetType="autosnippet"},
fmta([[\begin{figure}[htb]
	\centering
	\begin{subfigure}[t]{0.49\textwidth}
		\centering
		\caption{}
		\includegraphics[width=\textwidth,height=\textheight]{./figures<>}
		\label{fig:}
	\end{subfigure}
	\hfill
	\begin{subfigure}[t]{0.49\textwidth}
		\centering
		\caption{}
		\includegraphics[width=\textwidth,height=\textheight]{./figures<>}
		\label{fig:}
	\end{subfigure}
	\caption{
	    \textbf{<>.}
	    \subfig{A}~.
	    \subfig{B}~refers back to \subfigother{A}.
	}
	\label{fig:<>}
\end{figure}

]], {i(1), i(2), i(3, "TITLE"), i(4)})),

s({trig="wrapfig", dscr="Text wrapping around figure.", condition=conds.line_begin*lsu.line_end, snippetType='autosnippet'},
fmta([[\begin{wrapfigure}<>{<>}{0.5\textwidth}
	\centering
	\includegraphics[width=\linewidth]{./figures<>}
	\caption{
	    \textbf{TITLE.}
	    CAPTION.
	}
	\label{fig:<>}
\end{wrapfigure}

]], {c(1, {t"", t"[NUM LINES]"}), c(2, {t"r", t("R",virt("Uppercase=float"))}), i(3, "FILENAME"), i(4)})),


s({trig="tab", dscr="Table", condition=conds.line_begin*lsu.line_end, snippetType="autosnippet"},
-- @{} suppresses space between columns. @{.} would use "." as column separator.
fmta([[\begin{table}[ht]
    \caption{
        \textbf{<>.}
        <>.
    }
    \begin{tabularx}{\textwidth}{@{}<>@{}}
        \toprule
        <> \\
        \midrule
        <> \\
        \bottomrule
    \end{tabularx}
    \label{tab:<>}
\end{table}

]], {
    i(1, [[TITLE]]),
    i(3, [[CAPTION]]),
    i(4, [[lcX]]),
    i(5, [[left & right & filling]]),
    i(6, [[left & right & filling]]),
    i(2, [[REF]]),
})),

-- \cmidrule{2-4} makes a thin line covering column 2 to 4.
s({trig="\\cmidrule", dscr="rule spanning subset of columns", snippetType="autosnippet"},
fmta("\\cmidrule{<>-<>}", {i(1, "FROM"), i(2, "TO")})),

s({trig="\\multicolumn", dscr="cell spanning multiple columns", snippetType="autosnippet"},
fmta("\\multicolumn{<>}{<>}{<>}", {i(1, "3"), i(2, "c"), i(3, "TEXT")})),

s({trig="dark", dscr="Dark mode."},
{t{"% temp dark mode", [[\usepackage{xcolor}\pagecolor{black}\color{white}]]}},
{show_condition=lsu.line_end}),

s({trig="gif", dscr="Animation that works in e.g. pdfpc. Convert gif to mov online."},
fmta([[\href{run:./figures/<>.mov?autostart&loop}{\includegraphics[width=\textwidth]{./figures/<>.png}}]],
{i(1, "FILENAME"), rep(1)}),
{show_condition=lsu.line_end}),

s(
    {trig="(%d)tw", dscr="textwidth", trigEngine="pattern", snippetType="autosnippet"},
    {re(1), t"\\textwidth"}
),
s(
    {trig="(%d)lw", dscr="linewidth", trigEngine="pattern", snippetType="autosnippet"},
    {re(1), t"\\linewidth"}
),

-- if using the acro package
-- s({trig="acro", dscr="Define new acronym", snippetType="autosnippet", condition=conds.line_begin},
--     fmta([[\DeclareAcronym{<>}{short=<>,long=<>}]], {
--         i(1), dUpper(2, 1), dUpper(3, 1),
--     })
-- ),
-- s({trig="glos", dscr="Define new glossary entry", snippetType="autosnippet", condition=conds.line_begin},
--     fmta([[\DeclareAcronym{<>}{preset=glossary,short=<>,long=<>,list=<>}]], {
--         i(1), dUpper(2, 1), dUpper(3, 1), i(4),
--     })
-- ),
-- if using the glossaries or glossaries-extra package
s({trig="acro", dscr="Define new acronym", snippetType="autosnippet", condition=conds.line_begin},
    fmta([[\newacronym<>{<>}{<>}{<>}]], {
            c(4, {t"", {t"[see={[glossary:]{gls-", i(1), t"}}]"}}),
            i(1),
            lsu.upper_jump(2, 1),
            lsu.upper_jump(3, 1),
    })
),
s({trig="glos", dscr="Define new glossary entry", snippetType="autosnippet", condition=conds.line_begin},
    fmta([[\newglossaryentry{<>}{
    name={<>},
    description={<>},
}]], {
        i(1), lsu.rep_jump(2, 1), lsu.upper_jump(3, 1),
    })
),

s({trig="newcommand", dscr="Define new command."},
fmta(
    [[\newcommand{\<>}[<>]<>{<>}]],
    {
        i(1, "NEWNAME"),
        i(2, "NUMBER OF ARGS INCL OPTIONALS"),
        c(3, {t"", t"[OPTIONAL DEFAULT VAL]"}),
        i(4, "BODY WITH #1 (FIRST OPT), #2, ..."),
    }
)),

-- use @ to convert to acronym/glossary, either before or after identifier
s({trig=" @(%w+)(%W)", trigEngine="pattern", snippetType="autosnippet"}, { t" ", f(acro, {}), re(2)}),
s({trig="(%w+)@", trigEngine="pattern", snippetType="autosnippet"}, { f(acro, {})}),
-- use @ as postifx multiple times to toggle acronym version between long, short, and default.
s({trig="(\\[Aa]c)(p?{%w+})@", trigEngine="pattern", snippetType="autosnippet"}, { re(1), t"l", re(2) }),
s({trig="(\\[Aa]c)l(p?{%w+})@", trigEngine="pattern", snippetType="autosnippet"}, { re(1), t"s", re(2) }),
s({trig="(\\[Aa]c)s(p?{%w+})@", trigEngine="pattern", snippetType="autosnippet"}, { re(1), re(2) }),

s({
    trig="pdftex",
    dscr="Include figure with text in separate tex file.",
    -- show_condition=lsu.line_end,
},
fmta([[\pdftex[<>]{<>}{<>}
]], {
    i(1, "\\textwidth"),
    i(2, "figures/"),
    i(3, "filename.pdf"),
})),

-- https://latex-programming.fandom.com/wiki/Minipage_(LaTeX_environment)
s({trig="minipage", dscr="Minipage"},
fmta([[\begin{minipage}<><><>{<>}
    <>
\end{minipage}
]], {
    c(1, {t"", t"[c]", t"[t]", t"[b]"}),
    i(2, t"[\\textheight]"),
    c(3, {t"", t"[t]", t"[c]", t"[b]", t("[s]", virt("[s]tretch"))}),
    i(4, t"\\textwidth"),
    i(5),
})),

}
