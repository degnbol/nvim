-- should be unnessary but s keeps being understood as a bool set somewhere
local ls = require "luasnip"
local s = ls.s

return {
-- not auto

s({trig="cont", dscr="Continue."},
{t"continue"}),
s({trig="ret", dscr="Return."},
{t"return "}),

s({trig="ROOT", dscr="Git root.", show_condition=conds.line_end},
{
    t{[[ROOT = `git root` |> readchomp]], ""},
    c(1, {t"include", t"cd"}),
    t[[("$ROOT/]], i(2), t[[")]],
}),

-- replace dumb pr -> print("<>") from one of the snippet bundles
s({trig="pr", dscr="println"},
{t"println(", i(1), t")"}),

s({trig="amat", dscr="AbstractMatrix"}, {t"AbstractMatrix"}),
s({trig="avec", dscr="AbstractVector"}, {t"AbstractVector"}),
s({trig="astr", dscr="AbstractString"}, {t"AbstractString"}),
-- for val::AbstractMatrix or T<:AbstractMatrix
s({trig=":amat", dscr="AbstractMatrix", snippetType="autosnippet"}, {t":AbstractMatrix"}),
s({trig=":avec", dscr="AbstractVector", snippetType="autosnippet"}, {t":AbstractVector"}),
s({trig=":astr", dscr="AbstractString", snippetType="autosnippet"}, {t":AbstractString"}),

s({trig="SparseMatrix", dscr="Sparse matrix"}, {t"SparseMatrixCSC"}),


s({trig="agg", dscr="Aggregate DataFrame", condition=conds.line_end},
{t"combine(groupby(", i(1,"df"), t", ", i(2, ":COL"), t"), ", i(3, ":COL2 => sum"), t")"}),

}, {
-- auto

-- overwrite the one from https://github.com/honza/vim-snippets
s({trig="#!", dscr="shebang", condition=conds.line_begin * conds.line_end},
t{"#!/usr/bin/env julia", ""}),

s({
    trig="using df",
    dscr="Using DataFrames and CSV.",
    wordTrig=false,
    condition=conds.line_begin * conds.line_end,
},
{t{"using DataFrames, CSV", ""}}),

s({trig="DF", dscr="DataFrame"}, {t"DataFrame"}),

s({trig="CSV.read", dscr="CSV read DataFrame", condition=conds.line_end},
{t'CSV.read("', i(1, "FILENAME.tsv.gz"), t[[", DataFrame; delim='\t')]]}),

s({trig="CSV.write", dscr="CSV write DataFrame", condition=conds.line_end},
{
    t'CSV.write("',
    i(1),
    t'.tsv.gz", df',
    i(2),
    t"; delim='\\t'",
    c(3, {t", compress=true", t""}),
    t')'
}),

s({trig="<=", dscr="Less than."}, {t"≤"}),
s({trig=">=", dscr="Greater than."}, {t"≥"}),
-- we don't do this one, since it's similar to == which is long,
-- whereas >= is similar to > which is short.
-- s({trig="!=", dscr="Not equal."}, {t"≠"}),
s({trig=".in", dscr="Each in"}, {t".∈ Ref(", i(1), t")"}),

s({trig="use", dscr="Using.", condition=conds.line_begin * conds.line_end}, {t"using "}),

s({trig="func", dscr="Function.", condition=conds.line_begin * conds.line_end},
fmta([[function <>(<>)<>
end]], {i(1, "NAME"), i(2), i(3)})),

-- TODO: not in comment, not when PlotlyJS isn't loaded.
s({trig="={", wordTrig=false, dscr="PlotlyJS attributes, assuming a=attr has been called."},
{t"=…(", i(1), t")"}),


}
