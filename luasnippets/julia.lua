#!/usr/bin/env lua
return {
-- not auto

s({trig="cont", dscr="Continue."},
{t"continue"}),
s({trig="ret", dscr="Return."},
{t"return "}),

s({trig="ROOT", dscr="Git root."},
{
    t{[[ROOT = `git root` |> readchomp]], ""},
    c(1, {t"include", t"cd"}),
    t[[("$ROOT/]], i(2), t[[")]],
}),

-- replace dumb pr -> print("<>") from one of the snippet bundles
s({trig="pr", dscr="println"},
{t"println(", i(1), t")"}),


s({trig="usedf", dscr="Using DataFrames and CSV."},
{t{"using DataFrames, CSV", ""}}),

}, {
-- auto

-- overwrite the one from https://github.com/honza/vim-snippets
s({trig="#!", dscr="shebang", snippetType="autosnippet", condition=conds.line_begin},
t{"#!/usr/bin/env julia", ""}),

s({trig="DF", dscr="DataFrame", snippetType="autosnippet"},
{t"DataFrame"}),

s({trig="CSV.read", dscr="CSV read DataFrame", snippetType="autosnippet"},
{t'CSV.read("', i(1, "FILENAME.tsv.gz"), t[[", DataFrame; delim='\t', quoted=false)]]}),

s({trig="CSV.write", dscr="CSV write DataFrame", snippetType="autosnippet"},
{
    t'CSV.write("',
    i(1),
    t'.tsv.gz", df',
    i(2),
    t"; delim='\\t'",
    c(3, {t", compress=true", t""}),
    t')'
}),

}
