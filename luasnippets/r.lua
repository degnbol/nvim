#!/usr/bin/env lua
local lsu = require "utils/luasnip"
local virt = lsu.virt

return {
--
s({trig="lib", dscr="library", snippetType="autosnippet"},
{
    c(1, {t("", virt("^l to suppress messages")),
    t"suppressPackageStartupMessages("}),
    t"library(", i(2, "package"), t")",
    f(function (node_text)
        if node_text[1][1] == "" then
            return "" else return ")"
        end
    end, {1}), -- provide text from node 1 (choice node)
},
{condition=conds.line_begin}),

s({trig="root", dscr="root", snippetType="autosnippet"},
t{"suppressPackageStartupMessages(library(here));ROOT=here()", ""},
{condition=conds.line_begin}),

s({trig="melt", dscr="melt", snippetType="autosnippet"},
fmta([[<> = melt(<>, <>, <>, variable.name="<>", value.name="<>")]],
{i(1, "dtm"), i(2, "dt"), i(3, "id.vars"), i(4, "measure.vars"), i(5, "variable"), i(6, "value")}),
{condition=conds.line_begin}),


-- data.table

s({trig="DTupdate", dscr="Update joined rows"},
fmta([[<>[<>, on=c(<>="<>"), c("<>") := .(<>)]
]], {i(1, "DT_TRG"), i(2, "DT_SRC"), i(3, "ON_TRG"), i(4, "ON_SRC"), i(5, "COL_TRG"), i(6, "i.COL_DUP")})),

s({trig="DTmulti", dscr="Multiple column assignment"},
{t'c("', i(1, 'COL'), t'"):=.(', i(2, 'COL'), t")"}),

-- ggplot

s({trig="ggplot", dscr="ggplot", condition=conds.line_begin, snippetType="autosnippet"},
fmta([[ggplot(<>, aes(<>)) +
    <>
]], { i(1, "dt"), i(2, "x=x"), i(3, "geom_point()")})),

s({trig="expandy0", dscr="Remove the default empty space under y=0"},
t"scale_y_continuous(expand=expansion(mult=c(0,.05)))", -- .05 based on ?expansion
{condition=conds.line_begin}),

-- altho you can annotate anything, not just text.
s({trig="annotate", dscr="annotate text"},
fmta([[annotate("text", label="<>", x=<>, y=<>, size=<>, hjust=0) +]],
{i(1, "LABEL"), i(2, "x"), i(3, "y"), i(4, "4")})),

}
