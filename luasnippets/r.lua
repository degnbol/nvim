#!/usr/bin/env lua
local lsu = require"luasnip_util"
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

-- ggplot

s({trig="ggplot", dscr="ggplot", condition=conds.line_begin, snippetType="autosnippet"},
fmta([[ggplot(<>, aes(<>)) +
    <>
]], { i(1, "dt"), i(2, "x=x"), i(3, "geom_point()")})),

s({trig="expandy0", dscr="Remove the default empty space under y=0"},
t"scale_y_continuous(expand=expansion(mult=c(0,.05)))", -- .05 based on ?expansion
{condition=conds.line_begin}),


}
