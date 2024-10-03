local ls = require "luasnip"
local s = ls.s
local lsu = require "utils/luasnip"
local virt = lsu.virt

return {
    s({trig="lib", dscr="library", snippetType="autosnippet"},
    {
        c(1, {t(""),
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

    s({trig="DTagg", dscr="Aggregation with naming"},
    fmta([[<>[<>, .(<>=<>), by=c("<>")]
    ]], {i(1, "dt"), i(2), i(3, "NAME"), i(4, "sum(COL)"), i(5, "BY")})),

    s({trig="dcast", dscr="Opposite of melt. Long to wide table reshaping.", show_condition=conds.line_end},
        fmta([[dcast(<>, ... ~ <>, value.var="<>")]], {
            i(1, "dt"),
            i(2, "VARIABLE"),
            i(3, "VALUE"),
        })),

    -- ggplot

    s({trig="ggplot", dscr="ggplot", show_condition=conds.line_end},
    fmta([[ggplot(<>, aes(<>)) +
        geom_<>() +
        theme_bw() +
        theme(<>
        )]], {
        i(1, "dt"),
        i(2, "x=x, color=color, fill=fill"),
        c(3, {t"point", t"line", t"col", t"histogram", t"boxplot"}),
        i(4),
    })),

    s({trig="theme_", dscr="Complete themes"}, {
        t"theme_",
        c(1, {
            -- https://ggplot2.tidyverse.org/reference/ggtheme.html#ref-examples
            t("bw",       virt("Gray grid, black border, gray facet strip")),
            t("linedraw", virt("Thin black lines, black facet strip")),
            t("classic",  virt("No panel grid, white facet strip")),
            t("light",    virt("Gray grid, gray border, gray facet strip with white text")),
            t("minimal",  virt("No panel border")),
            t("void",     virt("No panel decorations nor axis")),
        }),
        t"(",
        i(2, {"", "    base_size=11,", "    base_family=\"\"", ""}),
        t")",
    }),

    -- s({trig="=ele", dscr="Elements", snippetType="autosnippet", wordTrig=false}, {
    s({trig="element_", dscr="Elements"}, {
        t"=element_",
        c(1, {t"blank", t"text", t"line", t"rect", t"part_rect"}),
        t"()",
    }),
    s({trig="blank", dscr="element_blank"}, { t"element_blank(),", }),

    s({trig="expandy0", dscr="Remove the default empty space under y=0"},
    t"scale_y_continuous(expand=expansion(mult=c(0,.05)))", -- .05 based on ?expansion
    {condition=conds.line_begin}),

    -- altho you can annotate anything, not just text.
    s({trig="annotate", dscr="annotate text"},
    fmta([[annotate("text", label="<>", x=<>, y=<>, size=<>, hjust=0) +]],
    {i(1, "LABEL"), i(2, "x"), i(3, "y"), i(4, "4")})),

    s({trig="facet_grid", dscr="Facet grid"},
    fmta([[facet_grid(rows=vars(<>), cols=vars(<>))]], {i(1, "GROUP1"), i(2, "GROUP2")})
    ),

    s({trig="facet_wrap", dscr="Facet wrap"},
        fmta([[facet_wrap(vars(<>), <>=<>)]], {
            i(1, "GROUP"),
            c(2, {t"nrow", t"ncol"}),
            i(3, "2"),
        })
    ),
}
