local ls = require "luasnip"
local s = ls.s

return {
-- not auto

s({trig="esc", dscr="Temp ESC"},
{t[[<C-\><C-O>]]}),

s({trig="syntax", dscr="Syntax capture.", trigEngine="pattern", show_condition=conds.line_end},
{
    t"syntax ",
    c(1, {t"match", t"region"}),
    t" ",
    i(2, "HLGROUP"),
    t" ",
    d(3, function (args)
        if args[1][1] == "match" then
            return sn(nil, { -- the returned snippetNode doesn't need a position; it's inserted "inside" the dynamicNode.
                t"/",
                -- jump-indices are local to each snippetNode, so restart at 1.
                i(1, "PATTERN"),
                t"/",
            })
        else
            return sn(nil, { -- the returned snippetNode doesn't need a position; it's inserted "inside" the dynamicNode.
                c(1, {t"", {t"matchgroup=", i(1, "HLSTARTEND"), t" "}}),
                t"start=/",
                -- jump-indices are local to each snippetNode, so restart at 1.
                i(2, "PATTERN"),
                t"/ end=/",
                i(3, "PATTERN"),
                t"/",
            })
        end
    end, {1}),
    c(4, {t"", {t" containedin=", i(1, "A,B,C")}})
}),

}, {
-- auto

-- Since it is a wordTrig it will not get corrected in valid use, such as pathogen#infect()
-- TODO: make a condition for being inside a comment and don't exclude that with condition=-in_comment
s({trig="#", dscr="Correct wrong comment char."},
{t'"'}),

}
