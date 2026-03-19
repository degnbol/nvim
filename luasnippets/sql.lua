---@diagnostic disable: unused-local
local ls = require "utils/luasnip"
local s, t, i, c, f, d, sn, fmta, conds, rep, ms = ls.s, ls.t, ls.i, ls.c, ls.f, ls.d, ls.sn, ls.fmta, ls.conds, ls.rep, ls.ms
return {
    s({trig="read", dscr="Read table"},
    fmta([[
    -- TBL should already be created
    \copy <> from '<>' DELIMITER <> CSV HEADER;
    ]], {i(1, "TBL"), i(2, "FILENAME.tsv"), i(3, "E'\\t'")})),
    s({trig="write", dscr="Write table"},
    fmta([[
    \copy <> to '<>' DELIMITER <> CSV HEADER;
    ]], {i(1, "TBL"), i(2, "FILENAME.tsv"), i(3, "E'\\t'")})),
}
