
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
