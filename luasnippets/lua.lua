#!/usr/bin/env lua
return {
-- meta. snippet to write snippets.

s({trig="util", dscr="req utils", snippetType="autosnippet"},
{t[[local lsu = require"luasnip_util"]]}),

s({trig="snip", snippetType="autosnippet"},
fmta([[s({trig="<>", dscr="<>", <>},
{t"<>"}),
]], {i(1, "trigger"), rep(1), i(2, "wordTrig/regTrig/snippetType"), i(3, "expansion")}),
{condition=conds.line_begin}),

s({trig="ausnip", snippetType="autosnippet"},
t[[snippetType="autosnippet"]]),

s({trig="condb", snippetType="autosnippet"},
t"condition=conds.line_begin"),

-- config snippets

s({trig="aucmd", snippetType="autosnippet" },
fmta([[
vim.api.nvim_create_autocmd(<>, {
    pattern = <>,
    callback = <>
})
]],
{i(1, "Event(s)"), i(2, '"*.EXT"'), i(3, '"string cmd or lua func"')}),
{condition=conds.line_begin}),
s({trig="augroup", snippetType="autosnippet" },
-- clear=true is the default. You need to supply at least empty {}.
fmta('local <> = vim.api.nvim_create_augroup("<>", {clear=true})',
{rep(1), i(1, "name")}),
{condition=conds.line_begin}),

s({trig="nmap", snippetType="autosnippet"},
fmta([[vim.keymap.set("n", "<>", "<>")]],
{i(1, "from"), i(2, "to")}),
{condition=conds.line_begin}),
s({trig="imap", snippetType="autosnippet"},
fmta([[vim.keymap.set("i", "<>", "<>")]],
{i(1, "from"), i(2, "to")}),
{condition=conds.line_begin}),
s({trig="vmap", snippetType="autosnippet"},
fmta([[vim.keymap.set("v", "<>", "<>")]],
{i(1, "from"), i(2, "to")}),
{condition=conds.line_begin}),
s({trig="xmap", snippetType="autosnippet"},
fmta([[vim.keymap.set("x", "<>", "<>")]],
{i(1, "from"), i(2, "to")}),
{condition=conds.line_begin}),

}