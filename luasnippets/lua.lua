#!/usr/bin/env lua
return {
-- meta. snippet to write snippets.

s({trig="util", dscr="req utils", condition=conds.line_begin},
{t[[local lsu = require "utils/luasnip"]]}),

s({trig="snip", snippetType="autosnippet"},
fmta([[s({trig="<>", dscr="<>"<>},
{t"<>"}),
]], {i(1, "TRIGGER"), i(2, "DESCRIPTION"),
c(3, {t"", t", snippetType='autosnippet'", t", wordTrig=true, regTrig=false, snippetType='snippet'"}),
i(4, "EXPANSION")}),
{condition=conds.line_begin}),

s({trig="ausnip", snippetType="autosnippet"},
t[[snippetType="autosnippet"]]),

s({trig="condb", snippetType="autosnippet"},
t"condition=conds.line_begin"),

-- config snippets

s({trig="req", dscr="require", condition=conds.line_begin, snippetType='autosnippet'},
{c(2, {
    f(function (import_name)
        local parts = vim.split(import_name[1][1], '.', true)
        return "local " .. (parts[#parts] or "") .. " = " end, {1}),
    t""
}), t'require "', i(1), t'"'}),

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

s({trig="root", dscr="Get neovim config root", condition=conds.line_begin, snippetType='autosnippet'},
{t"local rtp = vim.opt.runtimepath:get()[1]"}),

-- api

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

s({trig="cursor", dscr="get cursor position", condition=conds.line_begin},
{t"local r, c = unpack(vim.api.nvim_win_get_cursor(0))"}),

s({trig="[%a._]*set_cursor", dscr="set cursor", regTrig=true, snippetType="autosnippet", condition=conds.line_begin},
{t"vim.api.nvim_win_set_cursor(0, {", i(1, "r, c"), t"})"}),

s({trig="current", dscr="get current line text"},
{t"local line = vim.api.nvim_get_current_line()"}),

-- -1 since nvim_win_get_cursor is (1,0)-indexed and nvim_buf_set_text is 0-indexed.
s({trig="char", dscr="get char under cursor"},
{t"vim.api.nvim_buf_get_text(0, r-1, c-1, r-1, c, {})[1]"}),

s({trig="get_mark", dscr="get buffer mark", condition=conds.line_begin, snippetType="autosnippet"},
{t[[local r_mark, c_mark = unpack(vim.api.nvim_buf_get_mark(0, "]], i(1), t'"))'}),

s({trig="buftype", dscr="buffer type"},
t"vim.bo.buftype"),
s({trig="filetype", dscr="file type"},
t"vim.bo.filetype"),

s({trig="filepath", dscr="get filepath for current buffer"},
{t"vim.api.nvim_buf_get_name(0)"}),

s({trig="startswith", dscr="startswith"},
{t"vim.startswith(", i(1), t")"}),

s({trig="endswith", dscr="endswith"},
{t"vim.endswith(", i(1), t")"}),

}
