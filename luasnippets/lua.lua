#!/usr/bin/env lua
local util = require "utils/luasnip"

return {
-- core lang snippets similar to other bracket related ones found in all.lua
s({trig="[[[", dscr="Double square brackets", snippetType='autosnippet'},
{t"[[", i(1), t"]]"}),
-- TODO: make a version for pressing [[[[ that adds newline or make pressing 
-- enter like pressing escape then O

-- meta. snippet to write snippets.

s({trig="util", dscr="req utils", condition=conds.line_begin},
{t[[local lsu = require "utils/]], i(1, "luasnip"), t'"'}),

s({trig="snip", snippetType="autosnippet"},
fmta([[s({trig="<>", dscr="<>"<>},
{t"<>"}),
]], {i(1, "TRIGGER"), i(2, "DESCRIPTION"),
c(3, {t"", t", snippetType='autosnippet'", t", wordTrig=true, trigEngine='pattern', snippetType='snippet'"}),
i(4, "EXPANSION")}),
{condition=conds.line_begin}),

s({trig="ausnip", snippetType="autosnippet"},
t[[snippetType="autosnippet"]]),

s({trig="condb", snippetType="autosnippet"},
t"condition=conds.line_begin"),

s({trig="wordTrig", dscr="Word trigger?"},
{t"wordTrig=false"}),

s({trig="trigEngine", dscr="Set trigger engine."},
{t'trigEngine="', c(1, {t"pattern", t"vim", t"ecma", t"plain"}), t'"'}),

s({trig="fmta", dscr="Convenience function for formatting snippets, especially multi line snippets."},
fmta([=[fmta([[<>
]], {<>})]=], {i(1, "<>"), i(2, "i(1)")})),

-- config snippets

s({trig="req", dscr="require", condition=conds.line_begin, snippetType='autosnippet'},
{c(2, {
    t"",
    f(function (import_name)
        local parts = vim.split(import_name[1][1], '.', true)
        return "local " .. (parts[#parts] or "") .. " = " end, {1})
}), t'require "', i(1), t'"'}),

s({trig="nmap", snippetType="autosnippet"},
fmta([[vim.keymap.set('n', '<>', "<<Cmd>><><<CR>>", { desc="<>" })]],
{i(1, "<leader>"), i(2), i(3)}),
{condition=conds.line_begin}),
s({trig="imap", snippetType="autosnippet"},
fmta([[vim.keymap.set('i', '<>', "<<Cmd>><><<CR>>", { desc="<>" })]],
{i(1, "<leader>"), i(2), i(3)}),
{condition=conds.line_begin}),
s({trig="vmap", snippetType="autosnippet"},
fmta([[vim.keymap.set('v', '<>', "<<Cmd>><><<CR>>", { desc="<>" })]],
{i(1, "<leader>"), i(2), i(3)}),
{condition=conds.line_begin}),
s({trig="xmap", snippetType="autosnippet"},
fmta([[vim.keymap.set('x', '<>', "<<Cmd>><><<CR>>", { desc="<>" })]],
{i(1, "<leader>"), i(2), i(3)}),
{condition=conds.line_begin}),
-- for cmdline, e.g. abbreviations
s({trig="cmap", snippetType="autosnippet"},
fmta([[vim.keymap.set('c', '<>', "<<Cmd>><><<CR>>", { desc="<>" })]],
{i(1, "<leader>"), i(2), i(3)}),
{condition=conds.line_begin}),
s({trig="omap", snippetType="autosnippet"},
fmta([[vim.keymap.set('o', '<>', "<<Cmd>><><<CR>>", { desc="<>" })]],
{i(1, "<leader>"), i(2), i(3)}),
{condition=conds.line_begin}),
s({trig="nxmap", snippetType="autosnippet"},
fmta([[vim.keymap.set({'n', 'x'}, '<>', "<<Cmd>><><<CR>>", { desc="<>" })]],
{i(1, "<leader>"), i(2), i(3)}),
{condition=conds.line_begin}),

s({trig="root", dscr="Get neovim config root", condition=conds.line_begin, snippetType='autosnippet'},
{t"local rtp = vim.opt.runtimepath:get()[1]"}),

-- api

s({trig="aucmd", snippetType="autosnippet" },
fmta([[
vim.api.nvim_create_autocmd("<>", {
    pattern = "<>",
    group = grp,
    callback = <>
})
]],
{i(1, "BufEnter"), i(2, '*'), i(3, '"string cmd or lua func"')}),
{condition=conds.line_begin}),
s({trig="augroup", snippetType="autosnippet" },
-- clear=true is the default. You need to supply at least empty {}.
fmta([[local <> = vim.api.nvim_create_augroup("<>", {clear=true})
]],
{c(2, {t"grp", rep(1)}), i(1, "name")}),
{condition=conds.line_begin}),

s({trig="cursor", dscr="get cursor position", condition=conds.line_begin},
{t"local r, c = unpack(vim.api.nvim_win_get_cursor(0))"}),

s({trig="[%a._]*set_cursor", dscr="set cursor", trigEngine='pattern', snippetType="autosnippet", condition=conds.line_begin},
{t"vim.api.nvim_win_set_cursor(0, {", i(1, "r, c"), t"})"}),

s({trig="current", dscr="get current line text"},
{t"local line = vim.api.nvim_get_current_line()"}),
s({trig="get_lines", dscr="get all lines as list of strings"},
{t"local lines = vim.api.nvim_buf_get_lines(0, 0, -1, false)"}),

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

s({trig="set_hl", dscr="set highlight group values"},
fmta('vim.api.nvim_set_hl(0, "<>", {<>})',
{i(1, "Comment"), i(2, 'fg="gray"')})),

s({trig="get_hl", dscr="get highlight group values"},
fmta([[vim.api.nvim_get_hl(0, {name="<>", link=false})['<>']
]], {i(1, "Comment"), i(2, "fg")})),

s({trig="g:", dscr="Set global variable.", condition=conds.line_begin, snippetType='autosnippet'},
{t'vim.api.nvim_set_var("', i(1), t'", ', i(2), t')'}),

s({trig="startswith", dscr="startswith"},
{t"vim.startswith(", i(1), t")"}),

s({trig="endswith", dscr="endswith"},
{t"vim.endswith(", i(1), t")"}),

s({trig="cmd", dscr="set new command for neovim's cmdline."},
fmta([[vim.api.nvim_create_user_command("<>", function ()
    <>
end, {})
]], {i(1, "NAME"), i(2)})),

-- https://www.reddit.com/r/neovim/comments/104lc26/how_can_i_press_escape_key_using_lua/
s({
    trig="press",
    condition=conds.line_begin,
    dscr="Press key where codes such as <right> are converted to some coded version first.",
},
fmta([[local keys = vim.api.nvim_replace_termcodes('<<<>>>', true,false,true)
vim.api.nvim_feedkeys(keys, 'm', false)
]], {i(1, "right")})),

-- help text says to prefer vim.system for both, but it doesn't seem to exist.
-- returns stdout

s({
    trig="sync",
    dscr="call external cmd synchronously",
    snippetType='autosnippet',
    condition=conds.line_begin,
},
{t'vim.fn.system("', i(1, "ls"), t'")'}),
-- Returns a code (not exit code).
-- Use e.g. vim.fn.jobstart("cmd", {on_stdout=...})
s({
    trig="async",
    dscr="call external cmd asynchronously",
    snippetType='autosnippet',
    condition=conds.line_begin,
},
{t'vim.fn.jobstart("', i(1, "ls"), t'")'}),

s({trig="defer", dscr="Defer call", condition=conds.line_begin},
fmta([[vim.defer_fn(function ()
    <>
end, <>)
]], {i(1), i(2, "500")})),

}
