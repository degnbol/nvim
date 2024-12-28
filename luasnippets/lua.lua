local ls = require "luasnip"
local util = require "utils/luasnip"
local s = ls.s
local t = ls.t
local i = ls.i

return {
-- core lang snippets similar to other bracket related ones found in all.lua

s({trig="edn", dscr="Typo end", snippetType='autosnippet'},
{t"end"}),

s({trig="stdin", dscr="Read stdin"},
{t[[local text = io.read("*a")]]}),

s({trig="read", dscr="Read text of a file"},
    fmta([[local fh = io.open("<>")
local <> = fh:read("*a") -- *a or *all
fh:close()
]], {i(1), i(2, "content")})),

s({trig="popen", dscr="popen template"},
fmta([[local handle = io.popen("<>")
if handle ~= nil then
    local result = handle:read("*a")
    handle:close()
end
]], {i(1, "COMMAND")})),

-- keep misremembering and LSP doesn't show it.
s({trig=":lower", dscr="Lowercase text."},
{t":lower()"}),
s({trig=":upper", dscr="Uppercase text."},
{t":upper()"}),

s({
    trig='%[%[%[%]%]',
    dscr="Third bracket in trig pattern is typed last.",
    snippetType='autosnippet',
    trigEngine=function (trigger) return util.match_ahead(2) end,
    priority = 1000,
},
-- NOTE: there is no final bracket in replacement since one is left over from end of trigger.
fmta([=[[[
    <>

]=], {i(1)})),

s({trig="join", dscr="Join array of strings to a string."},
{t"table.concat(", i(1,"strings"), t", ", i(2,"[sep]"), t")"}),

-- or: TABLE[#TABLE+1] = VALUE
s({trig="append", dscr="Append to table"},
{t"table.insert(", i(1, "TABLE"), t", ", i(2, "VALUE"), t")"}),

s({trig="rep", dscr="String repeat"},
{t"string.rep(", i(1, "' '"), t", ", i(2, "N"), t")"}),

s({trig="strip", dscr="Strip whitespace"},
{i(1), t[[:match("^[\t%s]*(.-)[\t%s]*$")]]}),

-- meta. snippet to write snippets.
-- TODO: condition these on relevant path in the same way you did it for completions for configuring lazy.
-- Also add all the allowed args from https://github.com/L3MON4D3/LuaSnip/blob/master/DOC.md

s({trig="inspect", dscr="Print inspect", show_condition=conds.line_end},
{t"print(vim.inspect(", i(1), t"))"}),

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

s({trig="condition", dscr="Set condition for autosnippet."},
{t'condition=', i(1, 'conds.line_begin')}),

s({trig="show_condition", dscr="Set show_condition for non-auto snippet."},
{t'show_condition=', i(1, 'conds.line_end')}),

s({trig="fmta", dscr="Convenience function for formatting snippets, especially multi line snippets."},
fmta([=[-- NOTE: no surrounding curly brackets
fmta([[<>
]], {<>})]=], {i(1, "<>"), i(2, "i(1)")})),

-- config snippets

s({trig="req", dscr="require", condition=conds.line_begin, snippetType='autosnippet'},
{c(2, {
    t"",
    f(function (import_name)
        local parts = vim.split(import_name[1][1], '.', true)
        return "local " .. (parts[#parts] or "") .. " = " end, {1})
}), t'require "', i(1), t'"'}),

s({trig="nmap"},
fmta([[vim.keymap.set('n', '<>', "<<Cmd>><><<CR>>", { desc="<>" })]],
{i(1, "<leader>"), i(2), i(3)}),
{condition=conds.line_begin}),
s({trig="imap"},
fmta([[vim.keymap.set('i', '<>', "<<Cmd>><><<CR>>", { desc="<>" })]],
{i(1, "<leader>"), i(2), i(3)}),
{condition=conds.line_begin}),
s({trig="vmap"},
fmta([[vim.keymap.set('v', '<>', "<<Cmd>><><<CR>>", { desc="<>" })]],
{i(1, "<leader>"), i(2), i(3)}),
{condition=conds.line_begin}),
s({trig="xmap"},
fmta([[vim.keymap.set('x', '<>', "<<Cmd>><><<CR>>", { desc="<>" })]],
{i(1, "<leader>"), i(2), i(3)}),
{condition=conds.line_begin}),
-- for cmdline, e.g. abbreviations
s({trig="cmap"},
fmta([[vim.keymap.set('c', '<>', "<<Cmd>><><<CR>>", { desc="<>" })]],
{i(1, "<leader>"), i(2), i(3)}),
{condition=conds.line_begin}),
s({trig="omap"},
fmta([[vim.keymap.set('o', '<>', "<<Cmd>><><<CR>>", { desc="<>" })]],
{i(1, "<leader>"), i(2), i(3)}),
{condition=conds.line_begin}),
s({trig="nxmap"},
fmta([[vim.keymap.set({'n', 'x'}, '<>', "<<Cmd>><><<CR>>", { desc="<>" })]],
{i(1, "<leader>"), i(2), i(3)}),
{condition=conds.line_begin}),

s({trig="rtp", dscr="Get neovim config root", condition=conds.line_begin, snippetType='autosnippet'},
{t"local rtp = vim.opt.runtimepath:get()[1]"}),

s({trig="setlocal", dscr="Set local option.", condition=conds.line_begin*conds.line_end},
{t"vim.opt_local."}),

s({trig="count", dscr="Count prefix for functions."},
{t{"local count = vim.v.count", "if count == 0 then count = 1 end"}}),

s({trig="isfile", dscr="Is file readable neovim function."},
{t"vim.fn.filereadable(", i(1,"FILENAME"), t")"}),

-- api

s({trig="nvim", dscr="Api functions", trigEngine="pattern", snippetType='autosnippet', condition=conds.line_begin},
{t"vim.api.nvim"}),

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

s({trig="vim.api.[%a._]*set_cursor", dscr="set cursor", trigEngine='pattern', snippetType="autosnippet", condition=conds.line_begin},
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

s({trig="command", dscr='equivalent of "command" in vimscript.'},
{t'vim.api.nvim_create_user_command("', i(1,"Q"), t'", "', i(2,"q"), t'", {})'}),
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
vim.api.nvim_feedkeys(keys, '<>', false)
]], {i(1, "right"), c(2, {t"m", t"n", t"t", t"i", t"x", t"x!"})})),

-- returns stdout
s({
    trig="sync",
    dscr="call external cmd synchronously",
    snippetType='autosnippet',
    condition=conds.line_begin,
},
{t'local obj = vim.system(', i(1, "{'echo', 'hello'}"), t', {text=true}):wait()'}),
s({
    trig="async",
    dscr="call external cmd asynchronously",
    snippetType='autosnippet',
    condition=conds.line_begin,
}, fmta([[vim.system({<>}, {text=true}, function(obj)
    if obj.code ~= 0 then
        -- needs local util = require "utils/init"
        util.schedule_notify(obj)
    else
        <>
    end
end)
]], {i(1, "'echo', 'hello'"), i(2, {'print(obj.signal)',
                '        print(obj.stdout)',
                '        print(obj.stderr)'
            })})),

s({trig="cword", dscr="current word under cursor"},
{t'vim.fn.expand("<cword>")'}),

s({trig="defer", dscr="Defer call", condition=conds.line_begin},
fmta([[vim.defer_fn(function ()
    <>
end, <>)
]], {i(1), i(2, "500")})),

s({trig="mode", dscr="Get mode (normal or visual etc)", show_condition=conds.line_end},
{t"local mode = vim.api.nvim_get_mode().mode"}),

s({trig="dynamic", dscr="Dynamic node"},
fmta([[d(<>, function (args<>)
    -- https://github.com/L3MON4D3/LuaSnip/blob/master/DOC.md#dynamicnode
    -- the returned snippetNode doesn't need a position; it's inserted "inside" the dynamicNode.
    return sn(nil, {
        -- jump-indices are local to each snippetNode, so restart at 1.
        i(1, args[1][1])
    })
end, {<>}<>)
]], {
    i(1, "JUMP"),
    c(3, {t"", t", parent, old_state, user_args"}),
    i(2, "NODEREF"),
    c(4, {t"", t", {user_args={}}"}),
})),

}

