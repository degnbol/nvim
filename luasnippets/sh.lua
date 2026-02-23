local ls = require "luasnip"
local s = ls.s

return {
--

s({trig="#!", dscr="Zsh shebang.", snippetType="autosnippet", condition=conds.line_begin},
{t{"#!/usr/bin/env zsh", ""}}),

s({trig="exists", dscr="Does it exist?"},
{t"command -v ", i(1, "cargo"), t" > /dev/null"}),

s({trig="ifexists", dscr="If install tool is avail"},
fmta([[
if command -v <> >> /dev/null; then
    <> install <>
else
    <> <>
fi
]], {i(1, "uv"), rep(1), i(2, "PACKAGE"), i(3, "mamba install -yc conda-forge"), rep(2)})),


s({trig="ifmac", dscr="If MacOS"},
fmta([=[if [[ "$OSTYPE" == darwin* ]]; then
    <>
fi
]=], {i(1)})),

}
