#!/usr/bin/env lua

return {
--

s({trig="installed", dscr="is it installed?"},
{t"which ", i(1, "cargo"), t" > /dev/null"}),

s({trig="ifinstalled", dscr="is it installed?"},
fmta([[
if which <> >> /dev/null; then
    <> install <>
else
    <> <>
fi
]], {i(1, "pipx"), rep(1), i(2, "PACKAGE"), i(3, "mamba install -yc conda-forge"), rep(2)})),


}