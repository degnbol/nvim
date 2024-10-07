local ls = require "luasnip"
local s = ls.s

local lsu = require "utils/luasnip"
local re = lsu.re

return {
}, {
-- autosnippets

s({trig="(%s*)([%w-]+):", dscr="Auto quote unquoted key.",
trigEngine="pattern"},
{re(1), t'"', re(2), t'": '}),

s({trig=": ?{{", dscr="Opening bracket after a key.",
trigEngine="pattern"},
-- two spaces
{t{": {", "  "}, i(1), t{"", "}"}}),

}
