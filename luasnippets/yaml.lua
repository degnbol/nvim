local ls = require "luasnip"
local s = ls.snippet
local sn = ls.snippet_node
local isn = ls.indent_snippet_node
local t = ls.text_node
local i = ls.insert_node
local f = ls.function_node
local c = ls.choice_node
local d = ls.dynamic_node
local r = ls.restore_node
local events = require("luasnip.util.events")
local ai = require("luasnip.nodes.absolute_indexer")
local extras = require("luasnip.extras")
local l = extras.lambda
local rep = extras.rep
local p = extras.partial
local m = extras.match
local n = extras.nonempty
local dl = extras.dynamic_lambda
local fmt = require("luasnip.extras.fmt").fmt
local fmta = require("luasnip.extras.fmt").fmta
local conds = require("luasnip.extras.expand_conditions")
local postfix = require("luasnip.extras.postfix").postfix
local types = require("luasnip.util.types")
local parse = require("luasnip.util.parser").parse_snippet
local ms = ls.multi_snippet
local k = require("luasnip.nodes.key_indexer").new_key

local lsu = require "utils/luasnip"
local re = lsu.re

---https://github.com/L3MON4D3/LuaSnip/blob/master/DOC.md#dynamicnode
local function _dUpper(args)
    -- the returned snippetNode doesn't need a position; it's inserted
    -- "inside" the dynamicNode.
    return sn(nil, {
        -- jump-indices are local to each snippetNode, so restart at 1.
        i(1, args[1][1]:upper())
    })
end
local function dUpper(jump_index, node_reference)
    return d(jump_index, _dUpper, {node_reference})
end

return {
s({trig="glossy", dscr="New glossy entry", condition=conds.line_begin},
fmta([[<>:
  short: <>
  long: <>
  <>
]], {i(1, "key"), dUpper(2, 1), i(3, ""), c(4, {t"", t"description: "})})),
}

