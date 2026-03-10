---@diagnostic disable: unused-local
local lsu = require "utils/luasnip"
local s, t, i, c, f, d, sn, fmta, conds, rep, ms = lsu.s, lsu.t, lsu.i, lsu.c, lsu.f, lsu.d, lsu.sn, lsu.fmta, lsu.conds, lsu.rep, lsu.ms
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

