#!/usr/bin/env lua
local ls = require "luasnip"
local extras = require "luasnip/extras"
local s = ls.s
local t = ls.t
local f = ls.f
local c = ls.c
local fmta = extras.fmta

local lsu = require "utils/luasnip"
local re = lsu.re

local snippets = {}

for _, trig in ipairs({
    "^([\t ]*)(%w+) ?=",
    "^([\t ]*)set (%w+) ?=",
}) do
    table.insert(snippets, s({
        trig=trig,
        trigEngine="pattern",
        snippetType='autosnippet',
        dscr="Variable assignment",
    },
        {re(1), t"set ", re(2), t" to"}))
end

return snippets
