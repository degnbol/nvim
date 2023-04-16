#!/usr/bin/env lua
local M = {}

local ls = require("luasnip")
local sn = ls.snippet_node
local i = ls.insert_node
local f = ls.f

function M.in_math()
  -- requires the VimTeX plugin
  return vim.fn['vimtex#syntax#in_mathzone']() == 1
end
function M.in_text() return not in_math() end
function M.in_comment()
  return vim.fn['vimtex#syntax#in_comment']() == 1
end
function M.in_env(name)
    local is_inside = vim.fn['vimtex#env#is_inside'](name)
    return (is_inside[1] > 0 and is_inside[2] > 0)
end
function M.in_itemize() return in_env('itemize') end

-- Summary: When `SELECT_RAW` is populated with a visual selection, the function
-- returns an insert node whose initial text is set to the visual selection.
-- When `SELECT_RAW` is empty, the function simply returns an empty insert node.
function M.get_visual(args, parent)
  if (#parent.snippet.env.SELECT_RAW > 0) then
    return sn(nil, i(1, parent.snippet.env.SELECT_RAW))
  else  -- If SELECT_RAW is empty, return a blank insert node
    return sn(nil, i(1))
  end
end

function M.regGroup(_, snip, i) return snip.captures[i] end
function M.re(i) return f(M.regGroup, nil, {user_args={i}}) end

return M
