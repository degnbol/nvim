#!/usr/bin/env lua
local M = {}

function M.in_math()
  return vim.fn['vimtex#syntax#in_mathzone']() == 1
end
function M.in_text() return not in_math() end
function M.in_comment()
  return vim.fn['vimtex#syntax#in_comment']() == 1
end
M.is_inside = vim.fn['vimtex#env#is_inside']
function M.in_env(name)
    local r, c = unpack(is_inside(name))
    return (r > 0 and c > 0)
end
function M.in_itemize() return in_env('itemize') end

return M
