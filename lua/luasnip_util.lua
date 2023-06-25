#!/usr/bin/env lua
local M = {}

local ls = require "luasnip"
local sn = ls.snippet_node
local i = ls.insert_node
local f = ls.f

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

-- usage example: c(1, {t("", virt("^l for alt choice")), t"alt choice"})
function M.virt(text)
    return {node_ext_opts={passive={virt_text={{text, "@comment"}}}}}
end

return M
