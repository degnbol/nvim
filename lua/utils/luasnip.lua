#!/usr/bin/env lua
local M = {}

local ls = require "luasnip"
local extras = require "luasnip.extras"
-- otherwise we get errors in other scripts that says s is defined as a bool
local s = ls.snippet
local sn = ls.snippet_node
local i = ls.insert_node
local f = ls.f
local d = ls.d

function M.get_visual(args, parent)
    -- Summary: When `SELECT_RAW` is populated with a visual selection, the function
    -- returns an insert node whose initial text is set to the visual selection.
    -- When `SELECT_RAW` is empty, the function simply returns an empty insert node.
    if (#parent.snippet.env.SELECT_RAW > 0) then
        return sn(nil, i(1, parent.snippet.env.SELECT_RAW))
    else  -- If SELECT_RAW is empty, return a blank insert node
        return sn(nil, i(1))
    end
end

function M.regGroup(_, snip, i) return snip.captures[i] end
function M.re(i) return f(M.regGroup, nil, {user_args={i}}) end

function M.virt(text)
    -- usage example: c(1, {t("", virt("^l for alt choice")), t"alt choice"})
    return {node_ext_opts={passive={virt_text={{text, "@comment"}}}}}
end

function M.match_ahead(n)
    return function (line_to_cursor, trigger)
        -- look for match which ends n chars ahead of the cursor.
        -- Note that it doesn't replace the chars ahead, it just requires it in the pattern.
        -- modified version of default match_pattern
        -- https://github.com/L3MON4D3/LuaSnip/blob/master/lua/luasnip/nodes/util/trig_engines.lua
        local line = vim.api.nvim_get_current_line()
        local _, col = unpack(vim.api.nvim_win_get_cursor(0))
        local find_res = { line:sub(1, col+n):find(trigger .. "$") }

        if #find_res > 0 then
            -- if there is a match, determine matching string, and the
            -- capture-groups.
            local captures = {}
            -- find_res[1] is `from`, find_res[2] is `to` (which we already know
            -- anyway).
            local from = find_res[1]
            local match = line_to_cursor:sub(from, #line_to_cursor)
            -- collect capture-groups.
            for i = 3, #find_res do
                captures[i - 2] = find_res[i]
            end
            return match, captures
        else
            return nil
        end
    end
end

---Get functionnode that uppercases the text in node i
---@i int the index of another node
---@return functionnode upper
function M.upper(i)
    return f(function (args, parent, user_args)
        return args[1][1]:upper()
    end, {i})
end
---https://github.com/L3MON4D3/LuaSnip/blob/master/DOC.md#dynamicnode
local function rep_jump(args)
    -- the returned snippetNode doesn't need a position; it's inserted
    -- "inside" the dynamicNode.
    return sn(nil, {
        -- jump-indices are local to each snippetNode, so restart at 1.
        i(1, args[1][1])
    })
end
---Like the extras.rep node, but allows a jump so the repeated text can be changed.
---@param jump_index integer when the cursor will jump here in the progression through a snippet.
---@param node_reference integer jump_index of another node for which to copy the text.
---@return dynamicnode
function M.rep_jump(jump_index, node_reference)
    return d(jump_index, rep_jump, {node_reference})
end
local function upper_jump(args)
    -- the returned snippetNode doesn't need a position; it's inserted
    -- "inside" the dynamicNode.
    return sn(nil, {
        -- jump-indices are local to each snippetNode, so restart at 1.
        i(1, args[1][1]:upper())
    })
end
---Uppercase node with jump.
---@param jump_index integer when the cursor will jump here in the progression through a snippet.
---@param node_reference integer jump_index of another node for which to copy the text and uppercase it.
---@return dynamicnode
function M.upper_jump(jump_index, node_reference)
    return d(jump_index, upper_jump, {node_reference})
end

-- luasnip line_end condition doesn't seem to currently work with blink.cmp
if pcall(require, 'blink.cmp') then
    M.line_end = function () return true end
else
    M.line_end = extras.conds.line_end
end

return M
