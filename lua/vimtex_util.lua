#!/usr/bin/env lua
local M = {}

function M.in_math()
  return vim.fn['vimtex#syntax#in_mathzone']() == 1
end
function M.in_text() return not M.in_math() end
function M.in_comment()
  return vim.fn['vimtex#syntax#in_comment']() == 1
end
M.is_inside = vim.fn['vimtex#env#is_inside']
function M.in_env(name)
    local r, c = unpack(M.is_inside(name))
    return (r > 0 and c > 0)
end
function M.in_itemize() return M.in_env('itemize') end

-- Use vimtex (1,1)-indexing!
-- only "name" is required.
function M.inside_cmd(name, r, c)
    local cur
    if r == nil or c == nil then
        cur = vim.fn['vimtex#cmd#get_current']()
    else
        cur = vim.fn['vimtex#cmd#get_at'](r, c)
    end
    for i = 1, 20, 1 do
        if vim.tbl_isempty(cur) then return end
        local rc1, rc2 = cur['pos_start'], cur['pos_end']
        if cur['name'] == '\\' .. name then
            return {rc1['lnum'], rc1['cnum'], rc2['lnum'], rc2['cnum']}
        end
        -- minus 1 to go outside (before) the already looked at cmd
        cur = vim.fn['vimtex#cmd#get_at'](rc1['lnum'], rc1['cnum'])
    end
end

-- use vimtex (1,1)-indexing!
function M.inside_cmd_range(name, r1, c1, r2, c2)
    local cmds = {}
    local cmd1 = M.inside_cmd(name, r1, c1)
    local cmd2 = M.inside_cmd(name, r2, c2)
    if cmd1 ~= nil then table.insert(cmds, cmd1) end
    -- insert cmd2 last ...
    
    -- minus 1 on rows and columns for 0-index, not for c2 since end-exclusive.
    local text = vim.api.nvim_buf_get_text(0, r1-1, c1-1, r2-1, c2, {})
    for i, line in ipairs(text) do
        local j = 0
        while true do
            j = string.find(line, '\\' .. name, j+1, true)
            if j == nil then break end
            local cmd = vim.fn['vimtex#cmd#get_at'](r1+i-1, c1+j)
            local rc1 = cmd['pos_start']
            local rc2 = cmd['pos_end']
            if (cmd1 == nil or rc1['cnum'] ~= cmd1[2]) and (cmd2 == nil or rc1['cnum'] ~= cmd2[2]) then
                table.insert(cmds, {rc1['lnum'], rc1['cnum'], rc2['lnum'], rc2['cnum']})
            end
        end
    end
    
    -- insert cmd2
    -- make sure not to add the same cmd range twice
    if cmd2 ~= nil and (cmd1 == nil or (cmd2[1] ~= cmd1[1] or cmd2[2] ~= cmd1[2])) then
        table.insert(cmds, cmd2)
    end

    return cmds
end


return M
