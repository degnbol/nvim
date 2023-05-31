#!/usr/bin/env lua
local util = require "util"
local before = util.before
local end_visual = util.end_visual
local get_visual_range = util.get_visual_range
local vtu = require "vimtex_util"
local in_math = vtu.in_math
local cmd = vim.cmd

-- bold and italics. Others could be underline, color text etc.
local textbf = "textbf"
local textit = "textit"
-- unlike mathbf that only bolds latin, symbf with unicode-math works with greek.
local mathbf = "symbf"
-- math text is italic by default, so we toggle roman instead.
-- Unlike mathrm that only un-italizes latin, symrm with unicode-math works with greek.
-- There's also \text but it does a bit more than just un-italize. It 
-- also makes text a bit less strong, so it should maybe be an explicit 
-- choice.
local mathit = "symrm"

local insert_or_del_cmd = function(name)
    cmd = vtu.inside_cmd(name)
    if cmd ~= nil then
        vim.fn['vimtex#cmd#delete'](cmd[1], cmd[2])
    else
        -- are we inside a cmd?
        local line = vim.api.nvim_get_current_line()
        local r, c = unpack(vim.api.nvim_win_get_cursor(0))
        -- +1 to include the char right after cursor bar
        local cmdStart = line:sub(1, c+1):match('\\%a+$')
        if cmdStart ~= nil then c=c-#cmdStart+1 end -- +1 to undo the +1
        vim.api.nvim_buf_set_text(0, r-1, c, r-1, c, {'\\' .. name .. '{'})
        if cmdStart == nil then
            -- move cursor to after the insertion
            vim.api.nvim_win_set_cursor(0, {r, c+ #name + 2})
        end
    end
end

-- TODO: if range start or end in any command, expand the range to include the command. Just like above
-- surround a selection with a cmd that is joined if overlapping, e.g. textbf 
-- where it isn't meaningful to bold text twice.
local function surround_visual(name)
    end_visual()
    local r1, c1, r2, c2 = get_visual_range()
    
    -- if there are unbalanced {} then we will break things
    local text = table.concat(vim.api.nvim_buf_get_text(0, r1-1, c1, r2-1, c2+1, {}), '\n')
    text = text:gsub('%b{}', '')
    text = text:gsub('\\' .. name .. '{', '') -- not a problem if it is the surround cmd
    local _, brs = text:gsub("[{}]", "") -- count
    if brs > 0 then
        -- TODO: add not implemented warning
        cmd "normal gv" -- un-cancel select
        return
    end
    
    c1, c2 = c1+1, c2+1 -- use vimtex (1,1)-index
    -- find cmds already present.
    -- We expand the "search" by 1 so that we consume (merge with) adjacent cmds as well.
    -- the max on c1 is needed. It would be easy enough to search wrapping lines as well, e.g. if c1==1 then r1-=1, c1=EOL end
    local cmds = vtu.inside_cmd_range(name, r1, math.max(1,c1-1), r2, c2+1)
    -- if any go beyond the selection we expand the selection to merge them
    for _, cmd in ipairs(cmds) do
        _r1, _c1, _r2, _c2 = unpack(cmd)
        if before(_r1, _c1, r1, c1) then
            r1, c1 = _r1, _c1
        end
        if before(r2, c2, _r2, _c2) then
            r2, c2 = _r2, _c2
        end
    end
    -- easier to modify and replace the selection text rather than the buffer 
    -- when changes might affect indexes for other cmds.
    -- -1 for 0-index, end-exclusive so not on r2
    local lines = vim.api.nvim_buf_get_lines(0, r1-1, r2, true)
    local shift, r_last = 0, r1 -- keep track of the shift in the next indexes
    for iCmd, cmd in ipairs(cmds) do
        _r1, _c1, _r2, _c2 = unpack(cmd)
        if _r1 ~= r_last then shift = 0 else -- reset on new line
            _c1=_c1-shift
        end
        -- remove beginning
        lines[_r1-r1+1] = lines[_r1-r1+1]:sub(1, _c1-1) .. lines[_r1-r1+1]:sub(_c1 + #name + 2)
        if _r2 ~= _r1 then shift = 0 else
            shift=shift+#name+2
            _c2=_c2-shift
        end
        -- remove end
        lines[_r2-r1+1] = lines[_r2-r1+1]:sub(1, _c2-1) .. lines[_r2-r1+1]:sub(_c2+1)
        -- keep note for next iteration
        shift=shift+1
        r_last = _r2
    end
    -- then surround the whole selection with the cmd
    lines[1] = lines[1]:sub(1, c1-1) .. '\\' .. name .. '{' .. lines[1]:sub(c1)
    -- if we are adding to a single line then the end will be shifted opposite
    if #lines == 1 then c2=c2 - shift + #name + 2
    elseif r2 == r_last then -- we may be multiple lines but have already modified the last line
        c2=c2-shift
    end
    lines[#lines] = lines[#lines]:sub(1, c2) .. '}' .. lines[#lines]:sub(c2+1)
    vim.api.nvim_buf_set_lines(0, r1-1, r2, true, lines)
end

vim.keymap.set("i", "<C-b>", function ()
    if in_math() then
        insert_or_del_cmd(mathbf)
    else
        insert_or_del_cmd(textbf)
    end
end)
vim.keymap.set("i", "<C-i>", function ()
    if in_math() then
        insert_or_del_cmd(mathit)
    else
        insert_or_del_cmd(textit)
    end
end)
vim.keymap.set("x", "<C-b>", function () 
    if in_math() then
        surround_visual(mathbf)
    else
        surround_visual(textbf)
    end
end)
vim.keymap.set("x", "<C-i>", function ()
    if in_math() then
        surround_visual(mathit)
    else
        surround_visual(textit)
    end
end)

