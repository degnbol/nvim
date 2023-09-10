#!/usr/bin/env lua

-- Changes the desc of a keymap. Useful e.g. to set desc for keymaps set with 
-- vimscript or by plugins that doesn't provide a desc.
-- Taken from set_mapping_desc from
-- https://github.com/echasnovski/mini.clue/blob/main/lua/mini/clue.lua
-- Modified to not throw errors, e.g. if mapping doesn't exist.
-- Instead, return whether it was successful.
function set_keymap_desc(mode, lhs, desc)
    -- allow mode to be either e.g. 'n' or {'n', 'x'}
    if type(mode) ~= 'string' then
        local ret = true
        for _, m in ipairs(mode) do
            ret = ret and set_keymap_desc(m, lhs, desc)
        end
        return ret
    else
        local map_data = vim.fn.maparg(lhs, mode, false, true)
        if vim.tbl_count(map_data) == 0 then return false end
        map_data.desc = desc
        vim.fn.mapset(mode, false, map_data)
        return true
    end
end
