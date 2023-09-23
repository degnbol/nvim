local util = require "utils"

local function isblock()
    local mode = vim.fn.mode()
    return mode ~= "v" and mode ~= "V"
end

vim.keymap.set('v', '<C-q>', function ()
    if not isblock() then return end
    util.end_visual()
    local rVis1, cVis1, rVis2, cVis2 = util.get_visual_range()
    local rCur1, cCur1 = unpack(vim.api.nvim_win_get_cursor(0))
    vim.keymap.set({'n', 'i'}, '<Esc>', function ()
        vim.keymap.del({'n', 'i'}, '<Esc>')
        vim.cmd.stopinsert()
        vim.cmd.norm "q"
        local rCur2, cCur2 = unpack(vim.api.nvim_win_get_cursor(0))
        for i = rVis1, rVis2 do
            if i ~= rCur1 then
                vim.api.nvim_win_set_cursor(0, {i, cCur1})
                -- group the macro eval to the one before, 
                -- and all the way back to the initial macro recording, so we 
                -- can undo the whole multi line edit with one u.
                vim.cmd "undojoin"
                vim.cmd.norm "@b"
            end
        end
        vim.cmd.stopinsert()
        vim.api.nvim_win_set_cursor(0, {rCur2, cCur2})
    end, { desc="Run macro." })
    vim.cmd.norm "qb"
end, { desc="Block mode improved. Apply macro along block selection." })

