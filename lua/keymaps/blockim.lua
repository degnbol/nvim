local util = require "utils/init"
local map = require "utils/keymap"

-- blockim(proved)
map.cv('<C-q>', function()
    local rVis1, cVis1, rVis2, cVis2 = util.get_visual_range()
    local rCur1, cCur1 = unpack(vim.api.nvim_win_get_cursor(0))
    vim.keymap.set('n', '<Esc>', function()
        vim.keymap.del('n', '<Esc>')
        vim.cmd.stopinsert()
        vim.cmd.norm "q"
        -- if nothing was recorded then just exit, otherwise we get an undojoin related error.
        if not vim.fn.getreg("b") then return end
        local rCur2, cCur2 = unpack(vim.api.nvim_win_get_cursor(0))
        for i = rVis1, rVis2 do
            if i ~= rCur1 then
                vim.api.nvim_win_set_cursor(0, { i, cCur1 })
                -- group the macro eval to the one before,
                -- and all the way back to the initial macro recording, so we
                -- can undo the whole multi line edit with one u.
                vim.cmd "undojoin"
                vim.cmd.norm "@b"
            end
        end
        vim.cmd.stopinsert()
        vim.api.nvim_win_set_cursor(0, { rCur2, cCur2 })
    end, { desc = "Run block mode macro on each selected line at the column of the cursor." })
    vim.cmd.norm "qb"
end, "Block mode improved. Apply macro along block selection.")
