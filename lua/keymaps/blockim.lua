local util = require "utils/init"
local map = require "utils/keymap"

-- blockim(proved)
-- Record a macro on one line, replay it on all other lines in the visual selection.
-- Uses RecordingLeave autocmd (press q to finish) instead of dynamic Esc mapping.
map.cv('<C-q>', function()
    local rVis1, _, rVis2, _ = util.get_visual_range()
    local rCur1, cCur1 = unpack(vim.api.nvim_win_get_cursor(0))
    vim.api.nvim_create_autocmd("RecordingLeave", {
        once = true,
        callback = function()
            vim.schedule(function()
                local reg = vim.fn.getreg("b")
                if reg == "" then return end
                local rCur2, cCur2 = unpack(vim.api.nvim_win_get_cursor(0))
                for i = rVis1, rVis2 do
                    if i ~= rCur1 then
                        vim.api.nvim_win_set_cursor(0, { i, cCur1 })
                        -- group all replays into one undo entry with the original edit.
                        vim.cmd "undojoin"
                        -- feedkeys with "nx": noremap + execute immediately.
                        -- norm @b doesn't work in vim.schedule after RecordingLeave.
                        vim.fn.feedkeys(reg, "nx")
                    end
                end
                vim.api.nvim_win_set_cursor(0, { rCur2, cCur2 })
            end)
        end,
    })
    vim.cmd.norm "qb"
end, "Block mode improved. Apply macro along block selection.")
