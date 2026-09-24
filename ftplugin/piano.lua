local util = require "utils/init"

local keys = {"a", "A", "a", "B", " ", "C", "d", "D", "d", "E", " ", "F", "f", "G"}

local function putKey()
    local r, c = util.get_cursor()
    local offset = 7 -- offset that depends on where the piano starts from
    local key = keys[(c+offset) % #keys]
    vim.api.nvim_buf_set_text(0, r, c, r, c+1, {key})
end

vim.keymap.set('n', '<C-k>', putKey, { buffer=true, desc="Put piano key at cursor" })


