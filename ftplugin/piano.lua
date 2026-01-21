

vim.schedule(function ()
    vim.cmd "syn match Comment /b #   # b #   b #   # b #   b #   # b #/"
end)

local keys = {"a", "A", "a", "B", " ", "C", "d", "D", "d", "E", " ", "F", "f", "G"}

local function putKey()
    local r, c = unpack(vim.api.nvim_win_get_cursor(0))
    local offset = 7 -- offset that depends on where the piano starts from
    local key = keys[(c+offset) % #keys]
    vim.api.nvim_buf_set_text(0, r-1, c, r-1, c+1, {key})
end

vim.keymap.set('n', '<C-k>', putKey, { desc="Put piano key at cursor" })


