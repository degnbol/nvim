local map = require "utils/keymap"
local util = require "utils/init"

map.n("gcA", function()
    local comment = vim.split(
        vim.bo.commentstring,
        "%s", { plain = true }
    )
    local r, _ = unpack(vim.api.nvim_win_get_cursor(0))
    local l_line = #vim.api.nvim_get_current_line()
    vim.api.nvim_win_set_cursor(0, { r, l_line })
    vim.api.nvim_put({ ' ' .. comment[1] }, 'c', true, true)
    vim.api.nvim_put({ comment[2] }, 'c', true, false)
    vim.cmd.startinsert()
    vim.api.nvim_win_set_cursor(0, { r, l_line + #comment[1] + 1 })
end, "New comment at EOL")

map.n("gco", function()
    local comment = vim.split(
        vim.bo.commentstring,
        "%s", { plain = true }
    )
    vim.cmd.normal 'o'
    vim.api.nvim_put({ comment[1] }, 'c', true, true)
    vim.api.nvim_put({ comment[2] }, 'c', true, false)
    vim.cmd.startinsert()
    util.set_col(#comment[1])
end, "New comment below")

map.n("gcO", function()
    local comment = vim.split(
        vim.bo.commentstring,
        "%s", { plain = true }
    )
    vim.cmd.normal 'O'
    vim.api.nvim_put({ comment[1] }, 'c', true, true)
    vim.api.nvim_put({ comment[2] }, 'c', true, false)
    vim.cmd.startinsert()
    util.set_col(#comment[1])
end, "New comment above")
