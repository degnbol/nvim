--- Hacky solution to set a window offset in a buffer so that a given integer 
--- row in the buffer appears the as topline in the window.
local function setWinTopline(win, topline)
    scrolloff = vim.opt_local.scrolloff:get()
    height = vim.api.nvim_win_get_height(win)
    lastline = vim.fn.line "$"
    vim.api.nvim_win_set_cursor(win, {math.min(topline+height-scrolloff, lastline), 0})
    vim.api.nvim_win_set_cursor(win, {math.min(topline+scrolloff, lastline), 0})
end

function vnewFurther()
    for i = 1, math.max(vim.v.count,1) do
        win = vim.api.nvim_get_current_win()
        vim.opt_local.scrollbind = true
        topline = vim.fn.line("w0")
        height = vim.api.nvim_win_get_height(win)
        vim.cmd.vnew "%"
        setWinTopline(0, topline+height)
        vim.opt_local.scrollbind = true
    end
end

vim.keymap.set('n', '<leader><leader>v', vnewFurther, { desc="vnew with one window height offset and scrollbind enabled. Accepts a count." })

