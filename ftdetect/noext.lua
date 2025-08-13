-- If a new file is opened without extension, assume it is an incomplete
-- tab-completion match.
-- Since we use BufNewFile this won't be a problem opening an already existing
-- file named without extension.
-- If you really wish to make a file without extension then simply use `touch`
-- first or specify filename after opening nvim.
vim.api.nvim_create_autocmd("BufNewFile", {
    pattern = "*",
    callback = function()
        local path = vim.api.nvim_buf_get_name(0)
        local fname = vim.fs.basename(path)
        -- We look for any occurrence of a literal dot, i.e. a file like
        -- ".gitignore" is also fine.
        if not fname:contains('.') then
            vim.cmd.quit()
        end
    end
})
