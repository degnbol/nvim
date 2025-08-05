local grp = vim.api.nvim_create_augroup("chmodx", { clear = true })

-- Auto chmod u+x some filetypes.
vim.api.nvim_create_autocmd("BufNewFile", {
    pattern = { "*.sh", "*.zsh", "*.r", "*.jl" },
    group = grp,
    -- Automatically do chmod u+x for a new file if it gets written.
    callback = function()
        vim.api.nvim_create_autocmd("BufWritePost", {
            buffer = 0,
            once = true,
            group = grp,
            callback = function()
                vim.system({ "chmod", "u+x", vim.api.nvim_buf_get_name(0) }):wait()
            end
        })
    end
})
