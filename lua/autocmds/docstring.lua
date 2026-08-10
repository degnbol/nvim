-- Register `:Docstring` (reconcile the cursor function's docstring against its
-- signature) as a buffer-local command on the filetypes the docstring library
-- supports. Buffer-local so the command only exists where it can act; the
-- implementation is required lazily on first FileType match.

vim.api.nvim_create_autocmd("FileType", {
    group = vim.api.nvim_create_augroup("my.docstring", { clear = true }),
    pattern = "python",
    callback = function(args)
        vim.api.nvim_buf_create_user_command(args.buf, "Docstring", function()
            require("docstring.command").reconcile_at_cursor()
        end, { desc = "Reconcile the cursor function's docstring with its signature" })
    end,
})
