-- Runs after vim-asciidoc's ftplugin (which clears 'comments' and provides
-- the Asciidoctor2PDF command + asciidoctor2pdf compiler). Wrapped in pcall
-- because the compiler file isn't available until the plugin has been
-- packadded — safe to fall through on the first load before lz.n triggers.

vim.opt_local.comments = "://"

pcall(vim.cmd, [[compiler asciidoctor2pdf]])

-- double <CR> to auto-close after successful compilation. Use :make to keep
-- the error list open.
vim.keymap.set("n", "<leader>cc", "<Cmd>Asciidoctor2PDF<CR><CR>",
    { buffer = true, desc = "Compile to PDF" })
vim.keymap.set("n", "<leader>oo", "<Cmd>AsciidoctorOpenPDF<CR><CR>",
    { buffer = true, desc = "Open compiled PDF" })
vim.keymap.set("n", "<leader>cC", function()
    vim.api.nvim_create_autocmd("BufWritePost", {
        buffer = 0,
        group = vim.api.nvim_create_augroup("asciidocCompile", { clear = true }),
        command = "silent Asciidoctor2PDF",
    })
end, { buffer = true, desc = "Compile on save" })
