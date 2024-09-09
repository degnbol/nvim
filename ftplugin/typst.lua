-- assume using typst.vim
vim.keymap.set('n', '<leader>cc', "<Cmd>TypstWatch<CR>", { desc="Compile continuously", buffer=true })
vim.keymap.set('n', '<leader>cv', "!open -a sioyek<CR>", { desc="Compile view", buffer=true, silent=true, })

-- toggle updating pdf on typing.
-- Doesn't work if running module file, hence we use TypstWatch instead by default.
vim.keymap.set('n', '<leader><leader>t', function ()
    local lsp = require"lspconfig".typst_lsp
    local exportPdf
    if lsp.manager.config.settings.exportPdf == "onType" then
        exportPdf = "never" -- use TypstWatch to compile
    else
        exportPdf = "onType"
    end
    vim.cmd.LspStop()
    require"lspconfig".typst_lsp.setup { settings = { exportPdf = exportPdf } }
    vim.cmd.LspStart()
    print(exportPdf)
end, { desc="Toggle update on type" })

vim.opt_local.wrap = true
vim.opt_local.sidescrolloff = 0
-- definitely don't break at @ in typst (when wrapping with linebreak). 
-- It is used right before citations and is not a place to break.
vim.opt_local.breakat:remove('@')
