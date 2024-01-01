-- assume using typst.vim
vim.keymap.set('n', '<leader>cc', "<Cmd>TypstWatch<CR>", { desc="Compile continuously", buffer=true })
vim.keymap.set('n', '<leader>cv', "!open -a sioyek<CR>", { desc="Compile view", buffer=true, silent=true, })

-- toggle updating pdf on typing.
-- Doesn't work if running module file, hence we use TypstWatch instead by default.
vim.keymap.set('n', '<leader><leader>t', function ()
    local lsp = require"lspconfig".typst_lsp
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

vim.api.nvim_set_hl(0, "typstMarkupHeading", {link="@text.title", default=true})
-- Semantic tokens.
-- captures \ at the end of lines, ~ between words and \* that let's you write a literal *. The last one is not ideal but worth it for the other two.
vim.api.nvim_set_hl(0, "@lsp.type.escape.typst", {link="@comment", default=true})
-- undo coloring from above on "..." with help from custom capture in after/queries/typst/highlights.scm
vim.api.nvim_set_hl(0, "@punct.ellipsis.typst", {link="Normal", default=true})
vim.api.nvim_set_hl(0, "@lsp.typemod.punct.emph.typst", {link="@comment", default=true})
vim.api.nvim_set_hl(0, "@lsp.typemod.punct.strong.typst", {link="@comment", default=true})
vim.api.nvim_set_hl(0, "@lsp.type.punct.typst", {link="Delimiter", default=true})
vim.api.nvim_set_hl(0, "@lsp.type.pol.typst", {link="@variable", default=true})
vim.api.nvim_set_hl(0, "@lsp.type.string.typst", {link="@string", default=true})

-- local grp = vim.api.nvim_create_augroup("typst", {clear=true})
-- vim.api.nvim_create_autocmd("Colorscheme", {
--     buffer = true,
--     group = grp,
--     callback = function ()
--        
--     end
-- })

vim.opt_local.wrap = true
vim.opt_local.sidescrolloff = 0

