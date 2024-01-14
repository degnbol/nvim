local hl = require "utils/highlights"

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

hl.def("typstMarkupHeading", "@text.title")
-- Semantic tokens.
-- captures \ at the end of lines, ~ between words and \* that let's you write a literal *. The last one is not ideal but worth it for the other two.
hl.def("@lsp.type.escape", "@comment")
-- undo coloring from above on "..." with help from custom capture in after/queries/typst/highlights.scm
hl.def("@punct.ellipsis", "Normal")
hl.def("@lsp.typemod.punct.emph", "@comment")
hl.def("@lsp.typemod.punct.strong", "@comment")
hl.def("@lsp.type.punct", "Delimiter")
hl.def("@lsp.type.pol", "@variable")
hl.def("@lsp.type.string", "@string")
hl.def("@lsp.type.heading", "@text.title")
hl.def("@lsp.type.number", "@number")

vim.opt_local.wrap = true
vim.opt_local.sidescrolloff = 0

