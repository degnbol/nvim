-- In-process LSP providing hover for glossy `@term` refs in typst. Thin
-- protocol shim; all resolution lives in lua/typst_glossary.lua (shared with
-- the `grd` goto-def path). tinymist resolves glossy refs to the package's
-- runtime `label()` placeholder and returns no hover for them, so this joins
-- the standard multi-client hover (`vim.lsp.buf.hover`) as the sole provider on
-- glossary keys.

return require("utils.lsp").hover_server {
    filetypes = { "typst" },
    hover = function(bufnr, row, col)
        return require("typst_glossary").hover(bufnr, row, col)
    end,
}
