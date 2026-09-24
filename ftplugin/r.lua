vim.cmd([[syntax region @variable matchgroup=Operator start=/`/ end=/`/]])

-- Errors.
vim.diagnostic.enable(false, { bufnr = 0 })

-- Fallback for <leader>cc when the script has no shebang.
vim.b.interpreter = 'Rscript'
