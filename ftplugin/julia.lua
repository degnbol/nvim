-- remove o, we want to continue comments while editing them only (r).
vim.opt.formatoptions = "jwcrql"

-- see ../after/syntax/julia.vim
vim.opt.conceallevel = 1
-- you have to be in insert mode to unconceal
vim.opt.concealcursor = "nvc"

-- A fallback if indent is bad in Julia.
-- Simply use the builtin smartindent instead of indentexpr.
-- UNLESS the additional_vim_regex_highlighting option is set for treesitter, then use GetJuliaIndent.
-- vim.opt.indentexpr = "" -- is set in after/syntax/julia.vim since this dir isn't late enough
-- vim.opt.smartindent = true

vim.opt.list = false

require"completion/plotlyjs".setup()

