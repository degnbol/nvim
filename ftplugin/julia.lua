-- julia-vim sets it to #= %s =# by default https://github.com/JuliaEditorSupport/julia-vim/blob/master/ftplugin/julia.vim
vim.api.nvim_buf_set_option(0, "commentstring", "# %s")
-- remove o, we want to continue comments while editing them only (r).
-- no t and having c+a means only comments are autoformatted.
vim.opt.formatoptions = "jwcrqla"
-- when treesitter runs for julia the indent is wrong after inline for loop,
-- but setting the indentexpr to the one from treesitter it works, but then it
-- breaks indents after end of function.
-- The problem with indent happening when writing a bracket is actually due to
-- indentkeys which are a set of keys that when written will trigger autoindent
-- of a line. Setting the indent of a line is ALWAYS done by calling
-- indentexpr which is GetJuliaIndent() by default. 
-- Setting the treesitter config enable flag sets indentexpr=nvim_treesitter#indent()
-- I think I finally found the solution to fix indent for julia.
-- Simply use the builtin smartindent instead of indentexpr.
-- UNLESS the additional_vim_regex_highlighting option is set for treesitter, then use GetJuliaIndent.
-- vim.opt.indentexpr = "" -- is set in after/syntax/julia.vim since this dir isn't late enough
-- vim.opt.smartindent = true

