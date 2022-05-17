-- julia-vim sets it to #= %s =# by default https://github.com/JuliaEditorSupport/julia-vim/blob/master/ftplugin/julia.vim
vim.api.nvim_buf_set_option(0, "commentstring", "# %s")
-- remove o, we want to continue comments while editing them only (r).
vim.opt.formatoptions = "tjwcrqla"
