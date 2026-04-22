-- Overrides for plugin-provided ftplugin settings that load after ours.
-- A user FileType autocmd runs after the built-in `filetypeplugin` callback
-- (which does the `runtime! ftplugin/*.{vim,lua}` for both main and after/
-- directories), so this is a reliable late hook.

local grp = vim.api.nvim_create_augroup("ftoverrides", { clear = true })

-- Some plugins (julia-vim, etc.) set `fo+=o` in their ftplugin despite our
-- preference of only continuing comments while editing them (r), not when
-- opening a new line below with `o`.
vim.api.nvim_create_autocmd("FileType", {
    group = grp,
    pattern = { "julia", "lua", "vim", "zsh" },
    callback = function() vim.opt_local.formatoptions:remove("o") end,
})

-- asciidoc plugin sets commentstring but leaves 'comments' unset.
vim.api.nvim_create_autocmd("FileType", {
    group = grp,
    pattern = "asciidoc",
    callback = function() vim.opt_local.comments = "://" end,
})

-- julia-vim defaults commentstring to `#=%s=#` — we prefer line comments.
-- https://github.com/JuliaEditorSupport/julia-vim/blob/master/ftplugin/julia.vim
vim.api.nvim_create_autocmd("FileType", {
    group = grp,
    pattern = "julia",
    callback = function() vim.opt_local.commentstring = "#%s" end,
})
