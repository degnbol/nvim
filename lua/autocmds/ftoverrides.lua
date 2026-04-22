-- Overrides for nvim's bundled runtime ftplugins that set formatoptions in
-- ways we disagree with. $VIMRUNTIME/ftplugin/{julia,lua,vim}.vim do a
-- hard-coded `fo+=croql`, which re-adds `o` on top of our global default.
-- A user FileType autocmd runs after the built-in `filetypeplugin` callback
-- finishes its runtime! traversal, so this is a reliable late hook.

local grp = vim.api.nvim_create_augroup("ftoverrides", { clear = true })

-- `o` continues a comment when opening a new line with `o` / `O`. We prefer
-- only continuing while editing an existing comment line (`r`).
vim.api.nvim_create_autocmd("FileType", {
    group = grp,
    pattern = { "julia", "lua", "vim" },
    callback = function() vim.opt_local.formatoptions:remove("o") end,
})

-- $VIMRUNTIME/ftplugin/julia.vim sets commentstring to `# %s`; we prefer no
-- space so `gcc` produces `#foo` rather than `# foo`.
vim.api.nvim_create_autocmd("FileType", {
    group = grp,
    pattern = "julia",
    callback = function() vim.opt_local.commentstring = "#%s" end,
})
