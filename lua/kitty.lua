local grp = vim.api.nvim_create_augroup("KittySetVar", { clear = true })
vim.api.nvim_create_autocmd({ "VimEnter", "VimResume" }, {
    group = grp,
    callback = function()
        io.stdout:write("\x1b]1337;SetUserVar=nvim=MQo\007")
    end,
})

vim.api.nvim_create_autocmd({ "VimLeave", "VimSuspend" }, {
    group = grp,
    callback = function()
        io.stdout:write("\x1b]1337;SetUserVar=nvim\007")
    end,
})
