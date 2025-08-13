require "autocmds/lastplace"
require "autocmds/templates"
require "autocmds/treesitter"
require "autocmds/chmodx"
require "autocmds/lsp"
require "autocmds/kitty"

vim.api.nvim_create_autocmd("Filetype", {
    pattern = "*",
    group = vim.api.nvim_create_augroup("cmp-syntax", { clear = true }),
    callback = function()
        if vim.opt_local.omnifunc:get() == "" then
            vim.opt_local.omnifunc = "syntaxcomplete#Complete"
        end
    end
})

-- Hide search highlight when changing mode.
vim.api.nvim_create_autocmd("ModeChanged", {
    -- Mode changed is indicated as <before>:<after> with the codes:
    -- i=insert, c=cmdline, n=normal, v=visual, V=line visual, \x16=block visual (guessing it means <C-V>), no=normal operator pending.
    -- We need to avoid responding to changes between normal and cmdline mode since that change is triggered frequently by plugins etc.
    -- These patterns should capture the same:
    -- pattern = "*:*[oivV\x16]*",
    pattern = "*:*[^nc]",
    group = vim.api.nvim_create_augroup("auto_hlsearch", { clear = true }),
    callback = function()
        vim.schedule(vim.cmd.nohlsearch)
    end
})
