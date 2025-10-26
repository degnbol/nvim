-- remove o, we want to continue comments while editing them only (r).
vim.opt_local.formatoptions = "jwcrql"

-- you have to be in insert mode to unconceal
vim.opt_local.concealcursor = "nvc"

vim.opt_local.list = false

-- using blink.cmp instead
-- require"completion.plotlyjs.cmp_plotlyjs".setup()

-- manually mark that plotlyjs is being used
vim.keymap.set('n', '<localleader>+', function()
    vim.g.loaded_plotly = true
end, { buffer = true, desc = "Manually load plotlyjs completion" })
-- ... or check if plotlyjs is loaded by scanning first 20 lines
for _, line in ipairs(vim.api.nvim_buf_get_lines(0, 0, 20, false)) do
    if line:lower():match("plotly") then
        vim.g.loaded_plotly = true
    end
end

-- Manual efforts. Install julia LSP as described on
-- https://github.com/neovim/nvim-lspconfig/blob/master/lsp/julials.lua
vim.lsp.enable("julials")
-- vim.lsp.enable("jetls")

-- Annoying bug where something keeps setting buftype=nofile
vim.api.nvim_create_autocmd("BufWritePre", {
    buffer = 0,
    group = vim.api.nvim_create_augroup("julia_bug", { clear = true }),
    callback = function() vim.bo.buftype = nil end
})

-- TEMP: there should be a better way to do this. syn must be getting reset somewhere.
vim.defer_fn(function ()
    vim.cmd [[syn keyword @variable.builtin stdin stdout stderr]]
end, 0)

