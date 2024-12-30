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

-- using blink.cmp instead
-- require"completion.plotlyjs.cmp_plotlyjs".setup()

-- manually mark that plotlyjs is being used
vim.keymap.set('n', '<leader><leader>+', function ()
    vim.g.loaded_plotly = true
end, { desc="Manually load plotlyjs completion" })
-- ... or check if plotlyjs is loaded by scanning first 20 lines
for _, line in ipairs(vim.api.nvim_buf_get_lines(0, 0, 20, false)) do
    if line:lower():match("plotly") then
        vim.g.loaded_plotly = true
    end
end
