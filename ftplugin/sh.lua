-- textobj iP and aP is for between pipes in e.g. zsh by default but also map
-- i| and a|. Note that textobj "|" means column down, which is useful for e.g.
-- <C-v>|<C-q> where <C-q> is my custom blockim keymap.
vim.keymap.set({ "o", "x" }, "i|",
    "<cmd>lua require('various-textobjs').shellPipe(true)<CR>",
    { buffer = true }
)
vim.keymap.set({ "o", "x" }, "a|",
    "<cmd>lua require('various-textobjs').shellPipe(false)<CR>",
    { buffer = true }
)

vim.opt_local.list = false

-- stop continuing comments with o/O
vim.opt_local.formatoptions:remove('o')

-- shfmt for gq via formatprg (conform.nvim handles grf)
if vim.fn.executable("shfmt") == 1 then
    local dialect = vim.bo.filetype:find("zsh") and "zsh" or "auto"
    vim.opt_local.formatprg = "shfmt -ln " .. dialect .. " -i 4 -bn -ci -sr"
end
