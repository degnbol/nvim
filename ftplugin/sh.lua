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
    local is_zsh = vim.bo.filetype:find("zsh") ~= nil
    local prg = "shfmt -ln " .. (is_zsh and "zsh" or "auto") .. " -i 4 -ci -sr"
    -- Chain the cosmetic mlr-invocation formatter for zsh (order-independent).
    if is_zsh and vim.fn.executable("fmt-zsh-mlr") == 1 then
        prg = prg .. " | fmt-zsh-mlr"
    end
    vim.opt_local.formatprg = prg
end
