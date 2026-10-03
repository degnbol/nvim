-- textobj iP and aP is for between pipes in e.g. zsh by default but also map
-- i| and a|. Note that textobj "|" means column down, which is useful for e.g.
-- <C-v>|<C-q> where <C-q> is my custom blockim keymap.
vim.keymap.set({ "o", "x" }, "i|",
    "<cmd>lua require('various-textobjs').shellPipe(true)<CR>",
    { buf = 0 }
)
vim.keymap.set({ "o", "x" }, "a|",
    "<cmd>lua require('various-textobjs').shellPipe(false)<CR>",
    { buf = 0 }
)

vim.opt_local.list = false

-- stop continuing comments with o/O
vim.opt_local.formatoptions:remove('o')

-- shfmt for gq via formatprg (conform.nvim handles grf)
if vim.fn.executable("shfmt") == 1 then
    local is_zsh = vim.bo.filetype:find("zsh") ~= nil
    local prg = "shfmt -ln " .. (is_zsh and "zsh" or "auto") .. " -i 4 -ci -sr"
    -- Chain the cosmetic mlr-invocation formatter for zsh (order-independent).
    -- A filter's stderr lands in the buffer ('shellredir'), so fmt-zsh-mlr's
    -- warnings go to a file and are notified after the filter.
    if is_zsh and vim.fn.executable("fmt-zsh-mlr") == 1 then
        local stderr_file = vim.fn.tempname()
        prg = prg .. " | fmt-zsh-mlr 2>" .. vim.fn.shellescape(stderr_file)
        local group = vim.api.nvim_create_augroup("fmt_zsh_mlr_stderr", { clear = false })
        vim.api.nvim_clear_autocmds({ group = group, buf = 0 })
        vim.api.nvim_create_autocmd("ShellFilterPost", {
            group = group,
            buf = 0,
            callback = function()
                local stderr = vim.trim(require("utils/init").readtext(stderr_file) or "")
                os.remove(stderr_file)
                if stderr ~= "" then vim.notify(stderr, vim.log.levels.WARN) end
            end,
        })
    end
    vim.opt_local.formatprg = prg
end
