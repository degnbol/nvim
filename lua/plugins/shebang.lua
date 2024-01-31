-- insert shebang on new file edit
return {
    "samirettali/shebang.nvim",
    config = function()
        vim.g.shebang_shells = {
            sh = 'zsh', -- non-standard use of .sh extension
            zsh = 'zsh',
            R = 'Rscript',
        }
    end
}
