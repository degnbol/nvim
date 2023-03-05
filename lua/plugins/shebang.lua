-- insert shebang on new file edit
return {"samirettali/shebang.nvim", config=function()
    -- non-standard use of .sh extension
    vim.g.shebang_commands = {
        sh = '/usr/bin/env zsh',
        zsh = '/usr/bin/env zsh',
        R = '/usr/bin/env Rscript'
    }
end
}
