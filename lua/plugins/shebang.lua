-- insert shebang on new file edit
return {
    "samirettali/shebang.nvim",
    event = "InsertEnter",
    config=function()
        vim.g.shebang_commands = {
            sh = '/usr/bin/env zsh', -- non-standard use of .sh extension
            zsh = '/usr/bin/env zsh',
            R = '/usr/bin/env Rscript'
        }
    end
}
