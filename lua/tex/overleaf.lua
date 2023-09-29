#!/usr/bin/env lua
local rtp = vim.opt.runtimepath:get()[1]

-- async overleaf sync on save and load
local overleaf = vim.api.nvim_create_augroup("overleaf", {clear=true})
local remote = vim.fn.system("git remote get-url origin")
if string.sub(remote, 1, 25) == "https://git.overleaf.com/" then
    -- When collaborating you probably don't want to insert linebreaks all over the place.
    vim.opt.wrap = true
    vim.opt.formatoptions=vim.opt.formatoptions- "a" -- no autoformat
    vim.api.nvim_del_augroup_by_name("foHack") -- don't switch fo+=a back on later
    -- print("Adding overleaf aucmd")
    vim.api.nvim_create_autocmd({"BufRead", "BufWritePost"}, {
        group = overleaf,
        buffer = 0, -- only for current buffer. Mutually exclusive with pattern arg.
        callback = function ()
            vim.fn.jobstart({rtp .. "/lua/tex/overleaf/gitsync.sh", vim.api.nvim_buf_get_name(0)})
        end
    })
end


