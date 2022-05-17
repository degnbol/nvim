api = vim.api
cmd = vim.cmd
fn = vim.fn
-- If a new file is opened without extension, assume it is an incomplete tab-completion match.
-- We deal with it by opening the parent folder and moving the cursor to the first match.
-- Since we use BufNewFile this won't be a problem opening an already existing file named without extension.
-- If you really wish to make a file without extension then simply use `touch` first.
api.nvim_create_autocmd({"BufNewFile"}, {
    pattern="*",
    callback=function()
        if fn.expand('%:e') == '' then
            fname = fn.expand('%:t')
            folder = fn.expand('%:h')
            cmd 'bdel'
            cmd("e " .. folder)
            -- we need to let the tree explorer plugin open the window before we can search.
            vim.wait(1)
            cmd('silent! /' .. fname)
        end
    end
})
