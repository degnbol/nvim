api = vim.api
cmd = vim.cmd
fn = vim.fn
-- If a new file is opened ending in . it is a mistake when auto completing and
-- There are multiple extension matches.
-- We deal with it by opening the parent folder and moving the cursor to the first match.
api.nvim_create_autocmd({"BufNewFile"}, {
    pattern={"*."},
    callback=function()
        fname = fn.expand('%:t')
        folder = fn.expand('%:h')
        cmd 'bdel'
        cmd("e " .. folder)
        -- we need to let the tree explorer plugin open the window before we can search.
        vim.wait(1)
        cmd('/' .. fname)
    end
})
