-- If a new file is opened without extension, assume it is an incomplete 
-- tab-completion match. We deal with it by opening the parent folder and 
-- moving the cursor to the first match. Since we use BufNewFile this won't be 
-- a problem opening an already existing file named without extension.
-- If you really wish to make a file without extension then simply use `touch` 
-- first or specify filename after opening nvim.
vim.api.nvim_create_autocmd("BufNewFile", {
    pattern="*", callback=function()
        local ext = vim.fn.expand('%:e')
        local fname = vim.fn.expand('%:t')
        -- make exception for filename that starts with '.', e.g. .gitignore, 
        -- since this is not a mistake even though it doesn't have an 
        -- extension.
        if ext == '' and fname:sub(1, 1) ~= '.' then
            vim.cmd "q"
        end
    end
})
