-- Theme switching commands and startup colorscheme.
-- The generated colorscheme (~/nvim/colors/generated.lua) is a hardlink
-- to ~/dotfiles/colors/generated/<theme>/colorscheme.lua, managed by switch.sh.

vim.api.nvim_create_user_command("Dark", function()
    vim.fn.system("~/dotfiles/dark.sh")
    vim.cmd("colorscheme generated")
end, {})

vim.api.nvim_create_user_command("Light", function()
    vim.fn.system("~/dotfiles/light.sh")
    vim.cmd("colorscheme generated")
end, {})

vim.api.nvim_create_autocmd("VimEnter", {
    callback = function()
        vim.schedule(function()
            vim.cmd("colorscheme generated")
        end)
    end,
})
