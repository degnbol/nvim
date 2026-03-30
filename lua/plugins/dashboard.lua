local grp = vim.api.nvim_create_augroup("stdin_read", { clear = true })
vim.api.nvim_create_autocmd("StdinReadPre", {
    pattern = "*",
    group = grp,
    callback = function()
        -- Set this to avoid dashboard opening on a pipe without dash, i.e. `echo hello | nvim`
        -- Relevant discussions
        -- https://github.com/nvimdev/dashboard-nvim/pull/513
        -- https://github.com/nvimdev/dashboard-nvim/issues/443
        vim.g.read_from_stdin = true
    end
})

-- open to a dashboard for vi without a file selection, requires telescope or an alternative installed.
return {
    "dashboard-nvim",
    event = 'VimEnter',
    after = function()
        require("dashboard").setup({
            change_to_vcs_root = true,
            config = {
                header = {
                    '⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣀⣠⣤⣤⣴⣦⣤⣤⣄⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀ ',
                    '⠀⠀⠀⠀⠀⠀⢀⣤⣾⣿⣿⣿⣿⠿⠿⠿⠿⣿⣿⣿⣿⣶⣤⡀⠀⠀⠀⠀⠀⠀ ',
                    '⠀⠀⠀⠀⣠⣾⣿⣿⡿⠛⠉⠀⠀⠀⠀⠀⠀⠀⠀⠉⠛⢿⣿⣿⣶⡀⠀⠀⠀⠀ ',
                    '⠀⠀⠀⣴⣿⣿⠟⠁⠀⠀⠀⣶⣶⣶⣶⡆⠀⠀⠀⠀⠀⠀⠈⠻⣿⣿⣦⠀⠀⠀ ',
                    '⠀⠀⣼⣿⣿⠋⠀⠀⠀⠀⠀⠛⠛⢻⣿⣿⡀⠀⠀⠀⠀⠀⠀⠀⠙⣿⣿⣧⠀⠀ ',
                    '⠀⢸⣿⣿⠃⠀⠀⠀⠀⠀⠀⠀⠀⢀⣿⣿⣷⠀⠀⠀⠀⠀⠀⠀⠀⠸⣿⣿⡇⠀ ',
                    '⠀⣿⣿⡿⠀⠀⠀⠀⠀⠀⠀⠀⢀⣾⣿⣿⣿⣇⠀⠀⠀⠀⠀⠀⠀⠀⣿⣿⣿⠀ ',
                    '⠀⣿⣿⡇⠀⠀⠀⠀⠀⠀⠀⢠⣿⣿⡟⢹⣿⣿⡆⠀⠀⠀⠀⠀⠀⠀⣹⣿⣿⠀ ',
                    '⠀⣿⣿⣷⠀⠀⠀⠀⠀⠀⣰⣿⣿⠏⠀⠀⢻⣿⣿⡄⠀⠀⠀⠀⠀⠀⣿⣿⡿⠀ ',
                    '⠀⢸⣿⣿⡆⠀⠀⠀⠀⣴⣿⡿⠃⠀⠀⠀⠈⢿⣿⣷⣤⣤⡆⠀⠀⣰⣿⣿⠇⠀ ',
                    '⠀⠀⢻⣿⣿⣄⠀⠀⠾⠿⠿⠁⠀⠀⠀⠀⠀⠘⣿⣿⡿⠿⠛⠀⣰⣿⣿⡟⠀⠀ ',
                    '⠀⠀⠀⠻⣿⣿⣧⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣠⣾⣿⣿⠏⠀⠀⠀ ',
                    '⠀⠀⠀⠀⠈⠻⣿⣿⣷⣤⣄⡀⠀⠀⠀⠀⠀⠀⢀⣠⣴⣾⣿⣿⠟⠁⠀⠀⠀⠀ ',
                    '⠀⠀⠀⠀⠀⠀⠈⠛⠿⣿⣿⣿⣿⣿⣶⣶⣿⣿⣿⣿⣿⠿⠋⠁⠀⠀⠀⠀⠀⠀ ',
                    '⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠉⠛⠛⠛⠛⠛⠛⠉⠉⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀ ',
                },
                footer = {},
                packages = { enable = true }, -- show how many plugins neovim loaded
                disable_move = true,
                shortcut = {
                    { desc = ' Recent', group = '@property', action = 'FzfLua oldfiles', key = 'r' },
                    { desc = '󰱼 File', group = '@property', action = 'FzfLua files', key = 'f' },
                    { desc = '󱎸 Word', group = '@property', action = 'FzfLua live_grep', key = 'w' },
                    { desc = ' Marks', group = '@property', action = 'FzfLua marks', key = 'm' },
                    { desc = ' New', group = '@property', action = 'enew', key = 'n' },
                    { desc = ' Update', group = '@property', action = 'lua vim.pack.update()', key = 'U' },
                    { desc = '󰢛 Mason', group = '@property', action = 'Mason', key = 'M' },
                },
                project = { action = 'FzfLua files cwd=' },
            }
        })
    end
}
